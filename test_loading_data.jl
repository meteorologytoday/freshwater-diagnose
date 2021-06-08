using NCDatasets

include("streamfunction.jl")
include("Grid.jl")
include("MatrixOperators.jl")
include("MatrixSpatialOperators.jl")


R  = 287.0
cp = 5/2 * R
p0 = 1e5 # Pa
θ0 = 300.0
H = cp * θ0 / g
ρ0 = p0 / (R * θ0)
κ = R / cp

println("θ0 = ", θ0)

calσ = (z,) -> ρ0 * ( 1.0 - z / H)^(1/κ - 1)
calΠ = (p,) -> ( p / p0 )^(2/5) 
cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))
conform = (a,) -> transpose(a)[end:-1:1, :]

println("Loading data...")


Dataset("avgzm.nc", "r") do ds

    global T_V    = ds["T"][:]     |> conform
    global U_V    = ds["U"][:]     |> conform
    global V_V    = ds["V"][:]     |> conform
    global VU_V   = ds["VU"][:]     |> conform
    global lev  = ds["lev"][:]  * 100.0
    global ilev = ds["ilev"][:] * 100.0
    global lat  = ds["lat"][:]
    global slat = ds["slat"][:]

    lev = lev[end:-1:1]
    ilev = ilev[end:-1:1]
end

z_W = H * ( 1.0 .- calΠ.(ilev))
ϕ_V = deg2rad.(lat)

print("Making grid...")
gd = Grid(;
    ϕ_V = ϕ_V,
    z_W = z_W,
    Ω  = 7.292e-5,
    R  = 6.4e6,
) 
println("done.")

@time sop = MatrixSpatialOperators(gd, ocn_or_atm=:atm)
mop = sop.op

reshapeVW = (a,) -> reshape(a, mop.VW_dim) 
reshapeW = (a,) -> reshape(a, mop.W_dim) 
reshapeV = (a,) -> reshape(a, mop.V_dim) 
reshapeT = (a,) -> reshape(a, mop.T_dim) 

# Stream function
ψ_VW = reshapeVW(calψ(ilev=ilev, lat=lat, v=convert(Array{Float64}, V_V)))

# Derive structure ABC
println("Deriving A, B, and C...")

# Convert grid
T_T = reshapeT(sop.T_interp_V * T_V[:])
U_T = reshapeT(sop.T_interp_V * U_V[:])

# Calculate potential temperature
p0 = 1000.0 # hPa
Π_T = calΠ.(repeat(reshape(lev, :, 1), outer=(1, gd.Ny)))
Π_V = calΠ.(repeat(reshape(lev, :, 1), outer=(1, gd.Ny+1)))
θ_T = T_T ./ Π_T
θ_V = T_V ./ Π_V

σ_W = calσ.(gd.z_W)
σ_V = calσ.(gd.z_V)
σ_T = calσ.(gd.z_T)

W_invσ_W = σ_W.^(-1) |> cvtDiagOp
V_invσ_V = σ_V.^(-1) |> cvtDiagOp
T_invσ_T = σ_T.^(-1) |> cvtDiagOp

A_W = reshapeW( g/θ0 * W_invσ_W * sop.W_∂z_T * θ_T[:] ) 
W_A_W = cvtDiagOp(A_W)

B_V = reshapeV( - g/θ0 * V_invσ_V * sop.V_∂y_T * θ_T[:]) 
B_W = reshapeW( - g/θ0 * W_invσ_W * sop.W_interp_T * sop.T_∂y_V * θ_V[:]) 

V_B_V = cvtDiagOp(B_V)
W_B_W = cvtDiagOp(B_W)

# Calculate vorticity
f_T = 2 * gd.Ω * sin.(gd.ϕ_T)
ζ_T = - reshapeT(sop.T_DIVy_V * U_V[:])
η_T = f_T + ζ_T

C_T = η_T .* ( f_T + 2 * U_T .* tan.(gd.ϕ_T) / gd.R) ./ σ_T
C_V = reshapeV( sop.V_interp_T * C_T[:] )
V_C_V = cvtDiagOp(C_V)


A_T = reshapeT( sop.T_interp_W * A_W[:] )
B_T = reshapeT( sop.T_interp_V * B_V[:] )

Δ_T = A_T .* C_T - B_T.^2

#=
A = 100.0 / 10000.0
C = (1e-2)^2
B = A * 0.1

V_C_V = C * ones(Float64, mop.V_pts) |> cvtDiagOp
V_B_V = B * ones(Float64, mop.V_pts) |> cvtDiagOp
W_B_W = B * ones(Float64, mop.W_pts) |> cvtDiagOp
=#

ϕ = gd.ϕ_T
Z = gd.z_T / gd.H
ϕc = deg2rad(45)
σϕ = deg2rad(2)
Q0 = 0.1 # K / sec
Q = Q0 * exp.( - ((ϕ .- ϕc) / σϕ ).^2 ) .* sin.(π*Z)


#=

ddy A ddy ψ
ddz C ddz ψ
ddy B ddz ψ
ddz B ddy ψ

=#

mask_VW = reshapeVW(collect(Integer, 1:mop.VW_pts))
mask_VW[  1,   :] .= 0
mask_VW[end,   :] .= 0
mask_VW[  :,   1] .= 0
mask_VW[  :, end] .= 0

mask_eVW = reshape(mask_VW[mask_VW .!= 0], :)
eVW_send_VW = sparse( mop.VW_I_VW[ mask_eVW , :] )
VW_send_eVW = transpose(eVW_send_VW)
println("eWV_send_WV constructed")

# Making Operators
VW_E_VW = (
    sop.VW_∂y_W * W_A_W * sop.W_DIVy_VW
    + sop.VW_∂z_V * V_C_V * sop.V_∂z_VW  
    + sop.VW_interp_T * sop.T_∂y_V  * V_B_V * sop.V_∂z_VW 
    + sop.VW_interp_T * sop.T_∂z_W  * W_B_W * sop.W_DIVy_VW
) 

eVW_E_eVW = eVW_send_VW * VW_E_VW * VW_send_eVW
eVW_Qop_T = eVW_send_VW * sop.VW_interp_V * sop.V_∂y_T

F_T = reshapeT( (
    - sop.T_DIVy_V
    + sop.T_interp_V * cvtDiagOp(tan.(gd.ϕ_V)) / gd.R
) * VU_V[:] )

eVW_Fop_T =  eVW_send_VW * sop.VW_interp_W * sop.W_∂z_T * cvtDiagOp(f_T + 2 * U_T .* tan.(gd.ϕ_T) / gd.R)


# Invert ψ
RHS_Q = VW_send_eVW * eVW_Qop_T * Q[:] |> reshapeVW
@time ψ1 = VW_send_eVW * ( eVW_E_eVW \ (eVW_send_VW * RHS_Q[:] ) ) |> reshapeVW

RHS_F = VW_send_eVW * eVW_Fop_T * F_T[:] |> reshapeVW
@time ψ2 = VW_send_eVW * ( eVW_E_eVW \ (eVW_send_VW * RHS_F[:] ) ) |> reshapeVW

Ψ2 = ψ2 .* ( 2π * gd.R * cos.(gd.ϕ_VW))


# Output solution

println("Outputting result...")
Dataset("output.nc", "c") do ds

    defDim(ds, "Nz", gd.Nz)
    defDim(ds, "Ny", gd.Ny)
    defDim(ds, "Nx", gd.Nx)
    
    defDim(ds, "Nzp1", gd.Nz+1)
    defDim(ds, "Nyp1", gd.Ny+1)
    defDim(ds, "Nxp1", gd.Nx+1)
    
    defDim(ds, "time", Inf)
    # coordinate variables 
    z_T    = defVar(ds, "z_T", Float64, ("Nz",))
    z_W    = defVar(ds, "z_W", Float64, ("Nzp1",))
    y_T    = defVar(ds, "y_T", Float64, ("Ny",))
    y_V    = defVar(ds, "y_V", Float64, ("Nyp1",))
    p_T    = defVar(ds, "p_T", Float64, ("Nz",))
    p_W    = defVar(ds, "p_W", Float64, ("Nzp1",))

    p_T[:] = lev
    p_W[:] = ilev
    z_T[:] = gd.z_T[:, 1, 1] 
    z_W[:] = gd.z_W[:, 1, 1]
    y_T[:] = gd.ϕ_T[1, :, 1] .|> rad2deg
    y_V[:] = gd.ϕ_V[1, :, 1] .|> rad2deg

    # simulation variables            

    defVar(ds, "psi",  ψ_VW,    ("Nzp1", "Nyp1", "Nx"))
    defVar(ds, "A_W",  A_W,     ("Nzp1", "Ny",   "Nx"))
    defVar(ds, "T_T",  T_T,     ("Nz",   "Ny",   "Nx"))
    defVar(ds, "theta_T",  θ_T, ("Nz",   "Ny",   "Nx"))
    defVar(ds, "B_W",  B_W,     ("Nzp1", "Ny",   "Nx"))
    defVar(ds, "B_V",  B_V,     ("Nz",   "Nyp1", "Nx"))
    defVar(ds, "C_V",  C_V,     ("Nz",   "Nyp1", "Nx"))
    defVar(ds, "zeta_T",  ζ_T,  ("Nz",   "Ny", "Nx"))
    defVar(ds, "f_T",  f_T,     ("Nz",   "Ny", "Nx"))
    defVar(ds, "Delta_T",  Δ_T, ("Nz",   "Ny", "Nx"))
    
    defVar(ds, "psi1",  ψ1,    ("Nzp1", "Nyp1", "Nx"))
    defVar(ds, "Q1",    Q,     ("Nz",   "Ny",   "Nx"))
    
    defVar(ds, "psi2",  ψ2,    ("Nzp1", "Nyp1", "Nx"))
    defVar(ds, "PSI2",  Ψ2,    ("Nzp1", "Nyp1", "Nx"))
    defVar(ds, "F2",    VU_V,  ("Nz",   "Nyp1"))
    defVar(ds, "RHS_F2",RHS_F, ("Nzp1", "Nyp1", "Nx"))

end
