using NCDatasets

include("streamfunction.jl")
include("Grid.jl")
include("MatrixOperators.jl")
include("MatrixSpatialOperators.jl")



p0 = 1000.0 # hPa
θ0 = 300.0
cp = 1004.0
H  = cp * θ0 / g


calΠ = (p,) -> ( p / p0 ).^(2/5) 
cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))
conform = (a,) -> transpose(a)[end:-1:1, :]

println("Loading data...")


Dataset("avgzm.nc", "r") do ds

    global T_V    = ds["T"][:]     |> conform
    global U_V    = ds["U"][:]     |> conform
    global V_V    = ds["V"][:]     |> conform
    global lev  = ds["lev"][:]  
    global ilev = ds["ilev"][:]
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

@time sop = MatrixSpatialOperators(gd)
mop = sop.op

reshapeVW = (a,) -> reshape(a, mop.VW_dim) 
reshapeW = (a,) -> reshape(a, mop.W_dim) 
reshapeT = (a,) -> reshape(a, mop.T_dim) 

# Stream function
ψ_VW = reshapeVW(calψ(ilev=ilev, lat=lat, v=convert(Array{Float64}, V_V)))

# Derive structure ABC
println("Deriving A, B, and C...")

# Convert grid
T_T = reshapeT(sop.T_interp_V * T_V[:])

# Calculate potential temperature
p0 = 1000.0 # hPa
Π = calΠ.(repeat(reshape(lev, :, 1), outer=(1, gd.Ny)))
θ_T = T_T .* Π

A_W = reshapeW( g/θ0 * sop.W_∂z_T * θ_T[:]) 
W_A_W = cvtDiagOp(A_W)

#=
A = 100.0 / 10000.0
C = (1e-2)^2
B = A * 0.1

V_C_V = C * ones(Float64, mop.V_pts) |> cvtDiagOp
V_B_V = B * ones(Float64, mop.V_pts) |> cvtDiagOp
W_B_W = B * ones(Float64, mop.W_pts) |> cvtDiagOp
=#

#=
ϕ = gd.ϕ_T
Z = gd.z_T / gd.H
ϕc = deg2rad(45)
σϕ = deg2rad(5)
Q0 = 0.1 # K / sec
Q = Q0 * exp.( - ((ϕ .- ϕc) / σϕ ).^2 ) .* sin.(π*Z)
=#


#=

ddy A ddy ψ
ddz C ddz ψ
ddy B ddz ψ
ddz B ddy ψ

=#

#=
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



Eψ_n = reshapeVW( VW_E_VW * ψ[:] )

∂ψ∂y = reshape( sop.W_∂y_VW * ψ[:], mop.W_dim )
DIVyψ = reshape( sop.W_DIVy_VW * ψ[:], mop.W_dim )

#op_RHS_Q =   mop.V_∂y_T
#op_RHS_F = - W_∂z_T



# Invert ψ
#@time ψ_n = evW_E_VW \ ( op_RHS_Q * Q + op_RHS_F * F  )
println("Inverting...")
@time ψ_n = VW_send_eVW * ( eVW_E_eVW \ ( eVW_send_VW * Eψ_n[:] ) ) |> reshapeVW



RHS_Q = VW_send_eVW * eVW_Qop_T * Q[:] |> reshapeVW
@time ψ1 = VW_send_eVW * ( eVW_E_eVW \ (eVW_send_VW * RHS_Q[:] ) ) |> reshapeVW
# Output solution
=#

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

    z_T[:] = gd.z_T[:, 1, 1] 
    z_W[:] = gd.z_W[:, 1, 1]
    y_T[:] = gd.ϕ_T[1, :, 1]
    y_V[:] = gd.ϕ_V[1, :, 1]

    # simulation variables            

    defVar(ds, "psi",  ψ_VW,    ("Nzp1", "Nyp1", "Nx"))
    defVar(ds, "A_W",  A_W,     ("Nzp1", "Ny",   "Nx"))
    defVar(ds, "T_T",  T_T,     ("Nz",   "Ny",   "Nx"))
    

end
