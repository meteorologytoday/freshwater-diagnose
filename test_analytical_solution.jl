using NCDatasets

include("Grid.jl")
include("MatrixOperators.jl")
include("MatrixSpatialOperators.jl")


cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))

print("Making grid...")
gd = Grid(;
    Ny = 60,
    Nz = 30,
    Ω  = 7.292e-5,
    θn = 80.0,
    θs = 10.0,
    H  = 1.5e4,
    R  = 6.4e6,
) 
println("done.")

@time sop = MatrixSpatialOperators(gd)
mop = sop.op

reshapeVW = (a,) -> reshape(a, mop.VW_dim) 

A = 100.0 / 10000.0
C = (1e-2)^2

W_A_W = A * ones(Float64, mop.W_pts) |> cvtDiagOp
V_C_V = C * ones(Float64, mop.V_pts) |> cvtDiagOp

#
# Construct analytical solution
# Function we are doing is 
#
# ψ = sin(ϕ) sin(πz/H)
# 
# After it is operated, the solution is
#
#   - A ψ / a^2 * ( 2 + 1/cos^2(ϕ) ) - C π^2/H^2 * ψ
# = - ψ ( A / a^2 * ( 2 + 1/cos^2(ϕ) ) + C π^2/H^2 )
#

θ_VW = deg2rad.(gd.θ_VW)

ψ = sin.(θ_VW) .* sin.(π * gd.z_VW ./ gd.H)
Eψ_a = - ψ .* ( A / gd.R^2 * (2 .+ 1 ./ ( cos.(θ_VW) ).^2 .+ C * π^2 / gd.H^2) )


# Load data...
println("Loading data...")

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

op_LHS = eVW_send_VW * (
    
    sop.VW_∂y_W * W_A_W * sop.W_DIVy_VW +
    sop.VW_∂z_V * V_C_V * sop.V_∂z_VW  
#    T_∂y_V  * V_B_V * V_∂z_VW   +
#    T_∂z_W  * W_B_W * W_DIVy_VW

) 

Eψ_n = reshapeVW( VW_send_eVW * op_LHS * ψ[:] )

#op_RHS_Q =   mop.V_∂y_T
#op_RHS_F = - W_∂z_T



# Invert ψ
#@time ψ = op_LHS \ ( op_RHS_Q * Q + op_RHS_F * F  )
# @time ψ = op_LHS \ RHS

# Output solution


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
    y_T[:] = gd.θ_T[1, :, 1]
    y_V[:] = gd.θ_V[1, :, 1]

    # simulation variables            

    defVar(ds, "psi",    ψ,    ("Nzp1", "Nyp1", "Nx"))
    defVar(ds, "Epsi_a", Eψ_a, ("Nzp1", "Nyp1", "Nx"))
    defVar(ds, "Epsi_n", Eψ_n, ("Nzp1", "Nyp1", "Nx"))

    

end
