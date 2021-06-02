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


A = 100.0 / 10000.0
C = (1e-2)^2

W_A_W = A * ones(Float64, mop.W_pts) |> cvtDiagOp
V_C_V = C * ones(Float64, mop.V_pts) |> cvtDiagOp

# Load data...
println("Loading data...")

#=

ddy A ddy ψ
ddz C ddz ψ
ddy B ddz ψ
ddz B ddy ψ

=#

mask_VW = reshape(collect(Integer, 1:mop.VW_pts), mop.VW_dim)
mask_VW[  1,   :] .= 0
mask_VW[end,   :] .= 0
mask_VW[  :,   1] .= 0
mask_VW[  :, end] .= 0

mask_eVW = reshape(mask_VW[mask_VW .!= 0], :)
eVW_send_VW = sparse( mop.VW_I_VW[ mask_eVW , :] )

println("eWV_send_WV constructed")

op_LHS = eVW_send_VW * (
    
    sop.VW_∂y_W * W_A_W * sop.W_DIVy_VW +
    sop.VW_∂z_V * V_C_V * sop.V_∂z_VW  
#    T_∂y_V  * V_B_V * V_∂z_VW   +
#    T_∂z_W  * W_B_W * W_DIVy_VW

) 

#op_RHS_Q =   mop.V_∂y_T
#op_RHS_F = - W_∂z_T

# Construct analytical solution
# Function we are doing is 




# Invert ψ
#@time ψ = op_LHS \ ( op_RHS_Q * Q + op_RHS_F * F  )
@time ψ = op_LHS \ RHS

# Output solution


