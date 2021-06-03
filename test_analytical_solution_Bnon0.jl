using NCDatasets

include("Grid.jl")
include("MatrixOperators.jl")
include("MatrixSpatialOperators.jl")


cvtDiagOp = (a,) -> spdiagm(0 => view(a, :))

print("Making grid...")
gd = Grid(;
    Ny = 600,
    Nz = 300,
    Ω  = 7.292e-5,
    θn = deg2rad(70),
    θs = deg2rad(20),
    H  = 1.5e4,
    R  = 6.4e6,
) 
println("done.")

@time sop = MatrixSpatialOperators(gd)
mop = sop.op

reshapeVW = (a,) -> reshape(a, mop.VW_dim) 

A = 100.0 / 10000.0
C = (1e-2)^2

B = A / 2

W_A_W = A * ones(Float64, mop.W_pts) |> cvtDiagOp
V_C_V = C * ones(Float64, mop.V_pts) |> cvtDiagOp
V_B_V = B * ones(Float64, mop.V_pts) |> cvtDiagOp
W_B_W = B * ones(Float64, mop.W_pts) |> cvtDiagOp

target_solution = 1

println("##### Target solution: ", target_solution, " #####")

if target_solution == 1

    # === analytical solution set 3 : both y and z ===
    #
    # Construct analytical solution
    # Function we are doing is 
    #
    # ψ = sin(π m (ϕ - ϕs) / Δϕ) sin(πk z/H)
    # 
    # where Δϕ = ϕn - ϕs
    #
    # After it is operated, the solution is
    #
    # Eψ = - C (πk/H)^2 ψ 
    #      - A / a^2 * (
    #           ψ/cos^2(ϕ)
    #         + tan(ϕ) (πm / Δϕ) cos(πm (ϕ-ϕs) / Δϕ) sin(πk z/H)
    #         - (πm/Δϕ)^2 ψ
    #       )
    #   - A ψ / a^2 * ( 2 + 1/cos^2(ϕ) ) - C π^2/H^2 * ψ
    #   + B / a * (
    #        (πm/Δϕ) (πk/H) cos(πm (ϕ-ϕs) / Δϕ) cos(πk z/H)
    #      - sin(ϕ) (πk/H) sin(πm (ϕ-ϕs) / Δϕ) cos(πk z/H) 
    #   )
    #
    #

    Δθ = gd.θn - gd.θs
    m = 2.0
    k = 2.0

    θ = gd.θ_VW
    M = m*π/Δθ
    K = k*π/gd.H
    Y = M * (θ .- gd.θs)
    Z = K * gd.z_VW


    ψ = sin.(Y) .* sin.(Z)

    Eψ_a = ( - C * K^2 * ψ 
            - A / gd.R^2 * (
                ψ ./ (cos.(θ)).^2.0
             .+ tan.(θ) * M .* cos.(Y) .* sin.(Z)
             .+ M^2 * ψ
            )
            + B / gd.R * (
            2 * M * K * cos.(Y) .* cos.(Z)
          - K * sin.(θ) .* sin.(Y) .* cos.(Z) 
          )
    )

else

    throw(ErrorException("Wrong solution target..."))
end

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
    
    sop.VW_∂y_W * W_A_W * sop.W_DIVy_VW
    + sop.VW_∂z_V * V_C_V * sop.V_∂z_VW  
    + sop.VW_interp_T * sop.T_∂y_V  * V_B_V * sop.V_∂z_VW 
    + sop.VW_interp_T * sop.T_∂z_W  * W_B_W * sop.W_DIVy_VW

) 

Eψ_n = reshapeVW( VW_send_eVW * op_LHS * ψ[:] )

∂ψ∂y = reshape( sop.W_∂y_VW * ψ[:], mop.W_dim )
DIVyψ = reshape( sop.W_DIVy_VW * ψ[:], mop.W_dim )

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
    defVar(ds, "Error_abs", Eψ_n - Eψ_a, ("Nzp1", "Nyp1", "Nx"))
    defVar(ds, "Error_rto", (Eψ_n - Eψ_a) ./ Eψ_a, ("Nzp1", "Nyp1", "Nx"))

    defVar(ds, "dpsidy",  ∂ψ∂y,    ("Nzp1", "Ny", "Nx"))
    defVar(ds, "DIVypsi",  DIVyψ,  ("Nzp1", "Ny", "Nx"))

    

end
