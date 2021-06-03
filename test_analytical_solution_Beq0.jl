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

W_A_W = A * ones(Float64, mop.W_pts) |> cvtDiagOp
V_C_V = C * ones(Float64, mop.V_pts) |> cvtDiagOp

target_solution = 3

println("##### Target solution: ", target_solution, " #####")

if target_solution == 1

    # === analytical solution set 1: z direction ===
    #
    # Construct analytical solution
    # Function we are doing is 
    #
    # ψ = sin(πk z/H)
    # 
    # After it is operated, the solution is
    #
    # Eψ = - C (πk/H)^2 ψ - A / a^2 * ψ/cos^2(ϕ)
    #

 
    k = 2.0
    ψ = sin.(π * k * gd.z_VW / gd.H)

    Eψ_a = - C * (π * k / gd.H)^2 * ψ - A / gd.R^2 * ψ ./ (cos.(gd.θ_VW)).^2 

elseif target_solution == 2

    # === analytical solution set 1: y direction ===
    Δθ = gd.θn - gd.θs
    m = 2.0
    ψ = sin.(π * m * (gd.θ_VW .- gd.θs) / Δθ)
    
    Eψ_a = (
            - A / gd.R^2 * (
                ψ ./ (cos.(gd.θ_VW)).^2.0
             .+ tan.(gd.θ_VW) * (π * m / Δθ) .* cos.(π * m * (gd.θ_VW .- gd.θs) / Δθ)
             .+ (π*m/Δθ)^2 * ψ
            )
    )
elseif target_solution == 3

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
    # = - ψ ( A / a^2 * ( 2 + 1/cos^2(ϕ) ) + C π^2/H^2 )
    #

    Δθ = gd.θn - gd.θs
    m = 2.0
    k = 2.0

    ψ = sin.(π * m * (gd.θ_VW .- gd.θs) / Δθ) .* sin.(π * k * gd.z_VW / gd.H)

    #Eψ_a = - ψ .* ( A / gd.R^2 * (2 .+ 1 ./ ( cos.(gd.θ_VW) ).^2 .+ C * π^2 / gd.H^2) )

    Eψ_a = ( - C * (π * k / gd.H)^2 * ψ 
            - A / gd.R^2 * (
                ψ ./ (cos.(gd.θ_VW)).^2.0
             .+ tan.(gd.θ_VW) * (π * m / Δθ) .* cos.(π * m * (gd.θ_VW .- gd.θs) / Δθ) .* sin.(π * k * gd.z_VW / gd.H)
             .+ (π*m/Δθ)^2 * ψ
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
#    T_∂y_V  * V_B_V * V_∂z_VW   +
#    T_∂z_W  * W_B_W * W_DIVy_VW

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