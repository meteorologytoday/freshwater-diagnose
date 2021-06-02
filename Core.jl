mutable struct Core

    sop_b         :: MatrixSpatialOperators
    sop_bb        :: MatrixSpatialOperators
    sop_bib       :: MatrixSpatialOperators

    es            :: EllipticSolver
    vds           :: VerticalDiffusionSolver
    hds           :: HorizontalDiffusionSolver
    wksp          :: Workspace

    op_mtx        :: Any         
    lu_mtx        :: Any 

    function Core(env)
        
        Nx = 2
        Ny = env.Ny_bsn
        Nz = env.Nz

        println("Making matricies")
        @time sop_b   = MatrixSpatialOperators(env.gd_b)
        @time sop_bb  = MatrixSpatialOperators(env.gd_bb)

        println("Making bib")
        @time sop_bib = MatrixSpatialOperators(env.gd_bib)


        es = EllipticSolver(env.gd_b, sop_b, sop_bb)
        vds = VerticalDiffusionSolver(sop_bb; K_iso=env.Kv, K_cva=env.Kv_cva)
        hds = HorizontalDiffusionSolver(sop_bib; K=env.Kh)
        
        wksp = Workspace(;
            Ny = Ny,
            Nz = Nz,
        ) 

        hdiv     = zeros(Float64, 1, Ny)

        forcing_coe =  - ones(Float64, Nz, Ny, Nx) / env.Δt_forcing
        forcing_coe[2:end, :, :] .= 0  # Only restore the surface
        
        # Some additional matrices for Newton method 

        # calculate projection matrix T-grid to V- and W-grid

        # TODO: make your own U_interp_T
        P = [
                sop_bb.U_interp_T ;
                sop_bb.V_interp_T ;
                sop_bb.W_interp_T
        ]

        ∇ = sparse([ sop_bb.T_DIVx_U sop_bb.T_DIVy_V  sop_bb.T_DIVz_W ])
        
        # calculate the dilute matrix which only
        # operates on the eastern (i.e. second slab)
        D = ones(Float64, sop_bb.op.T_pts)
        D[(1+Nz*Ny):(Nz*Ny*2)] .= env.Δλb / (env.Δλb + env.Δλi)  
        D = spdiagm(0 => D)
        
        D∇ = sparse(D * ∇)
        
        S, G, Q = calSGQ(es)


        #### Old formulation #####
        # project from U,V,W grid onto T
        P_old = blockdiag(spdiagm(sop_bb.op.T_pts, sop_bb.op.U_pts), sop_bb.T_interp_V, sop_bb.T_interp_W)
        Σ_old = spdiagm(0 => ones(Float64, sop_bb.op.T_pts))
        Σ_old = [ Σ_old Σ_old Σ_old ]

        ΣP_old = Σ_old * P_old
        ∇_old = sparse([ spdiagm(sop_bb.op.U_pts, sop_bb.op.T_pts) ; sop_bb.V_∂y_T ; sop_bb.W_∂z_T ])



        op_mtx = Dict( 
#                :zonal_diffusion => calOp_zonal_diff(sopz, env.Khx),
#                :∇       => ∇,
                :P       => P,
                :D∇      => D∇,
                :S       => S, 
                :G       => G, 
                :Q       => Q, 
                :hdiff   => calOp_hdiff(hds),
                :hdiff_zonal => calOp_zonal_diff(sop_bib, env.Khx, env.gd_bib.Δx_T[:, :, 1]),
                :forcing => spdiagm( 0 => forcing_coe[:] ),
                :∇_old   => ∇_old,
                :ΣP_old  => ΣP_old,
        )

        lu_mtx = Dict(
#            :EBM_zonal_diffusion => lu( sopz.op.T_I_T - env.Δt * op_mtx[:3][:zonal_diffusion] ),
        )

        new(
            sop_b,
            sop_bb,
            sop_bib,
#            sopz,
            es,
            vds,
            hds,
            wksp,
            op_mtx,
            lu_mtx,
        )
        
    end
end


