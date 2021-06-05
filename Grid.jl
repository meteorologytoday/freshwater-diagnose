#=
Grid is a cross-section of ocean along meridional line 
=#
mutable struct Grid

    Ω   :: Float64
    R   :: Float64    # Radius of planet
    H   :: Float64    # Atmosphere 

    Nx :: Int64
    Ny :: Int64
    Nz :: Int64
 
    Δx_T  :: AbstractArray{Float64, 3}
    Δy_T  :: AbstractArray{Float64, 3}
    Δz_T  :: AbstractArray{Float64, 3}

    Δy_U  :: AbstractArray{Float64, 3}
    Δz_U  :: AbstractArray{Float64, 3}
 
    Δx_V  :: AbstractArray{Float64, 3}
    Δy_V  :: AbstractArray{Float64, 3}
    Δz_V  :: AbstractArray{Float64, 3}
 
    Δx_W  :: AbstractArray{Float64, 3}
    Δy_W  :: AbstractArray{Float64, 3}
    Δz_W  :: AbstractArray{Float64, 3}

    Δx_VW :: AbstractArray{Float64, 3}
    Δy_VW :: AbstractArray{Float64, 3}
    Δz_VW :: AbstractArray{Float64, 3}
   
    ϕ_T   :: AbstractArray{Float64, 3}
    ϕ_U   :: AbstractArray{Float64, 3}
    ϕ_V   :: AbstractArray{Float64, 3}
    ϕ_W   :: AbstractArray{Float64, 3}
    ϕ_VW   :: AbstractArray{Float64, 3}

    z_T   :: AbstractArray{Float64, 3}
    z_U   :: AbstractArray{Float64, 3}
    z_V   :: AbstractArray{Float64, 3}
    z_W   :: AbstractArray{Float64, 3}
    z_VW   :: AbstractArray{Float64, 3}

    ∫Δv   :: Float64
    
    Δλ_W  :: AbstractArray{Float64, 3}
    
    mask_T  :: AbstractArray{Float64, 3}

    function Grid(;
        ϕ_V   :: Array{Float64, 1},
        z_W   :: Array{Float64, 1},
        Ω     :: Float64,
        R     :: Float64,
    )

        Δλ = [ 2π ]
        Nx = length(Δλ)
        Ny = length(ϕ_V) - 1
        Nz = length(z_W) - 1


        _z_W = copy(z_W)
        _z_T = (_z_W[1:end-1] + _z_W[2:end]) / 2.0
        _Δz_T = _z_W[2:end] - _z_W[1:end-1]  
        
        _Δz_W = similar(_z_W)
        _Δz_W[1:end-1] = _Δz_T
        _Δz_W[end] = _Δz_T[end]        

        _ϕ_V = copy(ϕ_V)
        _ϕ_T = (_ϕ_V[1:end-1] + _ϕ_V[2:end]) / 2.0
        _Δϕ_T = _ϕ_V[2:end] - _ϕ_V[1:end-1]  
        
        _Δϕ_V = similar(_ϕ_V)
        _Δϕ_V[1:end-1] = _Δϕ_T
        _Δϕ_V[end] = _Δϕ_T[end]


        _z_U  = copy(_z_T)
        _Δz_U = copy(_Δz_T)
        _ϕ_U  = copy(_ϕ_T)
        _Δϕ_U = copy(_Δϕ_T)

        z_makeMesh = (a, ny, nx) -> repeat( reshape(a, :, 1, 1), outer=(1, ny, nx) )
        y_makeMesh = (a, nz, nx) -> repeat( reshape(a, 1, :, 1), outer=(nz, 1, nx) )
        x_makeMesh = (a, nz, ny) -> repeat( reshape(a, 1, 1, :), outer=(nz, ny, 1) )

        z_T   = z_makeMesh(_z_T,  Ny,   Nx)
        z_U   = z_makeMesh(_z_T,  Ny,   Nx+1)
        z_V   = z_makeMesh(_z_T,  Ny+1, Nx)
        z_W   = z_makeMesh(_z_W,  Ny,   Nx)
        z_VW  = z_makeMesh(_z_W,  Ny+1, Nx)
    
        Δz_T  = z_makeMesh(_Δz_T, Ny,   Nx)
        Δz_U  = z_makeMesh(_Δz_U, Ny,   Nx+1)
        Δz_V  = z_makeMesh(_Δz_T, Ny+1, Nx)
        Δz_W  = z_makeMesh(_Δz_W, Ny,   Nx)
        Δz_VW = z_makeMesh(_Δz_W, Ny+1, Nx)


        ϕ_T   = y_makeMesh(_ϕ_T,  Nz,   Nx  )
        ϕ_U   = y_makeMesh(_ϕ_U,  Nz,   Nx+1)
        ϕ_V   = y_makeMesh(_ϕ_V,  Nz,   Nx  )
        ϕ_W   = y_makeMesh(_ϕ_T,  Nz+1, Nx  )
        ϕ_VW  = y_makeMesh(_ϕ_V,  Nz+1, Nx  )

        Δϕ_T   = y_makeMesh(_Δϕ_T, Nz,   Nx  )
        Δϕ_U   = y_makeMesh(_Δϕ_T, Nz,   Nx+1)
        Δϕ_V   = y_makeMesh(_Δϕ_V, Nz,   Nx  )
        Δϕ_W   = y_makeMesh(_Δϕ_T, Nz+1, Nx  )
        Δϕ_VW  = y_makeMesh(_Δϕ_V, Nz+1, Nx  )

        Δλ_T   = x_makeMesh(Δλ, Nz, Ny)
        Δλ_V   = x_makeMesh(Δλ, Nz, Ny+1)
        Δλ_W   = x_makeMesh(Δλ, Nz+1, Ny)
        Δλ_VW  = x_makeMesh(Δλ, Nz+1, Ny+1)

        # Calculating horizontal grid edges
        Δx_T = R * cos.(ϕ_T) .* Δλ_T;
        Δy_T = R * Δϕ_T;

        Δy_U = R * Δϕ_U;

        Δx_V = R * cos.(ϕ_V) .* Δλ_V;
        Δy_V = R * Δϕ_V;
 
        Δx_W = R * cos.(ϕ_W) .* Δλ_W;
        Δy_W = R * Δϕ_W;

        Δx_VW = R * cos.(ϕ_VW) .* Δλ_VW;
        Δy_VW = R * Δϕ_VW;


        ∫Δv = sum(Δx_T .* Δy_T .* Δz_T)

        mask_T = Δy_T * 0.0 .+ 1

        return new(
            
            Ω,
            R,
            H,

            Nx,
            Ny,
            Nz,
         
            Δx_T,
            Δy_T,
            Δz_T,

            Δy_U,
            Δz_U,

            Δx_V,
            Δy_V,
            Δz_V,
 
            Δx_W,
            Δy_W,
            Δz_W,
        
            Δx_VW,
            Δy_VW,
            Δz_VW,

            ϕ_T,
            ϕ_U,
            ϕ_V,
            ϕ_W,
            ϕ_VW,

            z_T,
            z_U,
            z_V,
            z_W,
            z_VW,

            ∫Δv,

            Δλ_W,

            mask_T,
        ) 
        
    end
end



function genHorizontalGrid(;
    Nϕ :: Int64,
    ϕs :: Float64,
    ϕn :: Float64,
)

    ϕ_V = collect(Float64, range(ϕs, ϕn, length=Nϕ+1))
    ϕ_T = (ϕ_V[1:end-1] + ϕ_V[2:end]) / 2.0

    Δϕ_T = similar(ϕ_T)
    Δϕ_V = similar(ϕ_V)

    δϕ = ϕ_V[2] - ϕ_V[1]

    Δϕ_T .= δϕ
    Δϕ_V .= δϕ
    
    return ϕ_T, ϕ_V, Δϕ_T, Δϕ_V

end


function genVerticalGrid(;
    Nz    :: Int64,
    H     :: Float64,
)

    #grid_shape = (η, ) -> tanh((1.0 - η) / s) / tanh(1.0/s) - 1.0
   
 
    η = collect(Float64, range(0, 1, length=Nz+1))

    #z_W  = H * grid_shape.(η)
    z_W  = H * collect(Float64, range(0, 1, length=Nz+1))

    Δz_W = similar(z_W)

    Δz_T = z_W[1:end-1] - z_W[2:end]
    Δz_W[2:end-1] = ( Δz_T[2:end] + Δz_T[1:end-1] ) / 2.0
    Δz_W[1] = Δz_W[2]
    Δz_W[end] = Δz_W[end-1]

    z_T = (z_W[1:end-1] + z_W[2:end]) / 2.0

    return z_T, z_W, Δz_T, Δz_W
end


