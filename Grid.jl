#=
Grid is a cross-section of ocean along meridional line 
=#
mutable struct Grid

    Ω   :: Float64
    R   :: Float64    # Radius of planet
    H   :: Float64    # Atmosphere 

    θn :: Float64
    θs :: Float64

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
   
    θ_T   :: AbstractArray{Float64, 3}
    θ_U   :: AbstractArray{Float64, 3}
    θ_V   :: AbstractArray{Float64, 3}
    θ_W   :: AbstractArray{Float64, 3}
    θ_VW   :: AbstractArray{Float64, 3}

    z_T   :: AbstractArray{Float64, 3}
    z_U   :: AbstractArray{Float64, 3}
    z_V   :: AbstractArray{Float64, 3}
    z_W   :: AbstractArray{Float64, 3}
    z_VW   :: AbstractArray{Float64, 3}

    ∫Δv   :: Float64
    
    Δλ_W  :: AbstractArray{Float64, 3}
    
    mask_T  :: AbstractArray{Float64, 3}

    function Grid(;
        Ny    :: Int64,
        Nz    :: Int64,
        Ω     :: Float64,
        θn    :: Float64, 
        θs    :: Float64,
        H     :: Float64,
        R     :: Float64,
    )

        Δλ = [ 2π ]
        Nx = length(Δλ)

        _z_T, _z_W, _Δz_T, _Δz_W = genVerticalGrid(Nz=Nz, H=H)
        _θ_T, _θ_V, _Δθ_T, _Δθ_V = genHorizontalGrid(Nθ=Ny, θs=θs, θn=θn)

        _z_U  = copy(_z_T)
        _Δz_U = copy(_Δz_T)
        _θ_U  = copy(_θ_T)
        _Δθ_U = copy(_Δθ_T)

        #println("_Δz_W: ", size(_Δz_W))
        #println("_Δz_W: ", _Δz_W)

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


        θ_T   = y_makeMesh(_θ_T,  Nz,   Nx  )
        θ_U   = y_makeMesh(_θ_U,  Nz,   Nx+1)
        θ_V   = y_makeMesh(_θ_V,  Nz,   Nx  )
        θ_W   = y_makeMesh(_θ_T,  Nz+1, Nx  )
        θ_VW  = y_makeMesh(_θ_V,  Nz+1, Nx  )

        Δθ_T   = y_makeMesh(_Δθ_T, Nz,   Nx  )
        Δθ_U   = y_makeMesh(_Δθ_T, Nz,   Nx+1)
        Δθ_V   = y_makeMesh(_Δθ_V, Nz,   Nx  )
        Δθ_W   = y_makeMesh(_Δθ_T, Nz+1, Nx  )
        Δθ_VW  = y_makeMesh(_Δθ_V, Nz+1, Nx  )

        Δλ_T   = x_makeMesh(Δλ, Nz, Ny)
        Δλ_V   = x_makeMesh(Δλ, Nz, Ny+1)
        Δλ_W   = x_makeMesh(Δλ, Nz+1, Ny)
        Δλ_VW  = x_makeMesh(Δλ, Nz+1, Ny+1)

        # Calculating horizontal grid edges
        Δx_T = R * cos.(θ_T) .* Δλ_T;
        Δy_T = R * Δθ_T;

        Δy_U = R * Δθ_U;

        Δx_V = R * cos.(θ_V) .* Δλ_V;
        Δy_V = R * Δθ_V;
 
        Δx_W = R * cos.(θ_W) .* Δλ_W;
        Δy_W = R * Δθ_W;

        Δx_VW = R * cos.(θ_VW) .* Δλ_VW;
        Δy_VW = R * Δθ_VW;


        ∫Δv = sum(Δx_T .* Δy_T .* Δz_T)

        mask_T = Δy_T * 0.0 .+ 1

        return new(
            
            Ω,
            R,
            H,

            θn,
            θs,

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

            θ_T,
            θ_U,
            θ_V,
            θ_W,
            θ_VW,

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
    Nθ :: Int64,
    θs :: Float64,
    θn :: Float64,
)

    θ_V = collect(Float64, range(θs, θn, length=Nθ+1))
    θ_T = (θ_V[1:end-1] + θ_V[2:end]) / 2.0

    Δθ_T = similar(θ_T)
    Δθ_V = similar(θ_V)

    δθ = θ_V[2] - θ_V[1]

    Δθ_T .= δθ
    Δθ_V .= δθ
    
    return θ_T, θ_V, Δθ_T, Δθ_V

end


function genVerticalGrid(;
    Nz    :: Int64,
    H     :: Float64,
)

    #grid_shape = (η, ) -> tanh((1.0 - η) / s) / tanh(1.0/s) - 1.0
   
 
    η = collect(Float64, range(0, 1, length=Nz+1))

    #z_W  = H * grid_shape.(η)
    z_W  = H * collect(Float64, range(0, -1, length=Nz+1))

    Δz_W = similar(z_W)

    Δz_T = z_W[1:end-1] - z_W[2:end]
    Δz_W[2:end-1] = ( Δz_T[2:end] + Δz_T[1:end-1] ) / 2.0
    Δz_W[1] = Δz_W[2]
    Δz_W[end] = Δz_W[end-1]

    z_T = (z_W[1:end-1] + z_W[2:end]) / 2.0

    return z_T, z_W, Δz_T, Δz_W
end


