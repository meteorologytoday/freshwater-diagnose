using NCDatasets
using Formatting

include("constants.jl")

function calψ(;
    ilev  :: Array{Float64, 1},
    lat   :: Array{Float64, 1},
    v     :: Array{Float64, 2},
)

    # Here we assume that ilev is monotonically decreaseing (from lower to upper atmosphere)
    # We also assume that dimension is (z, y)
    Nz, Ny = size(v)
    Nzp1 = length(ilev)
    dp = ilev[2:end] - ilev[1:end-1]
    
    #println("Nz, Ny = ", size(v)) 
    ψ = zeros(Float64, Nzp1, Ny)  
    ψ[1, :] .= 0.0  # Set the lowest layer to zero
    for j=1:Ny
        for k=2:Nzp1
            ψ[k, j] = ψ[k-1, j] + dp[k-1] * v[k-1, j]
        end

        ψ[:, j] .*= 2π * Re * cos(deg2rad(lat[j])) / g
    end

    return ψ
end
