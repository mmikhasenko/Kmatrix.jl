"""
    Kmatrix{N,V}

K-matrix representation for a scattering system with `N` channels and `V` poles.

# Fields
- `poles::SVector{V}`: Vector of K-matrix poles, each with mass M and couplings gs
- `nonpoles::SMatrix{N,N}`: Non-pole (constant) terms in the K-matrix

The K-matrix is constructed as:
``K_{ij} = \\sum_p \\frac{g_i^p g_j^p}{m_p^2-s} + c_{ij}``
where `g_i^p` are the couplings and `c_{ij}` are the non-pole terms.
"""
struct Kmatrix{N, V}
    poles::SVector{V, NamedTuple{(:M, :gs), Tuple{Float64, SVector{N, Float64}}}}
    nonpoles::SMatrix{N, N, Float64}
end

function Kmatrix(_poles)
    V, N = length(_poles), length(first(_poles).gs)
    poles = map(_poles) do p
        (; M = p.M, gs = SVector{N}(p.gs))
    end |> SVector{V}
    nonpoles = SMatrix{N, N}(fill(0.0, (N, N)))
    return Kmatrix(poles, nonpoles)
end

amplitude(K::Kmatrix, m) =
    sum((gs * gs') ./ (M^2 - m^2) for (M, gs) in K.poles) + K.nonpoles
npoles(X::Kmatrix{N, V}) where {N, V} = V
nchannels(X::Kmatrix{N, V}) where {N, V} = N

"""
    Tmatrix{N,V}

T-matrix representation connecting K-matrix to physical scattering amplitudes.

# Fields
- `K::Kmatrix{N,V}`: Underlying K-matrix
- `channels::SVector{N,TwoBodyChannel}`: Vector of decay channels

The T-matrix is calculated as:
``T = [1-iKρ]^{-1}K``
where ρ contains the phase space factors for each channel.
"""
struct Tmatrix{N, V}
    K::Kmatrix{N, V}
    channels::SVector{N, <:AbstractChannel}
end

"""
    Dmatrix(T::Tmatrix, m)

Calculate the denominator matrix D = 1-iKρ at mass `m`.
Used in the T-matrix construction T = D⁻¹K.
"""
function Dmatrix(T::Tmatrix{N, V}, m) where {N, V}
    𝕀 = Matrix(I, (N, N))
    iρv = iρ.(T.channels, m) .* 𝕀
    K = amplitude(T.K, m)
    D = 𝕀 - K * iρv
end

detD(T::Tmatrix, m) = det(Dmatrix(T, m))
amplitude(T::Tmatrix, m) = inv(Dmatrix(T, m)) * amplitude(T.K, m)
npoles(X::Tmatrix{N, V}) where {N, V} = V
nchannels(X::Tmatrix{N, V}) where {N, V} = N
channels(X::Tmatrix) = X.channels

"""
    productionpole(T::Tmatrix, m, iR::Int)

Calculate the production amplitude contribution from the `iR`-th K-matrix pole.
Returns the term proportional to gs/(M²-s) where gs are the pole couplings.
"""
function productionpole(T::Tmatrix, m, iR::Int)
    @unpack M, gs = T.K.poles[iR]
    P = gs ./ (M^2 - m^2)
    return inv(Dmatrix(T, m)) * P
end

"""
    productionnonpole(T::Tmatrix{N,K}, m)

Calculate the production amplitude contribution from non-pole (constant) terms.
"""
function productionnonpole(T::Tmatrix{N, K}, m) where {N, K}
    return inv(Dmatrix(T, m)) * ones(N)
end
