"""
    ProductionAmplitude{N,V}

Production amplitude representation for N channels and V poles.

# Fields
- `T::Tmatrix{N,V}`: Underlying T-matrix description
- `αpoles::SVector{V}`: Production couplings to each K-matrix pole
- `αnonpoles::SVector{N}`: Direct production couplings to each channel

The production amplitude is calculated as:
```math
A = [1-iKρ]^{-1}P
```
where P contains both pole and non-pole production terms.
"""
struct ProductionAmplitude{N, V}
    T::Tmatrix{N, V}
    αpoles::SVector{V, <:Number}
    αnonpoles::SVector{N, <:Number}
end

npoles(X::ProductionAmplitude{N, V}) where {N, V} = V
nchannels(X::ProductionAmplitude{N, V}) where {N, V} = N
detD(X::ProductionAmplitude, m; ϕ = -π / 2) = detD(X.T, m; ϕ)
channels(X::ProductionAmplitude) = channels(X.T)

ProductionAmplitude(T::Tmatrix{N, V}) where {N, V} =
    ProductionAmplitude(T, SVector{V}(ones(V)), SVector{N}(ones(N)))

"""
    amplitude(A::ProductionAmplitude, m; ϕ=-π/2)

Calculate the full production amplitude at mass `m`.
Includes contributions from both pole and non-pole terms.
"""
function amplitude(A::ProductionAmplitude, m; ϕ = -π / 2)
    T, αpoles, αnonpoles = A.T, A.αpoles, A.αnonpoles
    P = αnonpoles
    for (α, Mgs) in zip(αpoles, A.T.K.poles)
        M, gs = Mgs.M, Mgs.gs
        P += α .* gs ./ (M^2 - m^2)
    end
    D⁻¹ = inv(Dmatrix(T, m; ϕ))
    return D⁻¹ * P
end

"""
    production_pole(A::ProductionAmplitude, m, iR::Int; ϕ=-π/2)

Calculate production amplitude contribution from only the `iR`-th pole.
Useful for studying individual resonance contributions.
"""
function production_pole(A::ProductionAmplitude{N, V}, m, iR::Int; ϕ = -π / 2) where {N, V}
    αnonpoles = SVector{N}(zeros(N))
    αpoles = zeros(Complex{Float64}, V)
    αpoles[iR] = A.αpoles[iR]
    A = ProductionAmplitude(A.T, SVector{V}(αpoles), αnonpoles)
    return amplitude(A, m; ϕ)
end

"""
    production_nonpole(A::ProductionAmplitude, m; ϕ=-π/2)

Calculate production amplitude contribution from only non-pole terms.
"""
function production_nonpole(A::ProductionAmplitude{N, V}, m; ϕ = -π / 2) where {N, V}
    @unpack T, αnonpoles = A
    αpoles = SVector{V}(zeros(V))
    A = ProductionAmplitude(T, αpoles, αnonpoles)
    return amplitude(A, m; ϕ)
end