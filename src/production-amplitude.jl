"""
    ProductionAmplitude{N,V}

Production amplitude representation for N channels and V poles.

# Fields
- `T::TMatrix{N,V}`: Underlying T-matrix description
- `α_poles::SVector{V}`: Production couplings to each K-matrix pole
- `α_nonpoles::SVector{N}`: Direct production couplings to each channel

The production amplitude is calculated as:
```math
A = [1-iKρ]^{-1}P
```
where P contains both pole and non-pole production terms.
"""
struct ProductionAmplitude{N, V}
    T::TMatrix{N, V}
    α_poles::SVector{V, <:Number}
    α_nonpoles::SVector{N, <:Number}
end

npoles(X::ProductionAmplitude{N, V}) where {N, V} = V
nchannels(X::ProductionAmplitude{N, V}) where {N, V} = N
detD(X::ProductionAmplitude, m) = detD(X.T, m)
channels(X::ProductionAmplitude) = channels(X.T)

ProductionAmplitude(T::TMatrix{N, V}, α_poles, α_nonpoles = zeros(N)) where {N, V} =
    ProductionAmplitude(T, SVector{V}(α_poles), SVector{N}(α_nonpoles))

"""
    amplitude(A::ProductionAmplitude, m)

Calculate the full production amplitude at mass `m`.
Includes contributions from both pole and non-pole terms.
"""
function amplitude(A::ProductionAmplitude, m)
    T, α_poles, α_nonpoles = A.T, A.α_poles, A.α_nonpoles
    P = α_nonpoles
    for (α, Mgs) in zip(α_poles, A.T.K.poles)
        M, gs = Mgs.M, Mgs.gs
        P += α .* gs ./ (M^2 - m^2)
    end
    D⁻¹ = inv(DMatrix(T, m))
    return D⁻¹ * P
end

"""
    production_pole(A::ProductionAmplitude, m, iR::Int)

Calculate production amplitude contribution from only the `iR`-th pole.
Useful for studying individual resonance contributions.
"""
function production_pole(A::ProductionAmplitude{N, V}, m, iR::Int) where {N, V}
    α_nonpoles = SVector{N}(zeros(N))
    α_poles = zeros(Complex{Float64}, V)
    α_poles[iR] = A.α_poles[iR]
    A = ProductionAmplitude(A.T, SVector{V}(α_poles), α_nonpoles)
    return amplitude(A, m)
end

"""
    production_nonpole(A::ProductionAmplitude, m)

Calculate production amplitude contribution from only non-pole terms.
"""
function production_nonpole(A::ProductionAmplitude{N, V}, m) where {N, V}
    @unpack T, α_nonpoles = A
    α_poles = SVector{V}(zeros(V))
    A = ProductionAmplitude(T, α_poles, α_nonpoles)
    return amplitude(A, m)
end