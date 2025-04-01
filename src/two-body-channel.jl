abstract type AbstractChannel end

"""
    TwoBodyChannel(m1, m2; L=0)

Two-body channel representation with masses `m1`, `m2` and angular momentum `L`.

# Fields
- `m1::Complex{Float64}`: Mass of first particle
- `m2::Complex{Float64}`: Mass of second particle  
- `L::Int`: Angular momentum quantum number (currently only L=0 implemented)
"""
struct TwoBodyChannel <: AbstractChannel
    m1::Complex{Float64}
    m2::Complex{Float64}
    L::Int
end

TwoBodyChannel(m1, m2; L::Int = 0) = TwoBodyChannel(m1, m2, L)

threshold(ch::TwoBodyChannel) = real(ch.m1 + ch.m2)

"""
    ρ(ch::TwoBodyChannel, m)

Calculate the phase space factor for a two-body channel at mass `m`.

# Arguments
- `ch::TwoBodyChannel`: Two-body channel definition
- `m::Number`: Mass where to evaluate phase space

# Returns
Phase space value ρ(m) for the given channel and mass
"""
function iρ(ch::TwoBodyChannel, m)
    ϕ = -π / 2
    ch.L != 0 && error("not implemented")
    1im *
    sqrt(cis(ϕ) * (m - (ch.m1 + ch.m2))) * cis(-ϕ / 2) *
    sqrt(m + (ch.m1 + ch.m2)) *
    sqrt((m^2 - (ch.m1 - ch.m2)^2)) /
    m^2
end
