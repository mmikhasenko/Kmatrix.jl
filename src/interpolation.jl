struct InterpolatedChannel{T1, T2, Ti} <: AbstractChannel
    full::T1
    inter::Ti
    approximate::T2
end

function InterpolatedChannel(ch, ch_approx;
    m_grid_max = 5.5,
    n_grid = 100)
    ratio(m) = real(real_ρ(ch, m)) / imag(iρ(ch_approx, m))
    thr = threshold(ch)
    mv = range(thr, m_grid_max, n_grid)
    rv = ratio.(mv)
    inter = interpolate((mv,), rv, Gridded(Linear()))
    InterpolatedChannel(ch, inter, ch_approx)
end

function real_ρ(ch::InterpolatedChannel, m::Real)
    @unpack inter = ch
    scale = m < inter.knots[1][end] ? inter(m) : inter.coefs[end]
    imag(iρ(ch.approximate, m)) * scale
end

function once_subtracted_dispersion_relation(f_of_s, s::Complex, s_thr::Real)
    quadgk(s_thr, Inf) do s′
        f_of_s(s′) / (s′ * (s′ - s))
    end[1] * s / π
end

# subtracted at zero 
function _iρ(ch::InterpolatedChannel, m::Complex)
    m_thr = threshold(ch.full)
    N(s) = real_ρ(ch, sqrt(s))
    return once_subtracted_dispersion_relation(N, m^2, m_thr^2)
end

# subtraction point is adjusted
function iρ(
    ch::InterpolatedChannel{<:QuasiTwoBodyChannelBW}, m::Complex)
    ml = nominal_threshold(ch.full)
    _iρ0 = real(_iρ(ch, ml + 1e-7im))
    return _iρ(ch, m) - _iρ0
end
