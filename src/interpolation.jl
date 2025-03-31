struct interpolatedChannel{T1, T2, Ti}
    full::T1
    inter::Ti
    approximate::T2
end

function interpolatedChannel(ch, ch_approx;
    m_grid_max = 5.5,
    n_grid = 100)
    ratio(m) = real(real_ρ(ch, m)) / imag(iρ(ch_approx, m))
    thr = threshold(ch)
    mv = range(thr, m_grid_max, n_grid)
    rv = ratio.(mv)
    inter = interpolate((mv,), rv, Gridded(Linear()))
    interpolateChannel(ch, inter, ch_approx)
end

function real_ρ(ch::interpolatedChannel, m::Real)
    @unpack inter = ch
    scale = m < inter.knots[1][end] ? inter(m) : inter.coefs[end]
    imag(iρ(ch.approximate, m)) * scale
end

# subtracted at zero 
function _iρ(ch::interpolatedChannel, m::Complex)
    thr = threshold(ch.full)
    s = m^2
    quadgk(thr^2, Inf) do s′
        m′ = sqrt(s′)
        real_ρ(ch, m′) / (s′ * (s′ - s))
    end[1] * s / π
end

# subtraction point is adjusted
function iρ(
    ch::interpolatedChannel{<:QuasiTwoBodyChannelBW}, m::Complex)
    ml = nominal_threshold(ch.full)
    f = ch.inter(ml)
    _iρ0 = real(_iρ(ch, ml + iϵ))
    _iρ(ch, m) - _iρ0
end
