@with_kw struct QuasiTwoBodyChannel{T}
    m1::Complex{Float64}
    ξ_lineshape::T
    ξ_ch::TwoBodyChannel
    L::Int = 0
    d::Float64 = 1.0
end

const QuasiTwoBodyChannelBW = QuasiTwoBodyChannel{BreitWigner};

threshold(ch::QuasiTwoBodyChannelBW) = real(ch.m1 + ch.ξ_ch.m1 + ch.ξ_ch.m2)
nominal_threshold(ch::QuasiTwoBodyChannelBW) = real(ch.m1 + ch.ξ_lineshape.m)

function real_ρ(ch::QuasiTwoBodyChannelBW, m::Real)
    @unpack ξ_lineshape, ξ_ch, L = ch
    σ_thr = threshold(ξ_ch)^2
    σ_hthr = (m - ch.m1)^2
    io = 0.0im
    # 
    @unpack l = ch.ξ_lineshape
    d_ij = ch.ξ_lineshape.d
    d_Rk = ch.d
    # 
    quadgk(σ_thr, σ_hthr) do σ
        _m2b = sqrt(σ)
        _ρ_ij = iρ(ξ_ch, _m2b) |> imag
        _ρ_Rk = iρ(TwoBodyChannel(ch.m1 + io, _m2b + io), m) |> imag
        #
        _p = (_ρ_ij * 2_m2b) |> real
        _q = (_ρ_Rk * m) |> real
        #
        ff = BlattWeisskopf{l}(d_ij)(_p)^2 * BlattWeisskopf{L}(d_Rk)(_q)^2
        # 
        abs2(ξ_lineshape(σ)) * _ρ_ij * _ρ_Rk * ff
    end[1]
end
