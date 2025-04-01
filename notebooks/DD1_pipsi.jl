### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 8225ce1e-ee79-4606-bb22-fc46b52d1002
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__))
	Pkg.instantiate()
	# 
	using KMatrix.HadronicLineshapes
	using KMatrix.StaticArrays
	using KMatrix.Parameters
	using KMatrix.QuadGK
	using Interpolations
	using DelimitedFiles
	using KMatrix
	using Plots
end

# ╔═╡ a27fb311-1a58-4198-b992-fe16b46e7cb8
md"""
# Flatte amplitude with subchannel resonances
(Misha Mikhasenko, RUB, 1/04/2025)

In this notebook we construct two-channel K-matrix with a single pole.
An accurate construction with a dispersion integral is compared to an approximate formulas using complex mass `m-i\Gamma/2`.


## Main findings

We find that while shapes match well for resonance above both threshold, the shapes devitate in the Flatte regime. It's clear that an approximate expression does not give a good approximation for the phase space below the threshold. It biases the width of the resulting state, since

$\Gamma \approx \frac{1}{m} (g_1^2\rho_1 + g_2^2 \rho_2)$

So the width does not come the same in two settings for the Flatte regime.

## Nontheless

When dealing with a single channel, the bias in the width can well be compensated by the contribution of the channel with the lower threshold ($g_1^2 \rho_1$).
It seems to me that appriximate construction would span a similar class of curves, but might get biased values of the poles and residues. I'd still guess that one can start fitting with the approximate approach to get a feeling if data are sensitive to the second channel, and move to the advanced, accurate parametrization to quantify the outcomes.

## Technical

Since the QTB phase space has a different scale at high energy, couplings have to be renormalized.

"""

# ╔═╡ 1792d675-9235-4823-9d72-9912a5a0da88
theme(:boxed, fontfamily="Computer Modern")

# ╔═╡ 0e9bd099-cd96-4991-8b98-7859925540e3
begin
	const m_D2600 = 2.6
	const Γ_D2600 = 0.14
	const mD = 1.86
	const mπ = 0.14
	const mjψ = 3.09
end;

# ╔═╡ a5128948-0272-4e1b-a553-4cfa66a91125
const support = (mjψ+mπ, 5.5)

# ╔═╡ 43482761-ca9a-45bc-86a3-159785c21701
const iϵ = 1e-7im

# ╔═╡ 9336b016-0ed3-11f0-2871-0377facf61a3
A_approx = let
	channels = SVector(
        TwoBodyChannel(mjψ, mπ),
        TwoBodyChannel(mD, m_D2600-1im*Γ_D2600/2)
    )
    MG = [(M = 4.2, gs = [0.5, 4.4])]
    K = Kmatrix(MG)
    T = Tmatrix(K, channels)
	ProductionAmplitude(T, SVector(1.0), SVector(0.0, 0.0))
end

# ╔═╡ 8916df22-2961-4bd1-87f4-b32fad7473f5
A = let
	ch_Dπ = TwoBodyChannel(mD, mπ)
	ch_ctb = TwoBodyChannel(mD, m_D2600-0.5im*Γ_D2600)
	# 
	BW_D2600 = BreitWigner(; m=m_D2600, Γ=Γ_D2600, ma=mD, mb=mπ, l=0, d=1.0)
	ch_qtb = QuasiTwoBodyChannel(mD+0im, BW_D2600, ch_Dπ, 0, 1.0);
	ch_qtb_approx = InterpolatedChannel(ch_qtb, ch_ctb;
			m_grid_max = 5.9,
			n_grid = 200);
	#
	channels = SVector(
        TwoBodyChannel(mjψ, mπ),
        ch_qtb_approx
    )
	ρ_ratio = channels[2].inter.coefs[end]
	MG = [(M = A_approx.T.K.poles[1].M, gs = A_approx.T.K.poles[1].gs .* [1, 1/sqrt(ρ_ratio)])]
	K = Kmatrix(MG)
    T = Tmatrix(K, channels)
	ProductionAmplitude(T, SVector(1.0), SVector(0.0, 0.0))
end

# ╔═╡ 17a53e2b-af6e-4bc1-85b3-3bd2f2f27f0b
md"""
### Save look up tables to disc
"""

# ╔═╡ b119c99d-24d1-4699-9b50-97bc83e2d873
let # in three columns
	mv = range(support..., 100)
	Av = iρ.(Ref(A_approx.T.channels[2]), mv)
	data = [mv real.(Av) imag.(Av)]
	writedlm("D1D_orho_reim.txt", data)
end

# ╔═╡ 2597c296-3e91-4123-9fbf-fb0581c2e60e
md"""
## Plot spectra
"""

# ╔═╡ bca30e9e-c23c-40e5-9a2a-15a885ef953e
let
	plot(title = "approximate parametrization in Flatte regime", xlab="m [GeV]", ylab="\$\\mathcal{A}\$")
	plot!(m->abs2(amplitude(A_approx, m)[1]), support...)
	vline!([A_approx.T.K.poles[1].M], lab="bare mass", alpha=0.3)
	vline!([threshold(A_approx.T.channels[2])], lab="nominal threshold")
end

# ╔═╡ 7e64d6e4-2736-40e7-b7bf-262dfd89996f
let
	plot(title = "accurate parametrization in Flatte regime", xlab="m [GeV]", ylab="\$\\mathcal{A}\$")
	plot!(m->abs2(amplitude(A, m+iϵ)[1]), support...)
	vline!([A.T.K.poles[1].M], lab="bare mass", alpha=0.3)
	vline!([threshold(A_approx.T.channels[2])], lab="threshold")
end

# ╔═╡ Cell order:
# ╟─a27fb311-1a58-4198-b992-fe16b46e7cb8
# ╠═8225ce1e-ee79-4606-bb22-fc46b52d1002
# ╠═1792d675-9235-4823-9d72-9912a5a0da88
# ╠═0e9bd099-cd96-4991-8b98-7859925540e3
# ╠═a5128948-0272-4e1b-a553-4cfa66a91125
# ╠═43482761-ca9a-45bc-86a3-159785c21701
# ╠═9336b016-0ed3-11f0-2871-0377facf61a3
# ╠═8916df22-2961-4bd1-87f4-b32fad7473f5
# ╟─17a53e2b-af6e-4bc1-85b3-3bd2f2f27f0b
# ╠═b119c99d-24d1-4699-9b50-97bc83e2d873
# ╟─2597c296-3e91-4123-9fbf-fb0581c2e60e
# ╟─bca30e9e-c23c-40e5-9a2a-15a885ef953e
# ╟─7e64d6e4-2736-40e7-b7bf-262dfd89996f
