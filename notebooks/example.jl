### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ ac9d90ae-0e35-11f0-28b4-8f45c9645dc1
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__))
	Pkg.instantiate()
	# 
	using HadronicLineshapes
	using Interpolations
	using Parameters
	using KMatrix
	using QuadGK
	using Plots
end

# ╔═╡ 2184fd62-8a98-4db6-afc2-9387d67e70e8
theme(:boxed, fontfamily="Computer Modern")

# ╔═╡ b742b299-71fc-4262-ab64-2f27b02a365c
begin
	const m_D2600 = 2.6
	const Γ_D2600 = 0.1
	const mD = 1.86
	const mπ = 0.14
end;

# ╔═╡ 27a17d49-33b8-4c96-b857-dacb9bb68a21
const iϵ = 1e-7im

# ╔═╡ 8f0d2d81-18d5-4217-b763-463da2dfe8d7
const support = (2mD+mπ, 5.4);

# ╔═╡ c465501d-7880-433b-8660-0d89337323ad
ch_Dπ = TwoBodyChannel(mD, mπ)

# ╔═╡ c7947c0b-3d4b-47e7-8af9-1607d0c9b203
ch_tb = TwoBodyChannel(mD, m_D2600)

# ╔═╡ cd53e797-7a04-48f5-8b03-e0e9e390f3f8
ch_ctb = TwoBodyChannel(mD, m_D2600-0.5im*Γ_D2600)

# ╔═╡ e530f6bc-b85e-4e4c-a30a-3ab69458b828
const BW_D2600 = BreitWigner(; m=m_D2600, Γ=Γ_D2600, ma=mD, mb=mπ, l=0, d=1.0)

# ╔═╡ 7e2cf6bc-d33e-4e0b-a5d8-60875c0511a2
ch_qtb = QuasiTwoBodyChannel(mD+0im, BW_D2600, ch_Dπ, 0, 1.0);

# ╔═╡ 48234119-2887-476c-a159-4fedf388db16
ch_qtb_approx = InterpolatedChannel(ch_qtb, ch_ctb;
			m_grid_max = 5.9,
			n_grid = 100);

# ╔═╡ a424b1fa-fd2a-4cb2-8a2e-cc7174a31415
md"""
## Comparison
"""

# ╔═╡ 8bdf9a7c-702b-46b6-b4ab-7c8c5eba3b3c
let
	plot(title="Phase space volume", xlab="m [GeV]", ylab="rho", leg=:left)
	plot!(m->imag(iρ(ch_tb, m)), support..., lab="TB")
	plot!(m->imag(iρ(ch_ctb, m)), support..., lab="CTB")
	mh = 5.5
	n_qtb = imag(iρ(ch_tb, mh)) / real(real_ρ(ch_qtb, mh))
	# 
	plot!(m->imag(iρ(ch_qtb_approx, m+iϵ)) * n_qtb, support...,
		lab="approx QTB", ls=:dash, lw=3)
	plot!(m->imag(1im*real_ρ(ch_qtb, m)) * n_qtb, support..., lab="QTB")
	# 
end

# ╔═╡ 15fdb2e5-9892-4706-8f78-a5e47a553e4d
let
	plot(title="Phase space volume", xlab="m [GeV]", ylab="rho", leg=:left)
	plot!(m->real(iρ(ch_qtb_approx, m+iϵ)), support..., lab="real")
	plot!(m->imag(iρ(ch_qtb_approx, m+iϵ)), support..., lab="imag")
end

# ╔═╡ 92cb358f-499f-4bd4-ae7f-74e316d5ce13
let
	plot(title="Phase space volume", xlab="m [GeV]", ylab="rho", leg=:bottom)
	plot!(m->real(iρ(ch_tb, m+iϵ)), 0.7, 5.4, lab="\$\\Re\\,(i\\rho)\$")
	# 
	plot!(m->imag(iρ(ch_tb, m+iϵ)), 0.7, 5.4, lab="\$\\Im\\,(i\\rho)\$")
end

# ╔═╡ c09c1d2e-44ac-4789-a406-f921a0121aca


# ╔═╡ f3eeaaa5-591a-4731-8c97-ce54c6c0a438
begin
	animation = @animate for i in vcat(range(-4, 0, 10), range(0, -4, 10))
		eff = exp(i)
		BW_D2600 = BreitWigner(; m=m_D2600, Γ=Γ_D2600*eff, ma=mD, mb=mπ, l=0, d=1.0)
		ch_qtb = QuasiTwoBodyChannel(mD+0im, BW_D2600, ch_Dπ, 0, 1.0);
		ch_qtb_approx = InterpolatedChannel(ch_qtb, ch_ctb;
					m_grid_max = 5.9,
					n_grid = 200);
		# 
		plot(title="\$i\\rho\$ in K-matrix construction", xlab="m [GeV]", ylab="\$i \\rho\$", leg=:topleft, ylim=(-1,1))
		r = imag(iρ(ch_qtb_approx, 5.5+iϵ))
		mv = range(support..., 110)
		plot!(m->real(iρ(ch_qtb_approx, m+iϵ))/r, mv, lab="Real")
		plot!(m->imag(iρ(ch_qtb_approx, m+iϵ))/r, mv, lab="Imag")
	end
	gif(animation, "anim_fps15.gif", fps = 5)
end

# ╔═╡ Cell order:
# ╠═ac9d90ae-0e35-11f0-28b4-8f45c9645dc1
# ╠═2184fd62-8a98-4db6-afc2-9387d67e70e8
# ╠═b742b299-71fc-4262-ab64-2f27b02a365c
# ╠═27a17d49-33b8-4c96-b857-dacb9bb68a21
# ╠═8f0d2d81-18d5-4217-b763-463da2dfe8d7
# ╠═c465501d-7880-433b-8660-0d89337323ad
# ╠═c7947c0b-3d4b-47e7-8af9-1607d0c9b203
# ╠═cd53e797-7a04-48f5-8b03-e0e9e390f3f8
# ╠═7e2cf6bc-d33e-4e0b-a5d8-60875c0511a2
# ╠═e530f6bc-b85e-4e4c-a30a-3ab69458b828
# ╠═48234119-2887-476c-a159-4fedf388db16
# ╟─a424b1fa-fd2a-4cb2-8a2e-cc7174a31415
# ╠═8bdf9a7c-702b-46b6-b4ab-7c8c5eba3b3c
# ╠═15fdb2e5-9892-4706-8f78-a5e47a553e4d
# ╠═92cb358f-499f-4bd4-ae7f-74e316d5ce13
# ╠═c09c1d2e-44ac-4789-a406-f921a0121aca
# ╠═f3eeaaa5-591a-4731-8c97-ce54c6c0a438
