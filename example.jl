### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 2a5d2458-9958-4c9b-a183-d7029b6360c9
# ╠═╡ show_logs = false
begin
    using StaticArrays
    using Plots
    using LaTeXStrings
    using PlutoUI
    using Parameters
    using LinearAlgebra
    using PlutoTeachingTools
    using CalculusWithJulia
end

# ╔═╡ 9d9a89ed-76f3-4e76-a3cc-0d33c747fbb5
answer_box(md"""
In that case, the T-matrix would be simply:

$T =
\begin{pmatrix}
\frac{g_1^2}{m_{(1)}^2-s-ig_1^2\rho_1} & 0\\
0 & \frac{h_2^2}{m_{(2)}^2-s-ih_2^2\rho_2}
\end{pmatrix}$
""")

# ╔═╡ 12a615bb-97b8-4fde-bd66-ac7083970e0e
RobustLocalResource("", joinpath("..", "figures", "1x1_production.png"), cache = false)

# ╔═╡ 0360447c-c6c7-4ba6-8e5f-a20d5797995b
begin
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
    #
    amplitude(K::Kmatrix, m) =
        sum((gs * gs') ./ (M^2 - m^2) for (M, gs) in K.poles) + K.nonpoles
    # 
    npoles(X::Kmatrix{N, V}) where {N, V} = V
    nchannels(X::Kmatrix{N, V}) where {N, V} = N
    #
end

# ╔═╡ 25758574-5da4-4fbe-946d-d55f6210b7e2
theme(:wong2, frame = :box, grid = false, minorticks = true,
    guidefontvalign = :top, guidefonthalign = :right,
    foreground_color_legend = nothing,
    xlim = (:auto, :auto), ylim = (:auto, :auto),
    lab = "", xlab = "Mass of system [GeV]")

# ╔═╡ 0e899f67-adec-4837-9993-c9fe22f788d1
md"""
## Example 3: Single-Channel Production Amplitude with Interfering Resonances

The third example examines a single-channel problem featuring two prominent resonances. The focus is on the production amplitude and the investigation of interference effects between resonances. This scenario is critical for comprehending how overlapping resonances interact within the K-matrix formalism, affecting the overall scattering amplitude and observable patterns in the data.
"""

# ╔═╡ 8b53b1cc-520b-48e4-b2b1-6ad6ebe443e2
Markdown.MD(Markdown.Admonition("danger", "Question", [md"""
**E2.Q3:** Let's , finally, put $g_2=h_1=\epsilon$.

What will we see for the off-diagonal terms of T-matrix?

"""]))

# ╔═╡ 729c86ab-cf81-48bf-82be-b89cf28eaee6
md"""

### Demonstration

- g₂ = $(@bind g2_T2 Slider(range(-1,3,81), default=0.0, show_value=true)) GeV for the first channel,
- h₁ = $(@bind h1_T2 Slider(range(-1,3,81), default=0.0, show_value=true)) GeV for the second channel, and 
"""

# ╔═╡ 91db142a-109f-414a-8d2c-9d3cd92bae40
md"""
### Production amplitude

The difference between production and scattering is that now one replace K-matrix with vector:

$A = [1-iK\rho]^{-1}
\begin{pmatrix}
\frac{\alpha_1 g}{m_{(1)}^2-s}\\
\frac{\alpha_2 h}{m_{(2)}^2-s}
\end{pmatrix}$

Where $a_1$ and $a_2$ are production factors which might be complex.
"""

# ╔═╡ 474e5d19-b537-42d2-9e92-2a26996cee2d
const iϵ = 1e-7im

# ╔═╡ 0f2ade1e-6958-4ebd-942a-c844bf3dbb99
md"""
### The case of weak coupling for one of the resonances
"""

# ╔═╡ ad74e82b-b2f4-4b8d-99ea-2fe295bb018d
md"""
### Calculation of the width

We can relate the approximate width to coupling using the Breit-Wigner expression, $1/(m_0^2-s-im_0\Gamma_0)$, where $\Gamma_0$ gives the width of the peak.
Hence

```math
\begin{align}
\Gamma_0 = \frac{g_1^2 \rho_1 (m_0) + g_2^2 \rho_2(m_0)+g_3^2 \rho_3(m_0)}{m_0}\\
\end{align}
```
"""

# ╔═╡ f9dad52e-d6a2-46c4-a5b3-91a50c9425c1
md"""
Production couplings
- α₁ = $(@bind α1_E3 Slider(range(0.01,2,61), default=1.0, show_value=true))
- |α₂| = $(@bind α2_E3 Slider(range(-1,2,61), default=1.0, show_value=true))
- Arg(α₂) = $(@bind ϕ2_E3 Slider(range(0,2π,100), default=0.0, show_value=true))
"""

# ╔═╡ 437d3bf9-cf44-40a8-97bc-8aee4b62d069
RobustLocalResource("", joinpath("..", "figures", "1x1_scattering.png"))

# ╔═╡ 27f0ae17-b61c-49c5-b4fc-6de5d2ddda94
md"""
### Scattering amplitude
"""

# ╔═╡ 59ceb096-1cc0-4c69-b56f-476710dd698e
Markdown.MD(Markdown.Admonition("danger", "Question", [md"""
**E1.Q2:** Which parameter is needed to be changed in order to see difference between BW and K-matrix formalism?
"""]))

# ╔═╡ 1273bd41-9986-4b9f-9e06-b3bed7ab65f0
answer_box(
    md"$T_{ij} = \frac{g_i g_j}{m_0^2-s-ig_1^2 \rho_1-ig_2^2 \rho_2-ig_3^2 \rho_3}$

    This is expression known as the [Flatte formula](https://inspirehep.net/literature/108884)
    ")

# ╔═╡ edef417d-b0e4-4cad-bf63-462a8d7e861f
Markdown.MD(Markdown.Admonition("danger", "Question", [md"""
**E3.Q1:** Why does it have zero between the poles?
"""]))

# ╔═╡ 6d5acd0c-dbcf-4d0a-a94c-76ac59006fc8
md"""

### Comparison to Breit-Wigner

Let's compare to the BW parametrization

$\text{BW} = \frac{m_0\Gamma_0}{m_0^2-s-im_0 \Gamma_0}$

"""

# ╔═╡ cf9733b3-bdc9-4e58-a7f3-87845eb907da
Markdown.MD(Markdown.Admonition("danger", "Question", [md"""
**E3.Q2:** What can the complexity of a production factor lead to? 
How the plot above will change?
"""]))

# ╔═╡ ae68a0e2-46ce-4186-97b3-4b03b5f2d8ce
begin
    struct TwoBodyChannel
        m1::Complex{Float64}
        m2::Complex{Float64}
        L::Int
    end
    TwoBodyChannel(m1, m2; L::Int = 0) = TwoBodyChannel(m1, m2, L)
    # 
    function ρ(ch::TwoBodyChannel, m; ϕ = -π / 2)
        ch.L != 0 && error("not implemented")
        sqrt(cis(ϕ) * (m - (ch.m1 + ch.m2))) * cis(-ϕ / 2) *
        sqrt(m + (ch.m1 + ch.m2)) *
        sqrt((m^2 - (ch.m1 - ch.m2)^2)) /
        m^2
    end
end

# ╔═╡ 7b717c8f-1fb8-4892-a250-c77e5e088445
aside(Markdown.MD(Markdown.Admonition("warning", "Tip", [md"$T = [1-iK\rho]^{-1}K$"])))

# ╔═╡ 798f53e9-d871-42d1-a81c-d35adc7ece21
md"""
### Setup

Let's consider 3 coupled channels, the scattering amplitude T is a 3x3 matrix:

$T = \begin{pmatrix}
T_{1\to1} & T_{1\to2} & T_{1\to3}\\
T_{2\to1} & T_{2\to2} & T_{2\to3}\\
T_{3\to1} & T_{3\to2} & T_{3\to3}\\
\end{pmatrix}$

We model it with the simplest **one-pole** K-matrix.

$K = \frac{1}{m_0^2-s}
\begin{pmatrix}
g_1^2 & g_1 g_2 & g_1 g_3\\
g_2 g_1 & g_2^2 & g_2 g_3\\
g_3 g_1 & g_3 g_2 & g_3^2\\
\end{pmatrix}$ 

and $\rho=\text{Diag}(\rho_1,\rho_2,\rho_3)$
"""

# ╔═╡ fff6f9f5-b002-46ea-b6ff-1ecf32357ea9
md"""
### The case of weak coupling for both resonances
"""

# ╔═╡ cb9fdca3-db42-4501-a2fe-c3ff099b2d83
TableOfContents()

# ╔═╡ 87fc0818-273b-4d0b-814a-058365ee07a0
answer_box(
    md" $g_{i}\approx\sqrt\frac{Г_{0}m_{0}}{\rho_{i}(m_0)}$
    ")

# ╔═╡ 3019d77e-a41e-4fa5-a0bf-b91d3d72e96f
md"""
## Implementation
"""

# ╔═╡ e8ba4b4e-2d21-4bfe-bf32-02969f9b2970
md"""
The provided examples are designed to facilitate this understanding by illustrating the connection between the K-matrix parameters and observable quantities.


**Exapmle 1 (3x3 one pole):** The connection between the K-matrix and the multichannel Breit-Wigner model, demonstrating the generation of resonance widths.

**Example 2 (2x2 two poles):** The analysis of a 2x2 problem highlighting the dynamics of weakly coupled resonances.

**Example 3 (1x1 two poles):** The examination of a single-channel production amplitude involving two interfering resonances.

These examples provide insights into estimating parameters, understanding resonance behavior, and interpreting observable phenomena within the K-matrix framework.
"""

# ╔═╡ 9ade0bbc-0b04-4e6d-92cf-54062b638cfe
md"""
# Practical K-Matrix Applications

## Introduction

This educational material introduces the K-matrix formalism, focusing on its practical application in describing scattering amplitude observables. The K-matrix approach provides a robust framework for analyzing scattering processes, essential for understanding resonance phenomena in hadron physics. One of the main challenges in employing this formalism lies in the initial estimation of the bare parameters, which requires a thorough comprehension of the underlying principles and mechanics of the model.

On the figure below one can find illustration of the ab->cd scattering process with resonance which has mass $m_{0}$ and width $Г_{0}$. This process can be described with K-matrix:

$K_{ij} = \frac{g_i g_j}{m_0^2-s}$

where $g_{i} g_j$ indicates different possible channels of the scattering process and $s=m^{2}$. K-matrix has to be real and symmetric by construction.

Then, one can introduce transition T-matrix:

$T = [I -i  K \rho]^{-1} K$
where $I$ is identity matrix of three dimensions and $\rho$ being a diagonal matrix of phase space factors $\rho=Diag(\rho_i)$.
"""

# ╔═╡ c7e615fa-62aa-4de2-8ef4-2df8534b2c06
Markdown.MD(Markdown.Admonition("danger", "Question", [md"""
**E2.Q1:** Let's start with setting off-diagonal parameters described coupling to zeros $g_{2}=h_{1}=0$. Then K-matrix would be What would be 

$K =
\begin{pmatrix}
\frac{g_1^2}{m_{(1)}^2-s} & 0\\
0 & \frac{h_2^2}{m_{(2)}^2-s}
\end{pmatrix}$

How will the T-matrix look like?
"""]))

# ╔═╡ 1410b82a-0017-4ef3-adf8-1f4da66393a4
answer_box(
    md" $g_{2}$ and $g_{3}$
    ")

# ╔═╡ 915f987d-9bb5-4e0b-9cf0-f52e3937695a
md"""
## Example 2: The 2x2 Problem with Weakly Coupled Resonances

In the second example, we explore the 2x2 problem, focusing on scenarios with weakly coupled resonances. This case study sheds light on the interactions between resonances in a two-dimensional parameter space, emphasizing the effects of weak coupling. It is an essential exploration for understanding how resonance coupling influences scattering amplitudes and observable resonance characteristics.
"""

# ╔═╡ 1215bc85-4760-473b-93d0-5d6a8952e27e
Markdown.MD(Markdown.Admonition("danger", "Question", [md"""
**E1.Q1:** When $K$ is degenerate and has rank 1, the expression for T is simple. Figure it out for 3x3 matrix.

**Advanced option:** Can you prove it for general case?
"""]))

# ╔═╡ a1558b96-576f-4661-b357-c9f036c0167d
answer_box(md"""
The T-matrix in that case would be:

$T \approx
\begin{pmatrix}
\frac{g_1^2}{m_{(1)}^2-s-ig_1^2\rho_1} & \frac{g_1\epsilon}{m_{(1)}^2-s-ig_1^2\rho_1}\\
\frac{\epsilon g_1}{m_{(1)}^2-s-ig_1^2\rho_1} & \frac{h_2^2}{m_{(2)}^2-s-ih_2^2\rho_2}
\end{pmatrix}$

One can notice that:

(1) $T^{(0)}_{11}=T^{(\epsilon)}_{11}$ and $T^{(0)}_{22}=T^{(\epsilon)}_{22}$

(2) $\frac{T^{(\epsilon)}_{12}}{T^{(\epsilon)}_{11}}=\frac{\epsilon}{g_{1}}$

(here $T^{(0)}$ and $T^{(\epsilon)}$ indicate previous and this case respectively)

Hence one will see effects of channel couplings only when $g_{2}$ will be in order of $g_{1}$

""")

# ╔═╡ 1e7ebcaf-dd83-40fb-9bd8-2e48a1911bfa
answer_box(md"""

$T_{12}=T_{21}=\frac{g_1\epsilon}{m_{(1)}^2-s-ig_1^2\rho_1}+\frac{h_2\epsilon}{m_{(2)}^2-s-ih_2^2\rho_2}=(BW_{g_{1}}+BW_{h_{2}})\epsilon$

""")

# ╔═╡ a223cbff-88b0-4a28-af53-c139e7b9108a
Markdown.MD(Markdown.Admonition("warning", "Tip", [md"For single resonance:

$T =
\frac{1}{m_{(1)}^2-s-ig_1^2\rho_1-ig_2^2\rho_2}
\begin{pmatrix}
g_1^2 & g_1g_2\\
g_2g_1 & g_2^2
\end{pmatrix}$

"]))

# ╔═╡ fe35af83-4910-48e4-b9de-5b8a1f85fb72
md"""
Cross section of certain channel calculated as 

$\frac{\mathrm{d}\sigma}{\mathrm{d}m} = \frac{1}{J_i} |T_{ij}|^2 \frac{\mathrm{d}\Phi_j}{\mathrm{d}m}$

with
- the flux $J_i$ set to 1, and
- the $\mathrm{d}\Phi_j / \mathrm{d}m = 2m\rho_j$ is a phase space per certain energy, and
- the $\rho_j$ is zero below the threshold of the j channel.
"""

# ╔═╡ cf9bf7ac-905f-4112-a7f4-36c536d33918
Markdown.MD(Markdown.Admonition("danger", "Question", [md"""
**E1.Q3:** Let's assume that we know cross-section distribution and can measure FWHM of it. How is it possible to estimate K-matrix parameters $g_{i}$?
"""]))

# ╔═╡ 74e06991-79a2-4711-8bc7-c8656249641f
md"""
### The case of not coupled resonances
"""

# ╔═╡ 3fde6651-a704-4757-b282-3a7cfcd36f6e
md"""
Let's setup couplings for three channels, that we can adjust:
- g₁ = $(@bind g1_T1 Slider(range(0,4,101), default=1.2, show_value=true)) GeV for the first channel,
- g₂ = $(@bind g2_T1 Slider(range(0,4,101), default=0.5, show_value=true)) GeV for the second channel, and 
- g₃ = $(@bind g3_T1 Slider(range(0,4,101), default=1.6, show_value=true)) GeV for the third channel.
- m₀ = $(@bind M Slider(range(0,10,101), default=5.3, show_value=true)) GeV/c² pole of the K-matrix.
"""

# ╔═╡ fb45f5a8-15c4-4695-b36b-f21aab1e3d80
begin
    struct ProductionAmplitude{N, V}
        T::Tmatrix{N, V}
        αpoles::SVector{V, <:Number}
        αnonpoles::SVector{N, <:Number}
    end
    #
    npoles(X::ProductionAmplitude{N, V}) where {N, V} = V
    nchannels(X::ProductionAmplitude{N, V}) where {N, V} = N
    detD(X::ProductionAmplitude, m; ϕ = -π / 2) = detD(X.T, m; ϕ)
    channels(X::ProductionAmplitude) = channels(X.T)
    # 
    ProductionAmplitude(T::Tmatrix{N, V}) where {N, V} =
        ProductionAmplitude(T, SVector{V}(ones(V)), SVector{N}(ones(N)))
    # 
    function amplitude(A::ProductionAmplitude, m; ϕ = -π / 2)
        @unpack T, αpoles, αnonpoles = A
        P = αnonpoles
        for (α, Mgs) in zip(αpoles, A.T.K.poles)
            @unpack M, gs = Mgs
            P += α .* gs ./ (M^2 - m^2)
        end
        D⁻¹ = inv(Dmatrix(T, m; ϕ))
        return D⁻¹ * P
    end
end

# ╔═╡ 1072afe6-bad7-4de3-9ad2-6dcef2b924bb
begin
    struct Tmatrix{N, V}
        K::Kmatrix{N, V}
        channels::SVector{N, TwoBodyChannel}
    end
    #	
    function Dmatrix(T::Tmatrix{N, V}, m; ϕ = -π / 2) where {N, V}
        𝕀 = Matrix(I, (N, N))
        iρv = 1im .* ρ.(T.channels, m; ϕ) .* 𝕀
        K = amplitude(T.K, m)
        D = 𝕀 - K * iρv
    end
    detD(T::Tmatrix, m; ϕ = -π / 2) = det(Dmatrix(T, m; ϕ))
    amplitude(T::Tmatrix, m; ϕ = -π / 2) = inv(Dmatrix(T, m; ϕ)) * amplitude(T.K, m)
    # 
    npoles(X::Tmatrix{N, V}) where {N, V} = V
    nchannels(X::Tmatrix{N, V}) where {N, V} = N
    channels(X::Tmatrix) = X.channels
end

# ╔═╡ 0ff7e560-37ca-4016-bc01-741322402679
T1 = let
    # thhree channels
    channels = SVector(
        TwoBodyChannel(1.1, 1.1),
        TwoBodyChannel(2.2, 2.2),
        TwoBodyChannel(1.3, 1.3))
    # one bare pole
    MG = [(M, gs = [g1_T1, g2_T1, g3_T1])]
    # 
    K = Kmatrix(MG)
    T = Tmatrix(K, channels)
end;

# ╔═╡ ce2c280e-6a55-4766-a0f9-941b448c41c9
begin
    Γv_T1 = map(zip(T1.K.poles[1].gs, T1.channels)) do (g, ch)
        m0 = T1.K.poles[1].M
        Γi = g^2 * ρ(ch, m0) / m0
        round(Γi, digits = 2)
    end |> real
    Markdown.parse("The sum: $(round(sum(Γv_T1), digits=2)) GeV = $(join(string.(Γv_T1), " + ")*" GeV"), that corresponds to $(join(string.(round.(100*Γv_T1./sum(Γv_T1), digits=1)) .*"%", ", ")), respectively.")
end

# ╔═╡ 9a2e8210-c140-4689-bb16-2aab3c3b2aaa
productionnonpole(T::Tmatrix{N, K}, m; ϕ = -π / 2) where {N, K} =
    inv(Dmatrix(T, m; ϕ)) * ones(N)
#

# ╔═╡ a6fd628c-86db-4d3f-836b-ff376cac7f1d
T2 = let
    # two channels
    channels = SVector(
        TwoBodyChannel(1.1, 1.1),
        TwoBodyChannel(1.3, 1.3))
    # two bare pole
    g1_T2 = 2.1
    h2_T2 = 2.5
    MG = [
        (M = 4.3, gs = [g1_T2, g2_T2]),
        (M = 6.3, gs = [h1_T2, h2_T2])]
    # 
    K = Kmatrix(MG)
    T = Tmatrix(K, channels)
end;

# ╔═╡ ebf8b842-99ce-40e8-9bf4-931471879bf9
function productionpole(T::Tmatrix, m, iR::Int; ϕ = -π / 2)
    @unpack M, gs = T.K.poles[iR]
    P = gs ./ (M^2 - m^2)
    return inv(Dmatrix(T, m; ϕ)) * P
end

# ╔═╡ fbfc0e6c-775e-4025-ad06-e3f6e291ec52
let
    plot(title = ["Scattering cross section" ""], yaxis = nothing,
        layout = grid(2, 1, heights = (0.5, 0.5)), size = (700, 500))
    plot!(2.6, 8, sp = 1, lab = "Tₖ[1,1]") do e
        A = amplitude(T2, e)[1, 1]
        phsp = real(ρ(T2.channels[1], e)) * e
        abs2(A) * phsp
    end
    vline!(sp = 1, [T2.K.poles[1].M])
    plot!(2.6, 8, sp = 2, lab = "Tₖ[2,2]") do e
        A = amplitude(T2, e)[2, 2]
        phsp = real(ρ(T2.channels[1], e)) * e
        abs2(A) * phsp
    end
    vline!(sp = 2, [T2.K.poles[2].M])
    plot!(2.6, 8, sp = 2, lab = "Tₖ[1,2]") do e
        A = amplitude(T2, e)[1, 2]
        phsp = real(ρ(T2.channels[1], e)) * e
        abs2(A) * phsp
        abs2(A) * phsp
    end
end

# ╔═╡ 35b9a451-5fac-46e2-97bb-ebc19d4b3418
function productionpole(A::ProductionAmplitude{N, V}, m, iR::Int; ϕ = -π / 2) where {N, V}
    αnonpoles = SVector{N}(zeros(N))
    αpoles = zeros(Complex{Float64}, V)
    αpoles[iR] = A.αpoles[iR]
    A = ProductionAmplitude(A.T, SVector{V}(αpoles), αnonpoles)
    return amplitude(A, m; ϕ)
end

# ╔═╡ 8b92df7f-d97b-43fa-8ac3-fed8ee974f5f
md"""
For this case K-matrix will become just a function of s:

$K = \frac{g^2}{m_{(1)}^2-s}+
\frac{h^2}{m_{(2)}^2-s}$

$T = [1-iK\rho ]^{-1} K$

If K is zero for $s=s_\text{z}$, T is zero.
"""

# ╔═╡ 53ad2e4e-0268-4488-9891-815922d8a8db
aside(Markdown.MD(Markdown.Admonition("warning", "Tip", [md"For explanation let's have a look at K-matrix in these case."])))

# ╔═╡ a7a68629-bf6f-435e-9a97-de9a02a31160
Markdown.MD(Markdown.Admonition("danger", "Question", [md"""
**E2.Q2:** Let's find the first expansion term to reflect on how a weakly coupled resonances show up in the second channel. For that put $g_2=\epsilon$, $h_1 = 0$.

Then,

$K^{(1)} = \frac{1}{m_{(1)}^2-s}
\begin{pmatrix}
g_1^2 & g_1 \epsilon\\
\epsilon g_1 & \epsilon^2
\end{pmatrix}
\approx \frac{1}{m_{(1)}^2-s}
\begin{pmatrix}
g_1^2 & g_1 \epsilon\\
\epsilon g_1 & 0
\end{pmatrix}$

In this case the T-matrix will become:
"""]))

# ╔═╡ babfbc1f-7beb-44d1-b3c8-75309e8b817c
T3 = let
    # one channels
    channels = SVector(
        TwoBodyChannel(1.1, 1.1))
    # two bare pole
    MG = [
        (M = 4.3, gs = [2.1]),
        (M = 6.3, gs = [2.5])]
    # 
    K = Kmatrix(MG)
    T = Tmatrix(K, channels)
end;

# ╔═╡ 8556c8f4-89e4-4544-ad73-4da8b43a7051
p = let
    f1(x) = 1 / (T3.K.poles[1].M - x)
    f2(x) = 1 / (T3.K.poles[2].M - x)
    f3(x) = f1(x) + f2(x)
    plot(title = "K-matrix parameter")
    plot!(rangeclamp(f1), 3, 7, lab = "First pole, K1")
    plot!(rangeclamp(f2), 3, 7, lab = "Second pole, K2")
    plot!(rangeclamp(f3), 3, 7, lab = "Total")
    vline!([T3.K.poles[1].M, T3.K.poles[2].M], linestyle = :dash)
    hline!([0], color = :black, size = (400, 300), leg = :top)
end;

# ╔═╡ e9c0d423-793d-4997-a0fd-5c67b41fffb2
answer_box(
    TwoColumn(md"There are two poles because of which both K-matrix have asymptotics which interfere and results in 0 between them", plot(p)),
)

# ╔═╡ 7b4fa73c-1075-4271-8d2f-1668d98904ab
let
    plot(title = "Scattering cross section", yaxis = nothing)
    plot!(2.6, 8, lab = "Tₖ") do e
        A = amplitude(T3, e)[1, 1]
        phsp = real(ρ(T3.channels[1], e)) * e
        abs2(A) * phsp
    end
    vline!([T3.K.poles[1].M, T3.K.poles[2].M])
end

# ╔═╡ ff61d6e1-5681-4d15-8164-62baee4613ba
A_E3 = ProductionAmplitude(T3, SVector{2}(α1_E3, α2_E3 * cis(ϕ2_E3)), SVector{1}(0));

# ╔═╡ b90a1c10-8e57-4b2a-b48a-ccb78010a4f5
let
    plot()
    plot!(title = "Production cross section", 2.6, 8, sp = 1, lab = "Total production", yaxis = nothing) do e
        A = amplitude(A_E3, e)[1]
        phsp = real(ρ(T3.channels[1], e)) * e
        abs2(A) * phsp
    end
    map([1, 2]) do ind
        plot!(2.6, 8, sp = 1, fill = 0, alpha = 0.2, lab = "from R$(ind)") do e
            A = productionpole(A_E3, e, ind)[1]
            phsp = real(ρ(T3.channels[1], e)) * e
            abs2(A) * phsp
        end
    end
    vline!([T3.K.poles[1].M, T3.K.poles[2].M])
    plot!()
end

# ╔═╡ 2aef1bd5-7a81-417b-a090-77644fc5f640
function productionnonpole(A::ProductionAmplitude{N, V}, m; ϕ = -π / 2) where {N, V}
    @unpack T, αnonpoles = A
    αpoles = SVector{V}(zeros(V))
    A = ProductionAmplitude(T, αpoles, αnonpoles)
    return amplitude(A, m; ϕ)
en

# ╔═╡ 2486eb34-a858-4ea9-99e1-f17627589461
RobustLocalResource("", joinpath("..", "figures", "2x2_scattering.png"))

# ╔═╡ a4750e66-b448-479e-a5dd-b9aec0f3a857
aside(RobustLocalResource("", joinpath("..", "figures", "3x3_scattering.png")))

# ╔═╡ c2a0a8bc-c1e6-4a48-91dc-590ca79383ff
md"""
## Example 1: Multichannel Breit-Wigner Model and Resonance Widths

The first example explains the relationship between the K-matrix and the multichannel Breit-Wigner (BW) model. It demonstrates how resonance widths are generated within this framework, offering insight into the dynamic nature of resonances. This example is pivotal for appreciating how the K-matrix formalism can be used to describe complex resonance structures in multichannel scattering processes.
"""

# ╔═╡ 1663348b-ed67-4851-9365-9641e6379fcd
md"""

The K-matrix for such case consists of two terms:

$K =K^{(1)}+K^{(2)}
= \frac{1}{m_{(1)}^2-s}
\begin{pmatrix}
g_1^2 & g_1 g_2\\
g_2 g_1 & g_2^2
\end{pmatrix} + 
\frac{1}{m_{(2)}^2-s}
\begin{pmatrix}
h_1^2 & h_1 h_2\\
h_2 h_1 & h_2^2
\end{pmatrix}$
"""

# ╔═╡ 3e76dfec-e83c-46df-8838-b299a7aaa5e3
let
    m0 = T1.K.poles[1].M
    plot(title = "Scattering cross section 1 → 1", yaxis = nothing)
    # K-matrix
    plot!(2.6, 8, lab = "Tₖ[1,1]") do e
        A = amplitude(T1, e)[1, 1]
        phsp = real(ρ(T1.channels[1], e)) * e
        abs2(A) * phsp
    end
    vline!([m0], lab = "m₀ (K-matrix pole)")
    # BW
    N_peak = abs2(amplitude(T1, m0 + iϵ)[1, 1])
    plot!(2.6, 8, lab = "BW") do e
        Γ0 = sum(Γv_T1)
        A = m0 * Γ0 / (m0^2 - e^2 - 1im * m0 * Γ0)
        phsp = real(ρ(T1.channels[1], e)) * e
        N_peak * abs2(A) * phsp
    end
    vline!(map(T1.channels) do ch
            real(ch.m1 + ch.m2)
        end, lab = "thresholds", ls = :dash)
end

# ╔═╡ Cell order:
# ╠═8556c8f4-89e4-4544-ad73-4da8b43a7051
# ╠═1072afe6-bad7-4de3-9ad2-6dcef2b924bb
# ╠═9d9a89ed-76f3-4e76-a3cc-0d33c747fbb5
# ╠═12a615bb-97b8-4fde-bd66-ac7083970e0e
# ╠═0360447c-c6c7-4ba6-8e5f-a20d5797995b
# ╠═25758574-5da4-4fbe-946d-d55f6210b7e2
# ╠═0e899f67-adec-4837-9993-c9fe22f788d1
# ╠═8b53b1cc-520b-48e4-b2b1-6ad6ebe443e2
# ╠═729c86ab-cf81-48bf-82be-b89cf28eaee6
# ╠═91db142a-109f-414a-8d2c-9d3cd92bae40
# ╠═474e5d19-b537-42d2-9e92-2a26996cee2d
# ╠═0f2ade1e-6958-4ebd-942a-c844bf3dbb99
# ╠═0ff7e560-37ca-4016-bc01-741322402679
# ╠═ad74e82b-b2f4-4b8d-99ea-2fe295bb018d
# ╠═f9dad52e-d6a2-46c4-a5b3-91a50c9425c1
# ╠═fbfc0e6c-775e-4025-ad06-e3f6e291ec52
# ╠═7b4fa73c-1075-4271-8d2f-1668d98904ab
# ╠═437d3bf9-cf44-40a8-97bc-8aee4b62d069
# ╠═27f0ae17-b61c-49c5-b4fc-6de5d2ddda94
# ╠═ce2c280e-6a55-4766-a0f9-941b448c41c9
# ╠═59ceb096-1cc0-4c69-b56f-476710dd698e
# ╠═2a5d2458-9958-4c9b-a183-d7029b6360c9
# ╠═1273bd41-9986-4b9f-9e06-b3bed7ab65f0
# ╠═edef417d-b0e4-4cad-bf63-462a8d7e861f
# ╠═9a2e8210-c140-4689-bb16-2aab3c3b2aaa
# ╠═a6fd628c-86db-4d3f-836b-ff376cac7f1d
# ╠═6d5acd0c-dbcf-4d0a-a94c-76ac59006fc8
# ╠═b90a1c10-8e57-4b2a-b48a-ccb78010a4f5
# ╠═cf9733b3-bdc9-4e58-a7f3-87845eb907da
# ╠═35b9a451-5fac-46e2-97bb-ebc19d4b3418
# ╠═ae68a0e2-46ce-4186-97b3-4b03b5f2d8ce
# ╠═7b717c8f-1fb8-4892-a250-c77e5e088445
# ╠═e9c0d423-793d-4997-a0fd-5c67b41fffb2
# ╠═798f53e9-d871-42d1-a81c-d35adc7ece21
# ╠═fff6f9f5-b002-46ea-b6ff-1ecf32357ea9
# ╠═ff61d6e1-5681-4d15-8164-62baee4613ba
# ╠═cb9fdca3-db42-4501-a2fe-c3ff099b2d83
# ╠═87fc0818-273b-4d0b-814a-058365ee07a0
# ╠═3019d77e-a41e-4fa5-a0bf-b91d3d72e96f
# ╠═e8ba4b4e-2d21-4bfe-bf32-02969f9b2970
# ╠═9ade0bbc-0b04-4e6d-92cf-54062b638cfe
# ╠═c7e615fa-62aa-4de2-8ef4-2df8534b2c06
# ╠═ebf8b842-99ce-40e8-9bf4-931471879bf9
# ╠═1410b82a-0017-4ef3-adf8-1f4da66393a4
# ╠═915f987d-9bb5-4e0b-9cf0-f52e3937695a
# ╠═1215bc85-4760-473b-93d0-5d6a8952e27e
# ╠═a1558b96-576f-4661-b357-c9f036c0167d
# ╠═fb45f5a8-15c4-4695-b36b-f21aab1e3d80
# ╠═1e7ebcaf-dd83-40fb-9bd8-2e48a1911bfa
# ╠═a223cbff-88b0-4a28-af53-c139e7b9108a
# ╠═fe35af83-4910-48e4-b9de-5b8a1f85fb72
# ╠═cf9bf7ac-905f-4112-a7f4-36c536d33918
# ╠═74e06991-79a2-4711-8bc7-c8656249641f
# ╠═3fde6651-a704-4757-b282-3a7cfcd36f6e
# ╠═8b92df7f-d97b-43fa-8ac3-fed8ee974f5f
# ╠═53ad2e4e-0268-4488-9891-815922d8a8db
# ╠═a7a68629-bf6f-435e-9a97-de9a02a31160
# ╠═babfbc1f-7beb-44d1-b3c8-75309e8b817c
# ╠═2aef1bd5-7a81-417b-a090-77644fc5f640
# ╠═2486eb34-a858-4ea9-99e1-f17627589461
# ╠═a4750e66-b448-479e-a5dd-b9aec0f3a857
# ╠═c2a0a8bc-c1e6-4a48-91dc-590ca79383ff
# ╠═1663348b-ed67-4851-9365-9641e6379fcd
# ╠═3e76dfec-e83c-46df-8838-b299a7aaa5e3
