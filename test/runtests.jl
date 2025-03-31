using KMatrix
using Test
using KMatrix.StaticArrays


@testset "TwoBodyChannel" begin
    ch = TwoBodyChannel(1.1, 1.1)
    @test ch.m1 == 1.1 + 0.0im
    @test ch.m2 == 1.1 + 0.0im
    @test ch.L == 0

    # Test phase space calculation
    m = 3.0
    iρval = iρ(ch, m)
    @test iρval isa Complex
    @test abs(iρval) > 0

    # Test below threshold
    @test imag(iρ(ch, 1.0)) ≈ 0 atol = 1e-10
end

@testset "3x3 one-pole K-matrix" begin
    channels = SVector(
        TwoBodyChannel(1.1, 1.1),
        TwoBodyChannel(2.2, 2.2),
        TwoBodyChannel(1.3, 1.3),
    )
    MG = [(M = 5.3, gs = [1.2, 0.5, 1.6])]
    K = Kmatrix(MG)
    T = Tmatrix(K, channels)

    @test nchannels(T) == 3
    @test npoles(T) == 1

    # Test amplitude calculation
    m = 5.0
    A = amplitude(T, m)
    @test size(A) == (3, 3)
    @test A ≈ transpose(A) # Check hermiticity
end

@testset "2x2 two-pole K-matrix" begin
    channels = SVector(
        TwoBodyChannel(1.1, 1.1),
        TwoBodyChannel(1.3, 1.3),
    )
    MG = [
        (M = 4.3, gs = [2.1, 0.0]),
        (M = 6.3, gs = [0.0, 2.5]),
    ]
    K = Kmatrix(MG)
    T = Tmatrix(K, channels)

    @test nchannels(T) == 2
    @test npoles(T) == 2

    # Test decoupled case
    m = 5.0
    A = amplitude(T, m)
    @test abs(A[1, 2]) ≈ 0 atol = 1e-10
    @test abs(A[2, 1]) ≈ 0 atol = 1e-10
end

@testset "Production amplitude" begin
    channels = SVector(TwoBodyChannel(1.1, 1.1))
    MG = [
        (M = 4.3, gs = [2.1]),
        (M = 6.3, gs = [2.5]),
    ]
    K = Kmatrix(MG)
    T = Tmatrix(K, channels)

    # Test with default production couplings
    A = ProductionAmplitude(T)
    @test length(A.αpoles) == 2
    @test all(A.αpoles .== 1)

    # Test with custom production couplings
    α = SVector(1.0, 2.0 * cis(π / 4))
    A_custom = ProductionAmplitude(T, α, SVector(0.0))

    m = 5.0
    amp = amplitude(A_custom, m)
    @test length(amp) == 1

    # Test individual pole contributions
    amp1 = production_pole(A_custom, m, 1)
    amp2 = production_pole(A_custom, m, 2)
    @test amp1 ≈ [-0.3068709069551295 + 0.06943242478375794im]
    @test amp2 ≈ [0.2807585353071957 + 0.17715198048327294im]
end
