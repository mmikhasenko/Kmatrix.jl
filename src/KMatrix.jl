module KMatrix

using LinearAlgebra
using StaticArrays
using Parameters
using HadronicLineshapes
using QuadGK

export TwoBodyChannel, iρ
export threshold
include("two-body-channel.jl")

export real_ρ
export nominal_threshold
export QuasiTwoBodyChannel, QuasiTwoBodyChannelBW
include("quasi-two-body.jl")

export interpolatedChannel
include("interpolation.jl")

export Kmatrix, amplitude, npoles, nchannels
export Tmatrix, Dmatrix, detD, channels
include("t-matrix.jl")

export ProductionAmplitude
export production_pole, production_nonpole
include("production-amplitude.jl")

end # module KMatrix
