module KMatrix

using LinearAlgebra
using StaticArrays
using Parameters

export TwoBodyChannel, ρ
include("two-body-channel.jl")

export Kmatrix, amplitude, npoles, nchannels
export Tmatrix, Dmatrix, detD, channels
include("t-matrix.jl")

export ProductionAmplitude
export productionpole, productionnonpole
include("production-amplitude.jl")

end # module KMatrix
