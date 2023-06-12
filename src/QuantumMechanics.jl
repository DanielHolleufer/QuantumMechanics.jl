module QuantumMechanics

using LinearAlgebra

export Basis,
    GenericBasis, CompositeBasis, FockBasis, SpinBasis, tensor, ⊗, partialtrace, flatten
export Ket, Bra, copy, dagger, norm, normalize, normalize!

include("bases.jl")
include("states.jl")

end
