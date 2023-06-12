module QuantumMechanics

using LinearAlgebra

export Basis,
    GenericBasis, CompositeBasis, FockBasis, SpinBasis, tensor, ⊗, partialtrace, flatten
export Ket, Bra, dagger

include("bases.jl")
include("states.jl")

end
