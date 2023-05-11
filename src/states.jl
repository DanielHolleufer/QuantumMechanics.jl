abstract type QuantumState end
abstract type StateVector <: QuantumState end

mutable struct Ket
    basis::Basis
    data::Vector{ComplexF64}
end

Base.:(==)(u::Ket, v::Ket) = u.data == v.data && u.basis == v.basis