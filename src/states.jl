abstract type QuantumState end
abstract type StateVector <: QuantumState end

mutable struct Ket <: StateVector
    basis::Basis
    data::Vector{ComplexF64}
end
Base.:(==)(u::Ket, v::Ket) = u.data == v.data && u.basis == v.basis


mutable struct Bra <: StateVector
    basis::Basis
    data::Vector{ComplexF64}
end
Base.:(==)(u::Bra, v::Bra) = u.data == v.data && u.basis == v.basis

Base.:(*)(u::Bra, v::Ket) = dot(u.data, v.data)