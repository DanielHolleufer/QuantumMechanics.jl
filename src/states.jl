abstract type QuantumState end
abstract type StateVector <: QuantumState end

"""
    Ket(b::Basis, data::AbstractArray{Number})

Create a ket in the basis `b` with vector representation given by `data`.

The vector representation is converted to a vector whose elements are complex floats.
It also need not be normalized.

# Examples
```jldoctest
julia> Ket(FockBasis(4), [1, 0, 0, 0, 0])
Ket{FockBasis{Int64}, ComplexF64}(FockBasis{Int64}((5,), 4, 0), ComplexF64[1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im])

julia> Ket(SpinBasis(1 // 2), [1/sqrt(2), im/sqrt(2)])
Ket{SpinBasis{1//2, Int64}, ComplexF64}(SpinBasis{1//2, Int64}((2,), 1//2), ComplexF64[0.7071067811865475 + 0.0im, 0.0 + 0.7071067811865475im])
```
"""
mutable struct Ket{B,T} <: StateVector
    basis::B
    data::Vector{T}
    function Ket(b::B, data::AbstractArray{T}) where {B<:Basis,T<:Number}
        if length(b) !== length(data)
            error("Dimension of data does not match dimension of basis.")
        end
        return new{B,complex(float(T))}(b, vec(data))
    end
end
Base.:(==)(u::Ket, v::Ket) = u.data == v.data && u.basis == v.basis

"""
    Bra(b::Basis, data::AbstractArray{Number})

Create a bra in the basis `b` with vector representation given by `data`.

The vector representation is converted to a vector whose elements are complex floats.
It also need not be normalized.

# Examples
```jldoctest
julia> Bra(FockBasis(4), [1, 0, 0, 0, 0])
Bra{FockBasis{Int64}, ComplexF64}(FockBasis{Int64}((5,), 4, 0), ComplexF64[1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im])

julia> Bra(SpinBasis(1 // 2), [1/sqrt(2), im/sqrt(2)])
Bra{SpinBasis{1//2, Int64}, ComplexF64}(SpinBasis{1//2, Int64}((2,), 1//2), ComplexF64[0.7071067811865475 + 0.0im, 0.0 + 0.7071067811865475im])
```
"""
mutable struct Bra{B,T} <: StateVector
    basis::Basis
    data::Vector{T}
    function Bra(b::B, data::AbstractArray{T}) where {B<:Basis,T<:Number}
        if length(b) !== length(data)
            error("Dimension of data does not match dimension of basis.")
        end
        return new{B,complex(float(T))}(b, vec(data))
    end
end
Base.:(==)(u::Bra, v::Bra) = u.data == v.data && u.basis == v.basis

Base.:(*)(u::Bra, v::Ket) = dot(u.data, v.data)

"""
    dagger(ψ::StateVector)

Hermitian transpose of the state vector. Using `dagger` on a `Ket` or a `Bra` will
respectively return a `Bra` or a `Ket` in the same basis as the input state, but whose data
has been complex conjugated.
"""
dagger(ψ::Ket) = Bra(ψ.basis, ψ.data')
dagger(ψ::Bra) = Ket(ψ.basis, ψ.data')
Base.adjoint(ψ::StateVector) = dagger(ψ)