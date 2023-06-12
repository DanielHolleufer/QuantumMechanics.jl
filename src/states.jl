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
    function Ket{B,T}(
        b::B, data::AbstractArray{<:Number}
    ) where {B<:Basis,T<:Complex{<:AbstractFloat}}
        if length(b) !== length(data)
            error("Dimension of data does not match dimension of basis.")
        end
        return new{B,T}(b, vec(data))
    end
end
function Ket{B}(b::B, data::AbstractArray{<:Number}) where {B<:Basis}
    return Ket{B,ComplexF64}(b, vec(data))
end
function Ket{T}(b::Basis, data::AbstractArray{<:Number}) where {T<:Complex{<:AbstractFloat}}
    return Ket{typeof(b),T}(b, vec(data))
end
function Ket(b::B, data::AbstractArray{<:Number}) where {B<:Basis}
    return Ket{B,ComplexF64}(b, vec(data))
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
    basis::B
    data::Vector{T}
    function Bra{B,T}(
        b::B, data::AbstractArray{<:Number}
    ) where {B<:Basis,T<:Complex{<:AbstractFloat}}
        if length(b) !== length(data)
            error("Dimension of data does not match dimension of basis.")
        end
        return new{B,T}(b, vec(data))
    end
end
function Bra{B}(b::B, data::AbstractArray{<:Number}) where {B<:Basis}
    return Bra{B,ComplexF64}(b, vec(data))
end
function Bra{T}(b::Basis, data::AbstractArray{<:Number}) where {T<:Complex{<:AbstractFloat}}
    return Bra{typeof(b),T}(b, vec(data))
end
function Bra(b::B, data::AbstractArray{<:Number}) where {B<:Basis}
    return Bra{B,ComplexF64}(b, vec(data))
end
Base.:(==)(u::Bra, v::Bra) = u.data == v.data && u.basis == v.basis

# Arithmetic and vector operations, addition, subtraction, scalar multiplication and
# division, and inner product.
Base.:(+)(u::Ket{B}, v::Ket{B}) where {B<:Basis} = Ket(u.basis, u.data + v.data)
Base.:(+)(::Ket, ::Ket) = error("The bases of the two kets do not match.")
Base.:(-)(u::Ket{B}, v::Ket{B}) where {B<:Basis} = Ket(u.basis, u.data - v.data)
Base.:(-)(::Ket, ::Ket) = error("The bases of the two kets do not match.")
Base.:(*)(x::Number, v::Ket) = Ket(v.basis, x * v.data)
Base.:(*)(v::Ket, x::Number) = Ket(v.basis, v.data * x)
Base.:(/)(v::Ket, x::Number) = Ket(v.basis, v.data / x)

Base.:(+)(u::Bra{B}, v::Bra{B}) where {B<:Basis} = Bra(u.basis, u.data + v.data)
Base.:(+)(::Bra, ::Bra) = error("The bases of the two kets do not match.")
Base.:(-)(u::Bra{B}, v::Bra{B}) where {B<:Basis} = Bra(u.basis, u.data - v.data)
Base.:(-)(::Bra, ::Bra) = error("The bases of the two kets do not match.")
Base.:(*)(x::Number, v::Bra) = Bra(v.basis, x * v.data)
Base.:(*)(v::Bra, x::Number) = Bra(v.basis, v.data * x)
Base.:(/)(v::Bra, x::Number) = Bra(v.basis, v.data / x)

Base.:(*)(u::Bra{B}, v::Ket{B}) where {B<:Basis} = dot(u.data, v.data)
Base.:(*)(::Bra, ::Ket) = error("The bases of the two kets do not match.")

"""
    copy(ψ::StateVector)
"""
Base.copy(ψ::T) where {T<:StateVector} = T(ψ.basis, copy(ψ.data))

"""
    dagger(ψ::StateVector)

Hermitian transpose of the state vector. Using `dagger` on a `Ket` or a `Bra` will
respectively return a `Bra` or a `Ket` in the same basis as the input state, but whose data
has been complex conjugated.
"""
dagger(ψ::Ket) = Bra(ψ.basis, ψ.data')
dagger(ψ::Bra) = Ket(ψ.basis, ψ.data')
Base.adjoint(ψ::StateVector) = dagger(ψ)

"""
    norm(ψ::StateVector)

Compute the norm of the state vector `ψ`.
"""
LinearAlgebra.norm(ψ::StateVector) = norm(ψ.data)

"""
    normalize(ψ::StateVector)

Normalizes the state vector `ψ`.
"""
LinearAlgebra.normalize(ψ::StateVector) = ψ / norm(ψ)

"""
    normalize!(ψ::StateVector)

In-place normalization of the state vector `ψ`.
"""
LinearAlgebra.normalize!(ψ::StateVector) = normalize!(ψ.data)