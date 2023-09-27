abstract type QuantumState end
abstract type StateVector <: QuantumState end

"""
    Ket(b::Basis, data::AbstractVector{<:Number})

Create a ket in the basis `b` with vector representation given by `data`.

If the elements of `data` are not `ComplexF64` they will beconverted to this type.

# Examples
```jldoctest
julia> Ket(FockBasis(4), [1, 0, 0, 0, 0])
Ket{FockBasis{Int64}}(FockBasis{Int64}((5,), 4, 0), ComplexF64[1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im])

julia> Ket(SpinBasis(1 // 2), [1/sqrt(2), im/sqrt(2)])
Ket{SpinBasis{1//2, Int64}}(SpinBasis{1//2, Int64}((2,), 1//2), ComplexF64[0.7071067811865475 + 0.0im, 0.0 + 0.7071067811865475im])
```
"""
mutable struct Ket{B} <: StateVector
    basis::B
    data::Vector{ComplexF64}
    function Ket{B}(b::B, data::AbstractVector{<:Number}) where {B<:Basis}
        if length(b) !== length(data)
            error("Dimension of data does not match dimension of basis.")
        end
        return new{B}(b, Vector{ComplexF64}(data))
    end
end
function Ket(b::B, data::AbstractVector{<:Number}) where {B<:Basis}
    return Ket{B}(b, data)
end
Base.:(==)(u::Ket, v::Ket) = u.data == v.data && u.basis == v.basis

"""
    Ket(b::Basis, data::Adjoint{<:Number,<:AbstractVector{<:Number}})

Create ket bra in the basis `b` with vector representation given by the adjointed of `data`.

If the elements of `data` are not `ComplexF64` they will beconverted to this type.

Since `data` is given as an adjointed vector, it will be adjointed again before when
constructing the `Ket`, such that the data field will be a vector.
"""
function Ket{B}(b::B, data::Adjoint{<:Number,<:AbstractVector{<:Number}}) where {B<:Basis}
    return Ket{B}(b, complex(data)')
end
function Ket(b::Basis, data::Adjoint{<:Number,<:AbstractVector{<:Number}})
    return Ket(b, complex(data)')
end

"""
    Bra(b::Basis, data::AbstractVector{<:Number})

Create a bra in the basis `b` with vector representation given by `data`.

If the elements of `data` are not `ComplexF64` they will beconverted to this type.

Furthermore, the `data` field of a `Bra` is the adjoint of the `data` argiment

# Examples
```jldoctest
julia> Bra(FockBasis(4), [1, 0, 0, 0, 0])
Bra{FockBasis{Int64}}(FockBasis{Int64}((5,), 4, 0), ComplexF64[1.0 + 0.0im 0.0 + 0.0im … 0.0 + 0.0im 0.0 + 0.0im])

julia> Bra(SpinBasis(1 // 2), [1/sqrt(2), im/sqrt(2)])
Bra{SpinBasis{1//2, Int64}}(SpinBasis{1//2, Int64}((2,), 1//2), ComplexF64[0.7071067811865475 + 0.0im 0.0 + 0.7071067811865475im])
```
"""
mutable struct Bra{B} <: StateVector
    basis::B
    data::Adjoint{ComplexF64,Vector{ComplexF64}}
    function Bra{B}(b::B, data::AbstractVector{<:Number}) where {B<:Basis}
        if length(b) !== length(data)
            error("Dimension of data does not match dimension of basis.")
        end
        return new{B}(b, Vector{ComplexF64}(data)')
    end
end
function Bra(b::B, data::AbstractVector{<:Number}) where {B<:Basis}
    return Bra{B}(b, data)
end
Base.:(==)(u::Bra, v::Bra) = u.data == v.data && u.basis == v.basis

"""
    Bra(b::Basis, data::Adjoint{<:Number,<:AbstractVector{<:Number}})

Create a bra in the basis `b` with vector representation given by `data`.

If the elements of `data` are not `ComplexF64` they will beconverted to this type.

Here `data` is already adjointed, and so will not be adjointed agian when constructing the
`Bra`.
"""
function Bra{B}(b::B, data::Adjoint{<:Number,<:AbstractVector{<:Number}}) where {B<:Basis}
    return Bra{B}(b, complex(data)')
end
function Bra(b::Basis, data::Adjoint{<:Number,<:AbstractVector{<:Number}})
    return Bra(b, complex(data)')
end

# Arithmetic and vector operations, addition, subtraction, scalar multiplication, scalar
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
Base.:(*)(x::Number, v::Bra) = Bra(v.basis, x' * v.data)
Base.:(*)(v::Bra, x::Number) = Bra(v.basis, v.data * x')
Base.:(/)(v::Bra, x::Number) = Bra(v.basis, v.data / x')

Base.:(*)(u::Bra{B}, v::Ket{B}) where {B<:Basis} = u.data * v.data
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