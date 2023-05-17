abstract type Basis end

"""
    length(b::Basis)

Return the total dimension of the Hilbert space spanned by the basis `b`.
"""
Base.length(b::Basis) = prod(b.dimensions)

"""
    GenericBasis(v)

Create a generic basis with dimension `v` for a Hilbert space.

The input `v` can be either an integer or a vector of integers. If `v` is an integer, then
the basis will have one component with dimension `v`. If `v` is a vector, then the basis
will have components equal to `length(v)` each with dimension equal to the corresponding
component of `v`.

# Examples
```jldoctest
julia> GenericBasis(4)
GenericBasis{Int64}([4])

julia> GenericBasis([2, 2, 2])
GenericBasis{Vector{Int64}}([2, 2, 2])
```
"""
struct GenericBasis{T} <: Basis
    dimensions::Vector{T}
    function GenericBasis(v::Vector{T}) where {T<:Integer}
        if prod(v .> 0)
            new{T}(v)
        else
            error("Basis must have positive dimensions.")
        end
    end
end
function GenericBasis(M::Matrix{T}) where {T<:Integer}
    if size(M)[1] != 1 && size(M)[2] != 1
        error("Array of dimensions must be either row or coloumn vector/matrix.")
    end
    return GenericBasis(vec(M))
end
GenericBasis(M::Matrix{T}) where {T<:AbstractFloat} = GenericBasis(convert.(Integer, M))
GenericBasis(N::T) where {T<:Integer} = GenericBasis([N])
GenericBasis(N::AbstractFloat) = GenericBasis(convert(Integer, N))
Base.:(==)(b1::GenericBasis, b2::GenericBasis) = b1.dimensions == b2.dimensions

"""
    FockBasis(cutoff::Integer, offset::Integer=0)

Create a basis for the Fock space starting at the offset and ending at the cutoff.

# Examples
```jldoctest
julia> FockBasis(10)
FockBasis{Int64}([11], 10, 0)

julia> FockBasis(10, 2)
FockBasis{Int64}([9], 10, 2)
```
"""
struct FockBasis{T} <: Basis
    dimensions::Vector{T}
    cutoff::T
    offset::T
    function FockBasis(cutoff::T, offset::T=0) where {T<:Integer}
        if cutoff < 0 || offset < 0 || cutoff < offset
            error("Cufoff must be larger than or equal to offset, and they must both be \
                   positive.")
        end
        new{T}([cutoff - offset + 1], cutoff, offset)
    end
end
function FockBasis(cutoff::S, offset::T=0) where {S<:Number,T<:Number}
    FockBasis(convert(Integer, cutoff), convert(Integer, offset))
end
Base.:(==)(b1::FockBasis, b2::FockBasis) = b1.cutoff == b2.cutoff && b1.offset == b2.offset


struct SpinBasis{S} <: Basis
    dimensions::Vector{Integer}
    spin::Rational{Integer}
    function SpinBasis{S}(spin::Rational{<:Integer}) where {S}
        num = numerator(spin)
        den = denominator(spin)
        if !(den == 1 || den == 2)
            error("The spin number must be either integer or half-integer.")
        end
        if num ≤ 0
            error("The spin number must be positive.")
        end
        new{spin}([numerator(2 * spin + 1)], spin)
    end
end
SpinBasis(spin::Rational) = SpinBasis{spin}(spin)
SpinBasis(spin) = SpinBasis(convert(Rational{Int}, spin))
Base.:(==)(b1::SpinBasis, b2::SpinBasis) = b1.spin == b2.spin

