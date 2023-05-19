abstract type Basis end

"""
    length(b::Basis)

Return the total dimension of the Hilbert space spanned by the basis `b`.
"""
Base.length(b::Basis) = prod(b.dimension)

"""
    GenericBasis(d)

Create a generic basis with dimension `d` for a Hilbert space.

# Examples
```jldoctest
julia> GenericBasis(4)
GenericBasis{Int64}(4)
```
"""
struct GenericBasis{T} <: Basis
    dimension::T
    function GenericBasis(d::T) where {T<:Integer}
        if d > 0
            return new{T}(d)
        else
            error("Basis must have positive dimensions.")
        end
    end
end
GenericBasis(d::T) where {T<:Number} = GenericBasis(convert(Integer, d))
Base.:(==)(b1::GenericBasis, b2::GenericBasis) = b1.dimension == b2.dimension

"""
    CompositeBasis(b::Basis...)

Create a basis for the composite system consisting of the systems corresponding to the given
bases.

# Examples
```jldoctest
julia> CompositeBasis(GenericBasis(8), FockBasis(50, 20))
CompositeBasis{Tuple{GenericBasis{Int64}, FockBasis{Int64}}, Vector{Int64}}((GenericBasis{In
t64}(8), FockBasis{Int64}(31, 50, 20)), [8, 31])
```
"""
struct CompositeBasis{B,T} <: Basis
    bases::B
    dimension::T
    function CompositeBasis(bases::B) where {B<:Tuple{Vararg{Basis}}}
        return new{B, Vector{Int}}(bases, [length(b) for b ∈ bases])
    end
end
CompositeBasis(b::Basis...) = CompositeBasis((b...,))

"""
    FockBasis(cutoff::Integer, offset::Integer=0)

Create a basis for the Fock space starting at the offset and ending at the cutoff.

# Examples
```jldoctest
julia> FockBasis(10)
FockBasis{Int64}(11, 10, 0)

julia> FockBasis(10, 2)
FockBasis{Int64}(9, 10, 2)
```
"""
struct FockBasis{T} <: Basis
    dimension::T
    cutoff::T
    offset::T
    function FockBasis(cutoff::T, offset::T=0) where {T<:Integer}
        if cutoff < 0 || offset < 0 || cutoff < offset
            error("Cufoff must be larger than or equal to offset, and they must both be \
                   positive.")
        end
        return new{T}(cutoff - offset + 1, cutoff, offset)
    end
end
function FockBasis(cutoff::S, offset::T=0) where {S<:Number,T<:Number}
    FockBasis(convert(Integer, cutoff), convert(Integer, offset))
end
Base.:(==)(b1::FockBasis, b2::FockBasis) = b1.cutoff == b2.cutoff && b1.offset == b2.offset

"""
    SpinBasis(spin)

Create a basis for a spin system with total spin equal to `spin`.

The argument `spin` must be integer or half integer. It is then stored as a Rational.

# Examples
```jldoctest
julia> SpinBasis(1 // 2)
SpinBasis{1//2, Int64}(2, 1//2)

julia> SpinBasis(100)
SpinBasis{100//1, Int64}(201, 100//1)
```
"""
struct SpinBasis{S,T} <: Basis
    dimension::T
    spin::Rational{T}
    function SpinBasis(spin::Rational{T}) where {T<:Integer}
        num = numerator(spin)
        den = denominator(spin)
        if !(den == 1 || den == 2)
            error("The spin number must be either integer or half-integer.")
        end
        if num ≤ 0
            error("The spin number must be positive.")
        end
        return new{spin,T}(numerator(2 * spin + 1), spin)
    end
end
SpinBasis(spin::T) where {T<:Number} = SpinBasis(convert(Rational, spin))
Base.:(==)(b1::SpinBasis, b2::SpinBasis) = b1.spin == b2.spin
