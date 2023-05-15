abstract type Basis end

"""
    length(b::Basis)

Return the total dimension of the Hilbert space spanned by the given basis.
"""
Base.length(b::Basis) = prod(b.dimensions)

struct GenericBasis <: Basis
    dimensions::Vector{Integer}
    function GenericBasis(N::Integer)
        N > 0 ? new([N]) : error("GenericBasis must have a positive dimension, N > 0.")
    end
end
Base.:(==)(b1::GenericBasis, b2::GenericBasis) = b1.dimensions == b2.dimensions


struct FockBasis <: Basis
    dimensions::Vector{Integer}
    cutoff::Integer
    offset::Integer
    function FockBasis(cutoff::Integer, offset::Integer = 0)
        cutoff - offset ≥ 0 ? new([cutoff + 1 - offset], cutoff, offset) : error("Cufoff must be larger than or equal to offset.")
    end
end
Base.:(==)(b1::FockBasis, b2::FockBasis) = b1.cutoff == b2.cutoff && b1.offset == b2.offset
