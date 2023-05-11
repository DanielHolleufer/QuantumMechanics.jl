abstract type Basis end

Base.length(b::Basis) = prod(b.dimensions)

struct GenericBasis <: Basis
    dimensions::Vector{Int}
end
GenericBasis(N::Integer) = GenericBasis([N])


struct FockBasis <: Basis
    dimensions::Vector{Int}
    cutoff::Int
    offset::Int
    function FockBasis(cutoff, offset = 0)
        new([cutoff + 1 - offset], cutoff, offset)
    end
end
