## Hypergraph types/functions

"""
    Hypergraph{M}(hyperedges, nvertices)

Form a hypergraph stored as a `Set` of hyperedges,
each of which is a `Set` of at-most `M` vertex indices
from `1:nvertices`.
The hyperedges can be any collection of collections of vertex indices,
and are copied into a new `Set{Set{Int}}`.

# Fields
- `hyperedges::Set{Set{Int}}`: set of hyperedges
- `nvertices::Int`: number of vertices in the hypergraph
"""
struct Hypergraph{M}
    hyperedges::Set{Set{Int}}
    nvertices::Int

    function Hypergraph{M}(hyperedges, nvertices) where {M}
        # Copy into a new set of sets
        hyperedges = Set{Set{Int}}(Set{Int}(hyperedge) for hyperedge in hyperedges)

        # Check consistency with nvertices and M
        for hyperedge in hyperedges
            any(!in(1:nvertices), hyperedge) && throw(
                ArgumentError("Hyperedge $hyperedge has vertices outside of 1:$nvertices."),
            )
            length(hyperedge) > M && throw(
                ArgumentError("Hyperedge $hyperedge has cardinality greater than $M."),
            )
        end
        return new{M}(hyperedges, nvertices)
    end
end

"""
    Hypergraph(hyperedges, nvertices = maximum(maximum, hyperedges))

Form a hypergraph from a collection of hyperedges (each a collection of vertex indices),
where `nvertices` and the maximum cardinality are inferred.
Use `Hypergraph{M}(hyperedges, nvertices)` when `hyperedges` is empty.
"""
function Hypergraph(hyperedges, nvertices = maximum(maximum, hyperedges; init = 0))
    isempty(hyperedges) && throw(
        ArgumentError(
            "Cannot infer the maximum cardinality from an empty collection of hyperedges. " *
            "Use `Hypergraph{M}(hyperedges, nvertices)` instead.",
        ),
    )
    M = maximum(hyperedge -> length(Set(hyperedge)), hyperedges)
    return Hypergraph{M}(hyperedges, nvertices)
end

# TODO: pretty-printing for Hypergraph
# TODO: basic operations (max cardinality, nvertices, nedges, addedge!, etc)

"""
    AdjacencyTensor{M}(hypergraph) <: AbstractArray{Float64,M}

Adjacency tensor corresponding to the hypergraph `hypergraph`.
The array is not materialized in memory.

# Fields
- `hypergraph::Hypergraph{M}`: hypergraph
"""
struct AdjacencyTensor{M} <: AbstractArray{Float64,M}
    hypergraph::Hypergraph{M}
end
Base.size(A::AdjacencyTensor{M}) where {M} = ntuple(_ -> A.hypergraph.nvertices, Val(M))
Base.IndexStyle(::Type{<:AdjacencyTensor}) = IndexCartesian()
function Base.getindex(A::AdjacencyTensor{M}, I::Vararg{Int,M}) where {M}
    @boundscheck checkbounds(A, I...)
    hyperedges = A.hypergraph.hyperedges
    hyperedge = Set(I)
    if hyperedge in hyperedges
        c = length(hyperedge)
        α = factorial(c)*stirlings2(M, c)
        return c/α
    else
        return 0.0
    end
end

# TODO: more efficient Array method that walks through the edgelist and assigns to entries (if needed)
#       see the permutations and partitions functions in Combinatorics
