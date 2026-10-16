# Hypergraph types

@testitem "Hypergraph construction" begin
    edges = Set([Set([1, 2, 3]), Set([2, 4]), Set([1, 4])])

    # Inferred nvertices and max cardinality
    H = Hypergraph(edges)
    @test H isa Hypergraph{3}
    @test H.nvertices == 4
    @test H.hyperedges == edges

    # Explicit nvertices and max cardinality
    @test Hypergraph(edges, 6).nvertices == 6
    @test Hypergraph{4}(edges, 4) isa Hypergraph{4}

    # Validation
    @test_throws ArgumentError Hypergraph(edges, 3)
    @test_throws ArgumentError Hypergraph{2}(edges, 4)

    # Non-Set input
    @test Hypergraph{3}([[1, 2], [2, 3]], 3).hyperedges == Set([Set([1, 2]), Set([2, 3])])
    @test Hypergraph([(1, 2), (2, 3, 4)]) isa Hypergraph{3}
    @test Hypergraph([[1, 2, 2], [2, 3]]) isa Hypergraph{2}    # repeated vertices

    # Stored hyperedges do not alias the input
    e = Set([1, 2])
    S = Set([e])
    H = Hypergraph(S)
    push!(S, Set([7, 8, 9]))
    push!(e, 3)
    @test H.hyperedges == Set([Set([1, 2])])

    # Empty hyperedge set (max cardinality cannot be inferred)
    @test Hypergraph{2}(Set{Set{Int}}(), 3).nvertices == 3
    @test_throws "Hypergraph{M}" Hypergraph(Set{Set{Int}}(), 3)
    @test_throws ArgumentError Hypergraph(Set{Set{Int}}())
end

@testitem "AdjacencyTensor" setup = [Reference] begin
    for (edges, n) in [
        (Set([Set([1, 2, 3]), Set([2, 4]), Set([1, 4])]), 4),
        (Set([Set([1, 2, 3, 4]), Set([2, 5]), Set([1, 3, 5]), Set([4])]), 5),
    ]
        H = Hypergraph(edges, n)
        M = maximum(length, edges)
        A = AdjacencyTensor(H)

        @test size(A) == ntuple(_ -> n, M)
        Hm = Reference.hyperedges_to_incidence_jl(collect(edges); num_nodes = n)
        @test Array(A) ≈ Reference.adjacency(Hm, n, M, length(edges))

        # Symmetric along every pair of modes
        Ad = Array(A)
        @test permutedims(Ad, (2, 1, 3:M...)) == Ad
        @test permutedims(Ad, (M, 2:(M-1)..., 1)) == Ad
    end

    # Bounds checking
    A = AdjacencyTensor(Hypergraph(Set([Set([1, 2, 3])])))
    @test_throws BoundsError A[0, 0, 9]
    @test_throws BoundsError A[1, 1, 4]
end
