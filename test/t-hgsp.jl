# t-HGSP operations

@testitem "hgsym / hgantisym" setup = [Reference] begin
    for sz in [(3, 3, 3), (3, 3, 2, 4), (2, 2, 3, 3, 2)]
        M = length(sz)
        B = randn(sz...)
        A = B + permutedims(B, (2, 1, 3:M...))    # symmetric along modes 1 and 2

        As = hgsym(A)
        @test size(As) == (sz[1], sz[2], (2 .* sz[3:end] .+ 1)...)
        @test As ≈ Reference.symmetric_jl(A)
        @test tadjoint(As) == As
        @test hgantisym(As) ≈ A
        @test hgantisym(As) ≈ Reference.anti_symmetric_jl(As)
    end

    # Floating-point output for integer input
    A = ones(Int, 2, 2, 2)
    @test_broken try
        eltype(hgsym(A)) == Float64 && hgsym(A) ≈ hgsym(float(A))
    catch
        false
    end

    # Input must be symmetric along modes 1 and 2
    @test_broken try
        hgsym(randn(3, 3, 3))
        false
    catch err
        err isa ArgumentError
    end
end

@testitem "HGSignal" setup = [Reference] begin
    x = randn(4)
    for M in 2:5
        X = HGSignal(x, M)
        @test size(X) == (4, 1, ntuple(_ -> 4, M - 2)...)
        @test Array(X) ≈ Reference.HypergraphSignal_jl(reshape(x, :, 1), M)
    end

    # Bounds checking
    @test_broken try
        HGSignal(x, 3)[1, 5, 1]
        false
    catch err
        err isa BoundsError
    end
end

@testitem "hgsignalfit" begin
    # Exact hypergraph signals are recovered up to sign
    y = randn(4)
    @test hgsignalfit(reshape(y, 4, 1)) == y
    y = randn(4)
    ŷ = @test_logs (:warn, "M = 3 solution has ambiguous sign") hgsignalfit(
        Array(HGSignal(y, 3)),
    )
    @test ŷ ≈ y || ŷ ≈ -y

    @test_throws ErrorException hgsignalfit(randn(3, 2, 3))
    @test_throws ErrorException hgsignalfit(Array(HGSignal(randn(3), 4)))
end

@testitem "hgshift / hgfilter" setup = [Reference] begin
    edges = Set([Set([1, 2, 3]), Set([2, 4]), Set([1, 4]), Set([3, 4])])
    for M in 3:4
        H = Hypergraph{M}(edges, 4)
        A = Array(AdjacencyTensor(H))
        x = randn(4)
        @test hgshift(A, x) ≈ Reference.t_mult_jl(Reference.symmetric_jl(A), x)
    end

    A = Array(AdjacencyTensor(Hypergraph(edges)))
    x = randn(4)
    y = @test_logs (:warn,) hgfilter(A, x)
    Ys = Reference.t_mult_jl(Reference.symmetric_jl(A), x)
    ŷ = @test_logs (:warn,) hgsignalfit(Reference.anti_symmetric_jl(Ys))
    @test y ≈ ŷ || y ≈ -ŷ    # roundoff differences can flip the eigenvector sign
    @test (@test_logs (:warn,) hgfilter(A, -x)) ≈ y    # shift is even in x
end