# t-algebra operations

@testitem "tmul" begin
    using FFTW: fft

    for sz in [(2, 3, 4), (2, 3, 3, 4), (2, 3, 3, 2, 5)]
        M = length(sz)
        A = randn(sz...)
        B = randn(sz[2], 4, sz[3:end]...)
        C = randn(4, 2, sz[3:end]...)

        # Fourier-domain slices are matrix products
        AB = tmul(A, B)
        @test size(AB) == (sz[1], 4, sz[3:end]...)
        ABh, Ah, Bh = fft(AB, 3:M), fft(A, 3:M), fft(B, 3:M)
        @test all(
            ABh[:, :, I] ≈ Ah[:, :, I] * Bh[:, :, I] for I in CartesianIndices(sz[3:end])
        )

        # Associativity
        @test tmul(tmul(A, B), C) ≈ tmul(A, tmul(B, C))

        # Identity tensor
        Id = zeros(sz[2], sz[2], sz[3:end]...)
        Id[:, :, ntuple(_->1, M-2)...] .= [i == j for i in 1:sz[2], j in 1:sz[2]]
        @test tmul(A, Id) ≈ A
    end

    # Dimension checks
    @test_throws DimensionMismatch tmul(randn(2, 3, 4), randn(2, 3, 4))
    @test_throws DimensionMismatch tmul(randn(2, 3, 4), randn(3, 3, 5))
    @test_throws ArgumentError tmul(randn(3), randn(3))

    # Matrices (M = 2)
    A, B = randn(2, 3), randn(3, 4)
    @test_broken tmul(A, B) ≈ A * B

    # Mixed and complex element types
    @test_broken try
        tmul(ones(Int, 2, 2, 3), ones(2, 2, 3)) ≈ tmul(ones(2, 2, 3), ones(2, 2, 3))
    catch
        false
    end
    @test_broken try
        Az, Bz = randn(ComplexF64, 2, 3, 4), randn(ComplexF64, 3, 2, 4)
        Ch = fft(tmul(Az, Bz), 3)
        Ah, Bh = fft(Az, 3), fft(Bz, 3)
        all(Ch[:, :, k] ≈ Ah[:, :, k] * Bh[:, :, k] for k in 1:4)
    catch
        false
    end
end

@testitem "⋆" begin
    using HypergraphSignals: ⋆

    @test (⋆) === tmul
    @test isconst(HypergraphSignals, :⋆)
    if VERSION >= v"1.11-"  # public only declared on Julia 1.11+
        @test Base.ispublic(HypergraphSignals, :⋆)
    end

    A, B, C = randn(2, 3, 4), randn(3, 3, 4), randn(3, 2, 4)
    @test A ⋆ B ⋆ C ≈ tmul(tmul(A, B), C)
end

@testitem "tadjoint" setup = [Reference] begin
    for sz in [(2, 3), (2, 3, 4), (3, 3, 4), (2, 3, 3, 4), (3, 3, 4, 4), (2, 3, 3, 3, 5)]
        A = randn(sz...)
        B = randn(sz[2], sz[1], sz[3:end]...)

        @test tadjoint(A) ≈ Reference.t_tran_jl(A)
        @test tadjoint(tadjoint(A)) == A
        if length(sz) >= 3
            @test tadjoint(tmul(A, B)) ≈ tmul(tadjoint(B), tadjoint(A))
        end
    end

    # Complex element types
    @test_broken try
        Az = randn(ComplexF64, 2, 3, 4)
        tadjoint(Az) ≈ Reference.t_tran_jl(Az)
    catch
        false
    end
end

@testitem "teigen" setup = [Reference] begin
    using FFTW: fft
    using LinearAlgebra: diag

    # t-self-adjoint tensors (symmetrized adjacency tensors)
    for edges in [
        Set([Set([1, 2, 3]), Set([2, 4]), Set([1, 4])]),
        Set([Set([1, 2, 3, 4]), Set([2, 5]), Set([1, 3, 5]), Set([4])]),
    ]
        A = hgsym(Array(AdjacencyTensor(Hypergraph(edges))))
        M = ndims(A)
        Λ, V = teigen(A)

        # Fourier-domain slices are eigendecompositions
        Ah, Λh, Vh = fft(A, 3:M), fft(Λ, 3:M), fft(V, 3:M)
        @test all(
            Ah[:, :, I] * Vh[:, :, I] ≈ Vh[:, :, I] * Λh[:, :, I] for
            I in CartesianIndices(size(A)[3:end])
        )

        # Eigenvalues match the reference implementation
        S = Reference.t_eig_jl(A; full_output = true)[6]
        @test all(
            sort(real(diag(Λh[:, :, I]))) ≈ sort(real(diag(S[:, :, I]))) for
            I in CartesianIndices(size(A)[3:end])
        )

        # Real eigenvalues, sorted consistently across slices
        @test_broken eltype(Λ) <: Real
    end

    # Frontal slices must be square
    @test_throws DimensionMismatch teigen(randn(2, 3, 4))

    # Non-t-self-adjoint tensors are not supported
    @test_broken try
        teigen(randn(3, 3, 4))
        false
    catch err
        err isa ArgumentError
    end

    # Matrices (M = 2)
    @test_broken try
        A = randn(3, 3)
        A = A + A'
        Λ, V = teigen(A)
        A * V ≈ V * Λ
    catch
        false
    end
end
