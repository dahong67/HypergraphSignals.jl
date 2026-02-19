module HypergraphSignals

using Combinatorics, NPZ
using FFTW, LinearAlgebra

function hyperedges_to_incidence_jl(hyperedges; num_nodes = nothing, one_based = true)
    if isempty(hyperedges)
        throw(ArgumentError("hyperedges cannot be empty"))
    end

    offset = one_based ? 0 : 1
    all_nodes = [node for edge in hyperedges for node in edge]
    max_node_in_data = isempty(all_nodes) ? 0 : maximum(all_nodes)

    n =
        isnothing(num_nodes) ? (max_node_in_data + offset) :
        max(num_nodes, max_node_in_data + offset)
    e = length(hyperedges)

    H = zeros(Int, n, e)

    for (col, edge) in enumerate(hyperedges)
        for node in edge
            idx = one_based ? node : node + 1

            if !(1 <= idx <= n)
                throw(BoundsError(H, (idx, col)))
            end
            H[idx, col] = 1
        end
    end

    return H
end

function load_traffic_hypergraph_jl(file_path, ; M)
    data = npzread(file_path)

    N = Int(data["N"])
    W = data["W"]

    A = (W .> 0) .| (transpose(W) .> 0)

    hyperedge_sets = Set{Set{Int}}()

    for v in 1:N
        neighbors = findall(x -> x != 0, A[v, :])

        size_neighbors = length(neighbors) + 1

        if size_neighbors <= M
            push!(hyperedge_sets, Set(vcat(v, neighbors)))
        else
            for subset in combinations(neighbors, M - 1)
                push!(hyperedge_sets, Set(vcat(v, subset)))
            end
        end
    end
    hyperedges = [sort(collect(edge)) for edge in hyperedge_sets]

    H = hyperedges_to_incidence_jl(hyperedges; num_nodes = N, one_based = true)

    return H
end

function adjacency(H_jl, N_jl, M_jl, E_jl)
    A_jl = zeros(fill(N_jl, M_jl)...)
    for e in 1:E_jl
        location = findall(!iszero, H_jl[:, e])
        all_perms = generate_perms(location, maximum(sum(H_jl; dims = 1)))
        num_perms = length(all_perms)
        c = length(location)
        for perm in all_perms
            A_jl[perm...] = c / num_perms
        end
    end
    return A_jl
end

function generate_perms(lst, M)
    if length(lst) == M
        return collect(permutations(lst))
    else
        c = length(lst)
        M_minus_c = M - c

        complementary_lst = collect(with_replacement_combinations(lst, M_minus_c))

        final_lst = []

        for comp in complementary_lst
            full_perm_base = vcat(comp, lst)
            perms = unique(permutations(full_perm_base))
            append!(final_lst, perms)
        end
        return final_lst
    end
end

function degree(H_jl, N_jl, M_jl)
    D_jl = zeros(fill(N_jl, M_jl)...)
    d_jl = vec(sum(H_jl; dims = 2)) # degree of each nodes

    for i in 1:N_jl
        idx = ntuple(_ -> i, M_jl)
        D_jl[idx...] = d_jl[i]
    end
    return D_jl
end

function laplacian(H_jl, N_jl, M_jl, E_jl)
    return degree(H_jl, N_jl, M_jl) - adjacency(H_jl, N_jl, M_jl, E_jl)
end

function HypergraphSignal_jl(S::Matrix, M::Int)
    N, d = size(S)
    function one_dim_hypergraphsignal_jl(s::Vector, M::Int)
        X = s
        for _ in 1:M-2
            new_shape = (ones(Int, ndims(X))..., N)
            X = X .* reshape(s, new_shape...)
        end
        current_dims = size(X)
        return reshape(X, current_dims[1], 1, current_dims[2:end]...)
    end

    results = [one_dim_hypergraphsignal_jl(S[:, i], M) for i in 1:d]
    return cat(results...; dims = 2)
end

function sym_step(A, i)
    sz = collect(size(A))
    sz[i] = 1
    first_0 = zeros(eltype(A), sz...)

    reverse_A = reverse(A; dims = i)
    As = cat(first_0, A, reverse_A; dims = i)

    return As ./ 2
end

function symmetric_jl(A)
    p = ndims(A)
    for i in 3:p
        A = sym_step(A, i)
    end
    return A
end

function anti_sym_step(As, i)
    total_len = size(As, i)
    N = Int((total_len - 1) / 2)

    colons = Any[Colon() for _ in 1:ndims(As)]
    colons[i] = 2:(N+1)

    A = As[colons...]
    return 2 .* A
end

function anti_symmetric_jl(As)
    p = ndims(As)
    for i in p:-1:3
        As = anti_sym_step(As, i)
    end
    return As
end

function t_mult_jl(A, B)
    p = ndims(A)
    if size(A)[3:end] != size(B)[3:end]
        println("Trailing dimensions must match.")
    end
    num_slices = prod(size(A)[3:end])
    shape_C = (size(A, 1), size(B, 2), size(A)[3:end]...)
    flatten_shape_A = (size(A, 1), size(A, 2), num_slices)
    flatten_shape_B = (size(A, 2), size(B, 2), num_slices)

    A_fft = copy(A)
    B_fft = copy(B)
    for i in 3:p
        A_fft = fft(A_fft, i)
        B_fft = fft(B_fft, i)
    end

    flatten_A = reshape(A_fft, flatten_shape_A)
    flatten_B = reshape(B_fft, flatten_shape_B)

    slices = [flatten_A[:, :, k] * flatten_B[:, :, k] for k in 1:num_slices]
    flatten_C = cat(slices...; dims = 3)

    C = reshape(flatten_C, shape_C...)

    for i in p:-1:3
        C = ifft(C, i)
    end

    # Imaginary check
    if sum(abs2.(imag.(C))) < 1e-10
        C = real.(C)
    end

    return C
end

function t_eig_jl(A; order = "LR", full_output = false)
    if size(A, 1) != size(A, 2)
        error("Frontal slices must be square.")
    end

    p = ndims(A)
    num_slices = prod(size(A)[3:end])
    flatten_shape_D = (size(A, 1), size(A, 2), num_slices)
    shape_ten = size(A)

    D = ComplexF64.(A)
    for i in 3:p
        fft!(D, i)
    end

    flatten_D = D

    n             = flatten_shape_D[1]
    flatten_V_fft = zeros(ComplexF64, flatten_shape_D)
    flatten_S_fft = zeros(ComplexF64, flatten_shape_D)
    flatten_Yt    = zeros(ComplexF64, flatten_shape_D)
    flatten_Yt2   = zeros(ComplexF64, flatten_shape_D)

    for i in 1:num_slices
        decomp = eigen(Hermitian(flatten_D[:, :, i]))
        s = decomp.values
        V = decomp.vectors
        s, V = sorting_values_jl(s, V, order)

        maxEig = maximum(abs.(s)) == 0 ? 1.0 : maximum(abs.(s))
        S_diag = diagm(s)

        flatten_Yt[:, :, i] = V * (abs.(S_diag) ./ maxEig) * V'
        flatten_Yt2[:, :, i] = V * (S_diag ./ maxEig) * V'
        flatten_V_fft[:, :, i] = V
        flatten_S_fft[:, :, i] = S_diag
    end

    S = reshape(flatten_S_fft, shape_ten)
    V = reshape(flatten_V_fft, shape_ten)
    Shift_Operator_Abs_Norm = reshape(flatten_Yt, shape_ten)
    Shift_Operator_Norm = reshape(flatten_Yt2, shape_ten)

    S_fourier = copy(S)
    V_fourier = copy(V)

    for i in p:-1:3
        S = ifft(S, i)
        V = ifft(V, i)
        Shift_Operator_Abs_Norm = ifft(Shift_Operator_Abs_Norm, i)
        Shift_Operator_Norm = ifft(Shift_Operator_Norm, i)
    end

    if sum(abs2.(imag.(S))) < 1e-10
        S = real.(S)
    end
    if sum(abs2.(imag.(V))) < 1e-10
        V = real.(V)
    end

    if full_output
        return V,
        S,
        Shift_Operator_Abs_Norm,
        Shift_Operator_Norm,
        V_fourier,
        S_fourier,
        flatten_D
    else
        return V, S
    end
end

function sorting_values_jl(values, vectors, order)
    if order == "LR"
        idx = sortperm(values; rev = true)
    elseif order == "SR"
        idx = sortperm(values)
    else
        error("order must be 'LR' or 'SR'")
    end
    return values[idx], vectors[:, idx]
end

function t_tran_jl(A)
    function tran(X)
        n1, n2, n3 = size(X)
        Xt = zeros(ComplexF64, n2, n1, n3)
        Xt[:, :, 1] = X[:, :, 1]'
        for i in 2:n3
            Xt[:, :, i] = X[:, :, n3-i+2]'
        end
        return Xt
    end

    if ndims(A) <= 3
        return tran(A)
    end

    sz = size(A)
    sz2 = collect(sz)
    sz2[1], sz2[2] = sz[2], sz[1]
    At = zeros(ComplexF64, Tuple(sz2))
    for i in 1:sz[end]
        At[fill(:, ndims(A) - 1)..., i] = tran(selectdim(A, ndims(A), i))
    end
    return At
end

end
