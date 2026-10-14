## t-hgsp operations

"""
    hgsym(A::Array{T,M}) where {T<:Real,M}

t-HGSP symmetrization of dense tensor `A`.
"""
function hgsym(A::Array{T,M}) where {T<:Real,M}
    sz = size(A)

    # Initialize array with zeros
    As = zeros(T, ntuple(k -> k <= 2 ? sz[k] : 2*sz[k]+1, Val(M)))

    # Copy in blocks and their reversed versions one at a time
    for I in CartesianIndices(ntuple(k -> k <= 2 ? 1 : 2, Val(M)))
        dest_inds = ntuple(Val(M)) do k
            if k <= 2
                return sz[k]
            elseif I[k] == 1    # copy in same order
                return 2:(sz[k]+1)
            else                # copy in reverse order
                return (2*sz[k]+1):-1:(sz[k]+2)
            end
        end
        copyto!(As, CartesianIndices(dest_inds), A, CartesianIndices(A))
    end

    # Divide two for each reflection
    As ./= 2^(M-2)

    return As
end

# TODO: add note to devdocs about trying copying using `views` and finding it was slower
# src_inds = ntuple(k -> k <= 2 ? axes(A, k) : [axes(A, k); reverse(axes(A, k))], Val(M))
# dest_inds = ntuple(k -> k <= 2 ? axes(A, k) : (2:size(As, k)), Val(M))
# copyto!(view(As, dest_inds...), view(A, src_inds...))
# As[dest_inds...] .= view(A, src_inds...)
# As[dest_inds...] = view(A, src_inds...)

"""
    hgantisym(A::Array{T,M}) where {T<:Real,M}

t-HGSP anti-symmetrization of dense tensor `A`.
"""
function hgantisym(As::Array{T,M}) where {T<:Real,M}
    sz = ntuple(k -> k <= 2 ? size(As, k) : convert(Int, (size(As, k)-1)/2), Val(M))

    # Initialize array with zeros
    A = zeros(T, sz)

    # Accumulate blocks and their reversed versions one at a time
    for I in CartesianIndices(ntuple(k -> k <= 2 ? 1 : 2, Val(M)))
        src_inds = ntuple(Val(M)) do k
            if k <= 2
                return sz[k]
            elseif I[k] == 1    # copy in same order
                return 2:(sz[k]+1)
            else                # copy in reverse order
                return (2*sz[k]+1):-1:(sz[k]+2)
            end
        end
        A .+= view(As, CartesianIndices(src_inds))
    end

    return A
end

"""
    HGSignal{T,M}(signal) <: AbstractArray{T,M}

Hypergraph signal tensor corresponding to the signal vector `signal`.
The tensor is an `M`-way `N × 1 × N × ⋯ × N` tensor where `N` is the length of `signal`.
The array is not materialized in memory.

# Fields
- `signal::Vector{T}`: signal vector
"""
struct HGSignal{T<:Real,M} <: AbstractArray{T,M}
    signal::Vector{T}
end
HGSignal(signal::Vector{T}, M) where {T<:Real} = HGSignal{T,M}(signal)
Base.size(X::HGSignal{T,M}) where {T,M} = ntuple(k -> k == 2 ? 1 : length(X.signal), Val(M))
Base.IndexStyle(::Type{<:HGSignal}) = IndexCartesian()
Base.getindex(X::HGSignal{T,M}, I::Vararg{Int,M}) where {T,M} =
    prod(X.signal[I[k]] for k in 1:M if k != 2)

# TODO: more efficient Array method using tensor kernels (if needed)

"""
    hgsignalfit(X::Array{T,M}) where {T<:Real,M}

Fit a hypergraph signal to the tensor `X`.
"""
function hgsignalfit(X::Array{T,M}) where {T<:Real,M}
    N = size(X, 1)

    # Check dims
    size(X, 2) != 1 && error("size(X,2) must equal one")
    any(!=(N), size(X)[3:end]) && error("size(X)[3:end] must all match size(X,1)")

    # Handle different cases
    if M == 1        # just a signal vector
        return X
    elseif M == 2    # just a column matrix
        return vec(X)
    elseif M == 3    # reshaped matrix solved by eigen
        @warn "M = 3 solution has ambiguous sign"
        Xr = dropdims(X; dims = 2)
        λ, V = eigen(Symmetric((Xr + Xr')/2))
        return sqrt(λ[end])*V[:, end]
    else             # M > 3 becomes a tensor eigenvalue problem
        error("M > 3 not supported (yet!)")
    end
end

"""
    hgshift(A::Array{T,M}, x::Vector{T}) where {T<:Real,M}

Perform a t-HGSP shift of the `x` vector using the shift `A` tensor,
returning the output tensor.
"""
function hgshift(A::Array{T,M}, x::Vector{T}) where {T<:Real,M}
    # Expand signal into full tensor
    X = Array(HGSignal(x, M))

    # Symmetrize
    As = hgsym(A)
    Xs = hgsym(X)

    # Return t-product
    return tmul(As, Xs)
end

"""
    hgfilter(A::Array{T,M}, x::Vector{T}) where {T<:Real,M}

Perform a t-HGSP shift of the `x` vector using the shift `A` tensor,
then fit a hypergraph signal to the output tensor.
"""
function hgfilter(A::Array{T,M}, x::Vector{T}) where {T<:Real,M}
    # Perform t-HGSP shift
    Ys = hgshift(A, x)

    # Return best-fit hypergraph signal
    return hgsignalfit(hgantisym(Ys))
end
