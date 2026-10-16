## t-algebra operations

"""
    tmul(A::Array{T,M}, B::Array{T,M}) where {T<:Real,M}

t-product of dense tensors `A` and `B` implemented via FFT.
"""
function tmul(A::Array{T,M}, B::Array{T,M}) where {T<:Real,M}
    # Check dims
    M >= 2 || error("number of dims must be at least two")
    size(A, 2) == size(B, 1) || error("size(A,2) must match size(B,1)")
    size(A)[3:M] == size(B)[3:M] || error("size(A)[3:end] must match size(B)[3:end]")

    # Compute FFTs
    Ah = fft(A, 3:M)
    Bh = fft(B, 3:M)

    # Compute slice-wise products
    Th = promote_type(eltype(Ah), eltype(Bh))
    Ch = similar(Bh, Th, size(Ah, 1), size(Bh, 2), size(Ah)[3:M]...)
    for I in CartesianIndices(size(A)[3:M])
        mul!(view(Ch, :, :, I), view(Ah, :, :, I), view(Bh, :, :, I))
    end

    # Compute IFFT and truncate to real
    C = real(ifft(Ch, 3:M))

    return C
end

"""
    ⋆(A, B)

t-product of tensors `A` and `B` (an alias for `tmul`). Typed via `\\star<tab>`.
Infix `A⋆B⋆C` is left-associative, i.e., it computes `(A⋆B)⋆C`.
"""
const ⋆ = tmul

# TODO: time against python/matlab implementations
# TODO: experiment more with pre-computing FFTW plans
# TODO: add note to devdocs about trying stack with eachslice and finding it was slower

"""
    tadjoint(A::Array{T,M}) where {T<:Real,M}

t-algebra transpose of dense tensor `A`.
"""
function tadjoint(A::Array{T,M}) where {T<:Real,M}
    # Transpose each frontal slice
    P = PermutedDimsArray(A, (2, 1, 3:M...))

    # Reverse the transposed slices "centered" at the first slice
    # (corresponds to circshift + reverse)
    R = circshift(P, ntuple(k -> k <= 2 ? 0 : -1, Val(M)))
    reverse!(R; dims = (3:M...,))

    return R
end

# TODO: make this lazy instead like the usual transpose since we probably mostly
#       need this for doing tmul, for which we could define the overload
# TODO: note in devdocs that a nontrivial amount of the time is spent dealing
#       with the permutation along the first two modes. this could also be
#       addressed by returning a PermutedDimsArray at the end, but that just
#       defers that computation - we'd likely hit it soon after in an fft call.

"""
    teigen(A::Array{T,M}) where {T<:Real,M}

t-eigendecomposition of dense tensor `A` implemented via FFT.
"""
function teigen(A::Array{T,M}; eigenargs...) where {T<:Real,M}
    # Check dims
    size(A, 1) == size(A, 2) || error("size(A,1) must match size(A,2)")

    # Compute FFT
    Ah = fft(A, 3:M)

    # Compute and stack slice-wise eigendecompositions
    Eh = eigen!.(eachslice(Ah; dims = (3:M...,)); eigenargs...)
    Λh = stack(Diagonal.(getproperty.(Eh, :values)))
    Vh = stack(getproperty.(Eh, :vectors))

    # Compute IFFT
    Λ = ifft(Λh, 3:M)
    V = ifft(Vh, 3:M)

    return Λ, V
end

# TODO: note that t-eigendecomposition does not necessarily have real eigvals
#       when the adjacency tensor is not also symmetric along modes 1 and 2
