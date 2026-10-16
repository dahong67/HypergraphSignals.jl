module HypergraphSignals

# Imports
using Combinatorics: stirlings2
using Compat: @compat
using FFTW: fft, ifft
using LinearAlgebra: Diagonal, Symmetric, eigen, eigen!, mul!

# Hypergraph types
include("hypergraphs.jl")
export AdjacencyTensor, Hypergraph

# t-algebra operations
include("t-algebra.jl")
export tadjoint, teigen, tmul
@compat public (⋆)

# t-HGSP operations
include("t-hgsp.jl")
export HGSignal, hgantisym, hgfilter, hgshift, hgsignalfit, hgsym

end