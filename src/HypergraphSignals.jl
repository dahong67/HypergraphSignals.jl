module HypergraphSignals

using Combinatorics: stirlings2
using LinearAlgebra: eigen!, mul!, Diagonal, Symmetric
export AdjacencyTensor, Hypergraph
export tmul, tadjoint, teigen
export hgsym, hgantisym, HGSignal, hgsignalfit, hgshift, hgfilter
public ⋆

include("hypergraphs.jl")
include("t-algebra.jl")
include("t-hgsp.jl")

end