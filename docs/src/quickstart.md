# Quick start guide

Let's install HypergraphSignals
and run our first t-product hypergraph signal processing (t-HGSP) computations!

## Step 1: Install Julia

Go to [https://julialang.org/downloads](https://julialang.org/downloads)
and install the current stable release.

To check your installation,
open up Julia and try a simple calculation like `1+1`:
```@repl
1 + 1
```

More info: [https://docs.julialang.org/en/v1/manual/getting-started/](https://docs.julialang.org/en/v1/manual/getting-started/)

## Step 2: Install HypergraphSignals

HypergraphSignals can be installed using
Julia's excellent builtin package manager.
It is not yet registered,
so we install it directly from its GitHub repository.

```julia-repl
julia> import Pkg; Pkg.add(url = "https://github.com/dahong67/HypergraphSignals.jl")
```

This downloads, installs, and precompiles HypergraphSignals
(and all its dependencies).
Don't worry if it takes a few minutes to complete.

!!! tip "Tip: Interactive package management with the Pkg REPL mode"

    Here we used the [functional API](https://pkgdocs.julialang.org/v1/api/)
    for the builtin package manager.
    Pkg also has a very nice interactive interface (called the Pkg REPL)
    that is built right into the Julia REPL!

    Learn more here: [https://pkgdocs.julialang.org/v1/getting-started/](https://pkgdocs.julialang.org/v1/getting-started/)

!!! tip "Tip: Pkg environments"

    The package manager has excellent support
    for creating separate installation environments.
    We strongly recommend using environments to create
    isolated and reproducible setups.

    Learn more here: [https://pkgdocs.julialang.org/v1/environments/](https://pkgdocs.julialang.org/v1/environments/)

## Step 3: Run HypergraphSignals

Let's create a random hypergraph
with ``N = 15`` vertices and hyperedges of cardinality at most ``M = 3``,
along with a random signal on its vertices!
We seed the random number generator so the results are reproducible.

```@repl quickstart
using Random; Random.seed!(0);
N, M = 15, 3
hyperedges = Set(Set(rand(1:N, rand(2:M))) for _ in 1:1000);
x = randn(N);    # semicolon suppresses the output
```

Here, each hyperedge is a set of (at most ``M``) vertices from ``1, \dots, N``,
and ``x \in \mathbb{R}^N`` assigns a value to each vertex.

Now, load the package and form the hypergraph
along with its (order-``M``) adjacency tensor
``A \in \mathbb{R}^{N \times \cdots \times N}``.

```@repl quickstart
using HypergraphSignals
H = Hypergraph{M}(hyperedges, N);
A = Array(AdjacencyTensor(H));
size(A)
```

The vertex signal ``x`` is likewise expanded into
an order-``M`` hypergraph signal tensor
``X \in \mathbb{R}^{N \times 1 \times N \times \cdots \times N}``
with entries

```math
X[i_1, 1, i_3, \dots, i_M]
=
x[i_1] \cdot x[i_3] \cdots x[i_M]
.
```

```@repl quickstart
X = Array(HGSignal(x, M));
size(X)
```

The t-HGSP framework works with symmetrized versions of these tensors,
obtained by a symmetric expansion along modes ``3, \dots, M``.
The symmetrization can be undone by anti-symmetrization.

```@repl quickstart
As = hgsym(A);
size(As)
hgantisym(As) == A
```

To shift the signal ``x`` along the hypergraph,
simply run `hgshift`.
This symmetrizes both tensors
and returns their t-product ``A_s \star X_s``.

```@repl quickstart
Y = hgshift(A, x);
using HypergraphSignals: ⋆    # to write `⋆`, type `\star` then hit tab
Y ≈ hgsym(A) ⋆ hgsym(X)
```

Here, `⋆` is a public (but unexported) alias for `tmul`.

The (spectral) frequency components of the hypergraph
come from the t-eigendecomposition of the symmetrized adjacency tensor,
``A_s = V \star \Lambda \star V^{\mathsf{H}}``,
where ``V^{\mathsf{H}}`` denotes the t-adjoint.
These can be computed with `teigen` and `tadjoint`.

```@repl quickstart
Λ, V = teigen(As);
V ⋆ Λ ⋆ tadjoint(V) ≈ As
```

Finally, to filter the signal ``x`` with the hypergraph,
run `hgfilter`.
This performs the t-HGSP shift above,
anti-symmetrizes the output,
then fits a hypergraph signal to the result,
so we end up with a signal vector on the vertices again!

```@repl quickstart
y = hgfilter(A, x)
```

!!! note "Note: Current limitations of `hgfilter`"

    Fitting the hypergraph signal is currently only supported for ``M \leq 3``.
    For ``M = 3``, the fitted signal comes from a leading eigenvector,
    so it is only determined up to sign (hence the warning).

!!! tip "Congratulations!"

    Congratulations!
    You have successfully installed HypergraphSignals
    and run some t-HGSP computations!

## Next steps

Ready to learn more?
Check out the [API](@ref) page for all the available types and functions.
