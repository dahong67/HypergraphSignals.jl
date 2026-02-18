using HypergraphSignals
using Documenter

# Setup
DocMeta.setdocmeta!(
    HypergraphSignals,
    :DocTestSetup,
    :(using HypergraphSignals);
    recursive = true,
)
ENV["LINES"] = 9

# Make docs
makedocs(;
    modules = [HypergraphSignals],
    authors = "David Hong <hong@udel.edu>, Isabel Cano <isabelc@udel.edu>, and contributors",
    sitename = "HypergraphSignals.jl",
    pages = ["Home" => "index.md"],
    format = Documenter.HTML(;
        canonical = "https://dahong67.github.io/HypergraphSignals.jl",
        edit_link = "main",
        assets = String[],
    ),
)

# Deploy docs
deploydocs(; repo = "github.com/dahong67/HypergraphSignals.jl", devbranch = "main")
