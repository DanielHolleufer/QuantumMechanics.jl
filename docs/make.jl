using QuantumMechanics
using Documenter

DocMeta.setdocmeta!(QuantumMechanics, :DocTestSetup, :(using QuantumMechanics); recursive=true)

makedocs(;
    modules=[QuantumMechanics],
    authors="Daniel Holleufer <daniel.holleufer@gmail.com> and contributors",
    repo="https://github.com/DanielHolleufer/QuantumMechanics.jl/blob/{commit}{path}#{line}",
    sitename="QuantumMechanics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DanielHolleufer.github.io/QuantumMechanics.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DanielHolleufer/QuantumMechanics.jl",
    devbranch="master",
)
