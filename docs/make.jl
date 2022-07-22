using Interesso
using Documenter

makedocs(;
    modules=[Interesso],
    authors="Eduardo G Vila and Lucian Nita",
    repo="https://github.com/ImperialCollegeLondon/Interesso.jl/blob/{commit}{path}#L{line}",
    sitename="Interesso.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://imperialcollegelondon.github.io/Interesso.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ImperialCollegeLondon/Interesso.jl",
)
