push!(LOAD_PATH, "../src/")

using Documenter
using JuAFEM

makedocs(
    format=Documenter.HTML(),
    sitename="JuAFEM.jl",
    modules=[JuAFEM],
    pages=[
        ##############################################
        ## MAKE SURE TO SYNC WITH docs/src/index.md ##
        ##############################################
        "Basics" => [
            "index.md",
            "install.md",
            "get_started.md"
           ],
        
        "faq.md"

    ],
)

deploydocs(
    repo="https://github.com/Dibakar58/JuAFEM.jl",
    devbranch="main"
)