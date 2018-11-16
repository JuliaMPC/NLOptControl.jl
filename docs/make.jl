#Pkg.add("Ipopt");Pkg.build("Ipopt")
Pkg.add("PyPlot");Pkg.build("PyPlot")
#Pkg.add("PGFPlots");Pkg.build("PGFPlots")
#Pkg.add("GR");Pkg.build("GR")

using Documenter,NLOptControl
makedocs(modules=[NLOptControl],
        doctest=false, clean=true,
        format =:html,
        authors="Huckleberry Febbo",
        assets = ["assets/style.css"],
        sitename="NLOptControl.jl",
        pages = Any[
        "Home" => "index.md",
        "Tutorials"=>Any[
              "tutorials/Brachistochrone/main.md",
              "tutorials/MoonLander/main.md",
              "tutorials/BrysonDenham/main.md",
              "tutorials/Beam/main.md",
              "tutorials/HyperSensitive/main.md",
              "tutorials/RobotArm/main.md",
              "tutorials/Rocket/main.md",
              "tutorials/Unicycle/main.md",
               ],
       "MPC"=>Any[
           "mpc/index.md"
           ]
               ]
               )

deploydocs(
    deps=nothing,
    repo="github.com/JuliaMPC/NLOptControl.jl.git",
    target="build",
    osname="linux",
    julia="0.6",
    make=nothing)

#=

     "Exported Functions"=>Any[
     "functions/NLOptControl.md",
     "functions/PrettyPlots.md"
     ]
     ]=#
