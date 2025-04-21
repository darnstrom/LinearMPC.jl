push!(LOAD_PATH,"../src/")
using LinearMPC
using Documenter

makedocs(
         sitename = "LinearMPC.jl",
         modules  = [LinearMPC],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(;
           repo="github.com/darnstrom/LinearMPC.jl",
          )
