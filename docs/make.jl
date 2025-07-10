push!(LOAD_PATH,"../src/")
using LinearMPC
using Documenter

makedocs(sitename = "LinearMPC.jl",
         modules  = [LinearMPC],
         pages=["Home" => "index.md"
                "Manual" => ["Getting started" => ["Simple example" => "manual/simple.md",
                                                   "Models" => "manual/model.md",
                                                   "Objective function" => "manual/objective.md",
                                                   "Constraints" =>  "manual/constraints.md",
                                                   "Code generation" =>  "manual/codegen.md",
                                           ],
                             "Beyond the basics" => ["Prestabilization" => "manual/prestab.md",
                                                     "Move blocking " => "manual/moveblock.md",
                                                     "Reference Preview" => "manual/reference_preview.md",
                                                     "Explict MPC" => "manual/explicit.md",
                                                     "Complexity Certification" => "manual/compcert.md",
                                                     "Solver Settings" => "manual/solver.md",
                                                    ]
                            ]
                "Functions" => "functions.md"
               ])

deploydocs(;repo="github.com/darnstrom/LinearMPC.jl",)
