# src/AIRSilt.jl

module AIRSilt
using FromFile, Oscar, Pluto

function __init__()
    GAP.Packages.load("qpa")

    if isinteractive()
        printstyled(raw"""
        _______________________  
         /\  | |_| |_' | |   |    Documentation: https://thecatbus.github.io/airsilt.jl
        /  \ | | \ ._| | |_  |    
        -----------------------   Type "playground()" to launch interactive notebook 
          TOOLS FOR τ-TILTING   
        =======================   Version 0.2.3 (2026-06-24)
        """)
    end

end

function playground()
    packagedir = joinpath(@__DIR__, "..")
    template = joinpath(packagedir, "examples", "playground.jl")
    examplesenv = joinpath(packagedir, "examples")

    notebook = joinpath(tempdir(), "playground.jl")
    cp(template, notebook)
    chmod(template, 0o644)

    startup_expr = "begin import Pkg; Pkg.activate(\"$(escape_string(examplesenv))\"); Pkg.instantiate() end"

    if !isfile(notebook)
        error("Notebook not found at $(notebook)")
    end

    @info "Launching interactive notebook..."

    Pluto.run(notebook=notebook,
        workspace_use_distributed=false,
        workspace_custom_startup_expr=startup_expr)
end

export playground

1 < 0 && include("rigidmodules.jl") && include("tautilting.jl") && include("mutationposet.jl") && include("gfan.jl")

@from "rigidmodules.jl" import TauRigid, ShiftedProjective
@from "tautilting.jl" import TauTilting, gmatrix, isBongartzCompletionAt, isCoBongartzCompletionAt
@from "mutationposet.jl" import TauPoset, mutate, mutate_in_tauposet!, tauposet, tikzplot
@from "gfan.jl" import gfanplot

export TauRigid, ShiftedProjective
export TauTilting, gmatrix, isBongartzCompletionAt, isCoBongartzCompletionAt
export TauPoset, mutate, mutate_in_tauposet!, tauposet, tikzplot
export gfanplot

end
