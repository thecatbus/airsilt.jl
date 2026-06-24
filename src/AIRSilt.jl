# src/AIRSilt.jl

module AIRSilt
using FromFile, Oscar

function __init__()
    GAP.Packages.load("qpa")
end

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
