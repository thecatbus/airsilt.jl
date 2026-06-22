module AIRSilt
export enumerateTauTilts, IndecTauRigid, SuppTauTilting, gvectors, mutate, LeftMutationQuiver, plotQuiver, plotFan

using Oscar

function __init__()
    GAP.Packages.load("qpa")
end

include("stringutils.jl")
using .StringUtils

include("taurigids.jl")
using .TauRigids

include("tauposet.jl")
using .TauPoset

include("tauhasse.jl")
using .TauHasse

include("gfan.jl")
using .GFan

end
