module AIRSilt
export enumerateTauTilts, IndecTauRigid, SuppTauTilting, gvectors, rank, LeftMutationQuiver, plotQuiver

using Oscar

function __init__()
    GAP.Packages.load("qpa")
end

include("taurigids.jl")
using .TauRigids

include("tauposet.jl")
using .TauPoset

include("tauhasse.jl")
using .TauHasse

end
