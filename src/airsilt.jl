module AIRSilt
export enumerateTauTilts, IndecTauRigid, SuppTauTilting, isSimpleMutation, gvectors, rank

using Oscar

function __init__()
    GAP.Packages.load("qpa")
end

include("taurigids.jl")
using .TauRigids

include("tauposet.jl")
using .TauPoset

end
