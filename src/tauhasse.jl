module TauHasse
export LeftMutationQuiver

using Oscar, Graphs, ..TauRigids, ..TauPoset

struct LeftMutationQuiver
    quiver::Graphs.SimpleGraphs.SimpleDiGraph{Int64}
    vertexLabels::Dict{Int16, SuppTauTilting}
    edgeLabels::Dict{Tuple{Int16,Int16}, Union{IndecTauRigid,Int16}}

    function LeftMutationQuiver(tauposet::Set{SuppTauTilting})
        quiver = SimpleDiGraph()
        vertexLabels = Dict()
        edgeLabels = Dict()

        i = 1
        for mx in tauposet
            vertexLabels[i] = mx
            Graphs.add_vertex!(quiver)
            i += 1
        end

        for mx in tauposet
            mxid = findfirst(==(mx), vertexLabels)

            for (x, myuid) in mx.mutations
                if !hasproperty(x, :gvec) 
                    continue 
                end

                gmx = hcat(x.gvec, setdiff(gvectors(mx), (x.gvec,))...)
                normal_dir = inv(gmx)[1, :]
                normal_dir = sum(normal_dir)>0 ? normal_dir : -1 * normal_dir

                if sum(normal_dir .* x.gvec) > 0
                    my = first(filter(mz -> mz.uid == myuid, tauposet))
                    myid = findfirst(==(my), vertexLabels)

                    Graphs.add_edge!(quiver, mxid, myid)
                    edgeLabels[(mxid, myid)] = x
                end
            end
        end

        return new(quiver, vertexLabels, edgeLabels)
    end

    function LeftMutationQuiver(alg :: GAP.GapObj; timeout=10)
        tauposet = enumerateTauTilts(alg; timeout=timeout)
        return LeftMutationQuiver(tauposet)
    end
end
end
