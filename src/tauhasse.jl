module TauHasse
export LeftMutationQuiver, plotQuiver

using Oscar, Graphs, TikzGraphs, ..TauRigids, ..TauPoset

struct LeftMutationQuiver
    quiver::Graphs.SimpleGraphs.SimpleDiGraph{Int64}
    vertexLabels::Dict{Int64,SuppTauTilting}
    edgeLabels::Dict{Tuple{Int64,Int64},Union{IndecTauRigid,Int16}}

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
                normal_dir = sum(normal_dir) > 0 ? normal_dir : -1 * normal_dir

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

    function LeftMutationQuiver(alg::GAP.GapObj; timeout=10)
        tauposet = enumerateTauTilts(alg; timeout=timeout)
        return LeftMutationQuiver(tauposet)
    end
end

# Pretty print the graph

function texAsComplex(M::SuppTauTilting)
    summands = []

    for gvec in gvectors(M)
        degOne = []
        degZero = []
        for (i, a) in enumerate(gvec)
            if a == 1
                push!(degZero, "P_$(i)")
            elseif a > 1
                push!(degZero, "P_$(i)^{$(a)}")
            elseif a == -1
                push!(degOne, "P_$(i)")
            elseif a < -1
                push!(degOne, "P_$(i)^{$(-a)}")
            end
        end

        if isempty(degZero)
            push!(summands, join(degOne, "[1]\\oplus ") * "[1]")
        elseif isempty(degOne)
            push!(summands, join(degZero, "\\oplus "))
        else
            push!(summands, "(" * join(degOne, "\\oplus ") * "\\rightarrow " * join(degZero, "\\oplus ") * ")")
        end
    end

    return "\\scriptsize\$" * join(summands, "\\oplus ") * "\$"
end

function plotQuiver(hasse::LeftMutationQuiver; skeleton=false)
    vxLabels = hasse.vertexLabels
    vertices = [(skeleton ? "\$\\bullet\$" : texAsComplex(vxLabels[i]))
                for i in 1:length(vxLabels)]
    return TikzGraphs.plot(hasse.quiver, vertices)
end

end
