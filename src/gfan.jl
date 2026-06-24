# src/gfan.jl

using ColorSchemes, FromFile, LinearAlgebra, MetaGraphsNext, PlotlyBase
@from "utils.jl" import presentgvec, presentgmat
@from "mutationposet.jl" import TauPoset

1 < 0 && include("utils.jl") && include("rigidmodules.jl") && include("tautilting.jl") && include("mutationposet.jl")

function gconetraces(gmat::AbstractMatrix; color="gray", normalize=false)
    traces = []

    cone = collect(eachcol(gmat))
    conelabel = presentgmat(gmat)
    raylabels = [presentgvec(g) for g in cone]

    if length(cone) == 2
        if normalize
            v1, v2 = cone[1], cone[2]

            α1, α2 = atan(v1[2], v1[1]), atan(v2[2], v2[1])
            diff = α2 - α1
            if diff < 0
                diff += 2π
            end
            if diff > π
                v1, v2 = v2, v1
                α1, α2 = α2, α1
                diff = 2π - diff
            end
            angles = range(α1, α1 + diff, length=30)

            xs = [0.0; [cos(a) for a in angles]; 0.0]
            ys = [0.0; [sin(a) for a in angles]; 0.0]

            shoelace = xs[1:end-1] .* ys[2:end] .- xs[2:end] .* ys[1:end-1]
            centroid_x = sum((xs[1:end-1] .+ xs[2:end]) .* shoelace) /
                         (3 * sum(shoelace))
            centroid_y = sum((ys[1:end-1] .+ ys[2:end]) .* shoelace) /
                         (3 * sum(shoelace))
        else
            xs = [0; getindex.(cone, 1); 0]
            ys = [0; getindex.(cone, 2); 0]
            centroid_x = sum(xs) / (length(xs) - 1)
            centroid_y = sum(ys) / (length(ys) - 1)
        end

        push!(traces,
            scatter(x=xs, y=ys,
                fill="toself",
                fillcolor=color,
                mode="none",
                opacity=0.5,
                hoverinfo="skip"))
        for gvec in cone
            push!(traces,
                scatter(x=[0, gvec[1]], y=[0, gvec[2]],
                    mode="lines",
                    line=attr(color="black", width=0.3),
                    hoverinfo="skip"))
        end
        push!(traces,
            scatter(x=getindex.(cone, 1), y=getindex.(cone, 2),
                mode="markers",
                marker=attr(symbol="x", color="black"),
                hovertext=raylabels,
                hoverinfo="text",
                hoverlabel=attr(bgcolor="rgba(0,0,0,0.01)",
                    bordercolor="rgba(0,0,0,0.01)",
                    font=attr(color="black"))))
        push!(traces,
            scatter(x=[centroid_x], y=[centroid_y],
                mode="markers",
                marker=attr(opacity=0, size=30),
                hoverinfo="text",
                hovertext=conelabel,
                hoverlabel=attr(bgcolor="rgba(0,0,0,0.01)",
                    bordercolor="rgba(0,0,0,0.01)",
                    font=attr(color="black"))))

    elseif length(cone) == 3
        xs = Float64.(getindex.(cone, 1))
        ys = Float64.(getindex.(cone, 2))
        zs = Float64.(getindex.(cone, 3))

        if normalize
            norms = norm.(cone)
            xs ./= norms
            ys ./= norms
            zs ./= norms
        end

        centroid_x = sum(xs) / length(xs)
        centroid_y = sum(ys) / length(ys)
        centroid_z = sum(zs) / length(zs)

        push!(traces,
            mesh3d(x=xs .* 0.9, y=ys .* 0.9, z=zs .* 0.9,
                i=[0], j=[1], k=[2],
                color=color,
                opacity=1.0,
                hoverinfo="skip"))
        for gvec in cone
            push!(traces,
                scatter3d(x=[0, gvec[1]], y=[0, gvec[2]], z=[0, gvec[3]],
                    mode="lines",
                    line=attr(color="black", width=0.7),
                    hoverinfo="skip"))
        end
        push!(traces,
            scatter3d(x=getindex.(cone, 1), y=getindex.(cone, 2), z=getindex.(cone, 3),
                mode="markers",
                marker=attr(symbol="x", color="black", size=1.2),
                hovertext=raylabels,
                hoverinfo="text",
                hoverlabel=attr(bgcolor="rgba(0,0,0,0.01)",
                    bordercolor="rgba(0,0,0,0.01)",
                    font=attr(color="black"))))
        push!(traces,
            scatter3d(x=[centroid_x * 0.9], y=[centroid_y * 0.9], z=[centroid_z * 0.9],
                mode="markers",
                marker=attr(opacity=0, size=30),
                hoverinfo="text",
                hovertext=conelabel,
                hoverlabel=attr(bgcolor="rgba(0,0,0,0.01)",
                    bordercolor="rgba(0,0,0,0.01)",
                    font=attr(color="black"))))

    else
        throw(DomainError("Plotting only implemented for gfans of dimensions two and three!"))
    end

    return traces
end


"""
    gfanplot(tauposet::TauPoset; palette=ColorSchemes.pastel, normalize=true, width=800, height=800)

Returns a `PlotlyBase.Plot` object which plots the g-fan formed by all τ-tilting pairs in a `TauPoset`. For 3-dimensional fans, cones are represented by a slice. The default behaviour is for this slice to have endpoints on a unit sphere, this can be disabled by setting `normalize=false` (in which case the slice will be determined by g-vectors of summands).

The object is exposed to all `PlotlyBase` methods, for example `addtraces!(gfanplot(...), scatter(x=[1],y=[1]))` will add a point at `(1,1)` in the plot. See `PlotlyBase` documentation for more.
"""
function gfanplot(tauposet::TauPoset; palette=ColorSchemes.pastel, normalize=true, width=800, height=800)
    traces = GenericTrace[]
    for gmat in labels(tauposet)
        append!(traces, gconetraces(gmat, color=get(palette, rand()), normalize=normalize))
    end

    layout =
        any(t -> haskey(t, :z) ||
                (haskey(t, :type) && t[:type] in ["scatter3d", "mesh3d"]), traces) ?
        # 3D Layout
        Layout(
            width=width,
            height=height,
            plot_bgcolor="white",
            scene=attr(
                xaxis=attr(title="",
                    ticks="outside",
                    showticklabels=false,
                    tickcolor="gray",
                    ticklen=3,
                    tickwidth=1,
                    showspikes=false),
                yaxis=attr(title="",
                    ticks="outside",
                    showticklabels=false,
                    showgrid=false,
                    tickcolor="gray",
                    ticklen=3,
                    tickwidth=1,
                    showspikes=false),
                zaxis=attr(title="",
                    ticks="outside",
                    showticklabels=false,
                    showgrid=false,
                    tickcolor="gray",
                    ticklen=3,
                    tickwidth=1,
                    showspikes=false),
                xaxis_showbackground=false,
                yaxis_showbackground=false,
                zaxis_showbackground=false,
            ),
            hovermode="closest",
            showlegend=false,
        ) :
        # 2D Layout
        Layout(
            width=width,
            height=height,
            xaxis=attr(visible=false),
            yaxis=attr(visible=false, scaleanchor="x", scaleratio=1),
            plot_bgcolor="white",
            hovermode="closest",
            showlegend=false)
    push!(traces,
        scatter(x=[0], y=[0], hoverinfo="skip",
            marker=attr(symbol="cross", color="black")))
    return Plot(traces, layout)
end
