# src/mutationposet.jl

using FromFile, Graphs, MetaGraphsNext, Oscar, TikzGraphs
using MetaGraphsNext: add_vertex!, add_edge!
@from "utils.jl" import presentgmat_tex, matrix_to_tex
@from "rigidmodules.jl" import TauRigid, ShiftedProjective
@from "tautilting.jl" import TauTilting, gmatrix, isBongartzCompletionAt

1 < 0 && include("utils.jl") && include("rigidmodules.jl") && include("tautilting.jl")


"""
    TauPoset

A mutable and edge-labelled directed graph representing the (partially computed) Hasse quiver of maximal τ-rigid pairs for a finite dimensional algebra. 

Implemented as a `MetaGraph` whose
- vertices carry data of type `TauTilting`
- vertices are labelled by, and can be referred to using the g-matrix of the corresponding τ-tilting pair
- edges correspond to left mutation in indecomposable summands
- metadata (accessible via `graph[]`) is a dictionary, with one key `:algebra` that has value the underlying algebra as a `GAP`-object.

The graph object is then exposed to methods in the `MetaGraphs.jl` package and the `Graphs.jl` package, for example, `collect(labels(...))` lists all vertex labels. See `MetaGraphs.jl` and `Graphs.jl` documentation for more.
"""
TauPoset = MetaGraph{
    Int64,                                      # Internal vertex representation
    DiGraph{Int64},                             # Type of Graph
    Matrix{Int64},                              # Type of Vertex label
    TauTilting,                                 # Type of Vertex data
    Nothing,                                    # Type of Edge data
    Dict{Symbol,Any},                           # Type of Graph data
}


"""
    add_to_tauposet!(poset::TauPoset, mx::TauTilting)

Adds the τ-tilting pair `mx` to a (possibly partially computed) `poset` of τ-tilting pairs. Automatically finds all edges to and from already existing vertices to `mx` and updates the `poset` accordingly. Returns a pair whose first element is `true` if the `poset` was updated, and `false` if the module was already in the poset. The second element is the g-matrix of `mx`, which can be used as a key to refer to the corresponding vertex.

Use carefully, since this the addition is not validated and can pollute the poset.
"""
function add_to_tauposet!(poset::TauPoset, mx::TauTilting)
    gmx = gmatrix(mx)
    haskey(poset, gmx) && return (false, gmx)

    add_vertex!(poset, gmx, mx)
    for gmz in labels(poset)
        differingcols = [i for (i, col) in enumerate(eachcol(gmx))
                         if !(col in eachcol(gmz))]
        if length(differingcols) == 1
            isBongartzCompletionAt(mx, differingcols) ?
            add_edge!(poset, gmx, gmz) :
            add_edge!(poset, gmz, gmx)
        else
            continue
        end
    end

    return (true, gmx)
end


"""
    mutate(mx::TauTilting, n::Int64; algebra::Union{GAP.GapObj, Nothing}=nothing)

Computes the mutation of a support τ-tilting module in an indecomposable summand, by explicitly computing minimal approximations. Right mutation has not been implemented yet. 

The algebra can be optionally provided as a named argument.

In most situations, it is better to set up a `TauPoset` object carrying all τ-tilting pair and call `computemutation!` instead. This may find the mutation at a summand at which only right mutation is possible.

# Arguments 
- `mx::TauTilting` : A vector of `TauRigid` and `ShiftedProjective` objects wrapping indecomposable summands of a τ-tilting pair. 
- `n::Int64` : The index of the summand at which to mutate
- `algebra::Union{GAP.GapObj,Nothing} : The algebra for which `mx` is a module. If not provided as a `GAP.GapObj`, it will be extracted by a call to `GAP.Globals.ActingAlgebra(mx)`.

# Returns
- `Tautilting` : If the left mutation `my` of `mx` in the `n`th summand is successfully computed, returns a vector of `TauRigid` and `ShiftedProjective` objects wrapping indecomposable summands of `my`.
- `nothing::Nothing` : Right mutation has not been implemented yet.
"""
function mutate(mx::TauTilting, n::Int64; algebra::Union{GAP.GapObj,Nothing}=nothing)

    x = n in eachindex(mx) ? mx[n] :
        throw(DomainError("Given module has no summand at index $(n)!"))

    algebra = !isnothing(algebra) ? algebra : GAP.Globals.ActingAlgebra(x.obj)

    if isBongartzCompletionAt(mx, [n])   # Left mutate
        if length(mx) == 1
            return TauTilting([ShiftedProjective(algebra, 1)])
        end

        smds_of_m = [mx[i] for i in eachindex(mx) if i != n && mx[i] isa TauRigid]
        m = isempty(smds_of_m) ? GAP.Globals.ZeroModule(algebra) :
            GAP.Globals.DirectSumOfQPAModules(GapObj(getproperty.(smds_of_m, :obj)))
        f = GAP.Globals.MinimalLeftApproximation(x.obj, m)
        cokerf = GAP.Globals.CoKernel(f)

        y = !iszero(GAP.Globals.Dimension(cokerf)) ?
            TauRigid(cokerf) :
        begin
            cmatrix = reduce(hcat, [mx[i] isa ShiftedProjective ? mx[i].gvec : mx[i].cvec
                                    for i in eachindex(mx) if i != n])
            support = findfirst(i -> @views(all(iszero, cmatrix[i, :])), axes(cmatrix, 1))
            ShiftedProjective(algebra, support)
        end

        return TauTilting([i == n ? y : mx[i] for i in eachindex(mx)])

    else    # Right mutate
        return nothing
    end
end


"""
    mutate_in_tauposet!(poset::TauPoset, gmx::Matrix{Int64}, n::Int64)

Finds or computes the mutation of a τ-tilting pair in an indecomposable summand. The τ-tilting pair must be a vertex in some `TauPoset`. If the (left or right) mutated pair already exists in the poset, the function simply looks up this value and returns it with no side-effect. Otherwise, the mutation is computed via minimal approximations and the value is added to the `poset`. 

Returns a pair whose first element is `true` if `poset` was modified, and `false` if not. The second element is the g-matrix of the result of mutation (right-mutation is not implemented so this will return `nothing` if the only way to reach the result is via a right-mutation).

# Arguments 
- `poset::TauPoset`: A (possibly partially computed) poset of τ-tilting pairs in which gmx exists
- `gmx::Matrix{Int64}`: The matrix of g-vectors of the τ-tilting pair to be mutated 
- `n::Int64`: The index at which mutation should be performed
"""
function mutate_in_tauposet!(poset::TauPoset, gmx::Matrix{Int64}, n::Int64)
    mx = haskey(poset, gmx) ? poset[gmx] :
         throw(ArgumentError("A τ-tilting pair with g-matrix $(gmx) hasn't been computed yet!"))

    for gmy in (inneighbor_labels(poset, gmx) ∪ outneighbor_labels(poset, gmx))
        if all(col in eachcol(gmy) for (i, col) in enumerate(eachcol(gmx)) if i != n)
            return (false, gmy)
        end
    end

    my = mutate(mx, n, algebra=poset[][:algebra])
    isnothing(my) && return (false, nothing)  # Until right mutation is implemented

    return add_to_tauposet!(poset, my)
end


"""
    tauposet(algebra::GAP.GapObj; timeout::Int64 = 15)

Generates the `TauPoset` of τ-tilting pairs of a given `algebra` by iteratively propagating mutation from the trivial τ-tilting pairs. For very large (or infinite) posets, the computation runs until a `timeout` is hit, after which the result of partial computation is returned.

# Arguments 
- `algebra::GAP.GapObj`: A finite dimensional algebra
- `timeout::Int64=15`: Number of seconds to run the computation for
"""
function tauposet(algebra::GAP.GapObj; timeout::Int64=10)
    if GAP.Globals.Dimension(algebra) == GAP.Globals.infinity
        throw(ArgumentError("Method only implemented for finite dimensional algebras! Hint: taking the quotient by an ideal generated by central elements in the radical does not change the τ-tilting poset, but may make the algebra finite dimensional."))
    end

    result = MetaGraph(DiGraph{Int64}(),
        label_type=Matrix{Int64},
        vertex_data_type=TauTilting,
        edge_data_type=Nothing,
        graph_data=Dict{Symbol,Any}(:algebra => algebra))

    seed = collect(GAP.Globals.IndecProjectiveModules(algebra))
    maxelt = TauTilting(TauRigid.(seed))
    minelt = TauTilting(ShiftedProjective.(seed))
    _, gmax = add_to_tauposet!(result, maxelt)
    _, gmin = add_to_tauposet!(result, minelt)

    precompile(mutate_in_tauposet!, (typeof(result), typeof(gmax), Int64))

    try
        timer = run(`sh -c "sleep $(timeout) && kill -2 $(getpid())"`, wait=false)
        fromtop, frombot = [gmax], [gmin]

        while !isempty(fromtop) || !isempty(frombot)
            nextfromtop = []
            nextfrombot = []

            for gmx in fromtop
                for n in axes(gmx, 2)
                    modified, gmy = mutate_in_tauposet!(result, gmx, n)
                    modified ? push!(nextfromtop, gmy) : continue
                end
            end

            for gmx in frombot
                for n in axes(gmx, 2)
                    modified, gmy = mutate_in_tauposet!(result, gmx, n)
                    modified ? push!(nextfromtop, gmy) : continue
                end
            end

            fromtop = nextfromtop
            frombot = nextfrombot
        end

        kill(timer)
        return result

    catch error
        kill(timer)
        if error isa InterruptException
            @info "Computation interrupted, possibly because timeout was hit.\nReturning partial results..."
            return result
        else
            rethrow(error)
        end
    end
end

"""
    tikzplot(poset::TauPoset; vertexlabel::Symbol=:complex, vxfont::String="", options::String="")

Returns a `TikzPictures.TikzPicture` object that is a `Tikz`-drawing of the `TauPoset`. This object can then be saved to a file, or its raw `Tikz`-code can be accessed, see `TikzPictures.jl` documentation. Optional arguments specify how vertices and edges should be rendered.

# Arguments 
- `poset::TauPoset`: The poset of τ-tilting pairs to be plotted
- `vertexlabel::Symbol`: How vertices should be represented. Accepts one of three options: 
    1. `vertexlabel=:complex` (default): displays each τ-tilting pair as its corresponding 2-term silting complex 
    2. `vertexlabel=:gmatrix`: displays each τ-tilting pair as its g-matrix 
    3. `vertexlabel=:bullet`: displays each vertex of the poset as a bullet 
- `fontoption::String`: Can be set to `"\\scriptsize"` or `"\\footnotesize"` or any such string that should be prepended on each vertex. 
- `options::String`: A string of comma-separated options that is supplied directly to `TikzGraphs.plot()`. This can be used eg to set node distance by setting `options="scale=2"`.
"""
function tikzplot(poset::TauPoset; vertexlabel::Symbol=:complex, vxfont::String="", options::String="")

    vxlabelfunc(gmat) =
        if vertexlabel == :complex
            presentgmat_tex(gmat)
        elseif vertexlabel == :gmatrix
            matrix_to_tex(gmat)
        elseif vertexlabel == :bullet
            "\$\\bullet\$"
        end

    vxLabels = String[]
    for vcode in Graphs.vertices(poset)
        gmat = label_for(poset, vcode)
        push!(vxLabels, vxfont * vxlabelfunc(gmat))
    end

    return TikzGraphs.plot(poset.graph, Layouts.Layered(), vxLabels, options=options)
end
