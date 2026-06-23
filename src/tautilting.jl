# src/taultilting.jl

using FromFile, Oscar
@from "rigidmodules.jl" import TauRigid, ShiftedProjective

1 < 0 && include("rigidmodules.jl")

""" 
    TauTilting

A type alias for an (ordered) list of `TauRigid` and `ShiftedProjective` objects. If `(M,P)` is a τ-tilting (i.e. maximal τ-rigid) pair, then the list contains indecomposable summands of `M` and `P[1]`.The order is induced by the vertices of the algebra, that is, it extends the (unique and well-defined) order under which indecomposable projectives appear at indices corresponding to their vertices.
"""
TauTilting = Vector{Union{TauRigid,ShiftedProjective}}


"""
    gmatrix(mx::TauTilting)

Returns the (square) matrix of g-vectors of a τ-tilting pair. The columns correspond to summands and respect the order induced by vertices of the algebra.
"""
gmatrix(mx::TauTilting) = hcat(getproperty.(mx,:gvec)...)


"""
    isBongartzCompletionAt(mx::TauTilting, ns::Vector{Int64})

Checks if a τ-tilting pair is the Bongartz completion of its summands *outside* the list of indices `ns`. That is, if `(M,P)=(M',P')⊕(M", P")` where the summands of `(M', P')` are indexed over `ns`, then the function returns `true` if `(M,P)` is the Bongartz completion of `(M", P")`.
"""
function isBongartzCompletionAt(mx::TauTilting, ns::Vector{Int64})
    if !all(ns .>= 1) || !all(ns .<= length(mx))
        throw(DomainError("Invalid list of indices!"))
    end

    cmatrix = inv(gmatrix(mx))'
    @views return all(cmatrix[:, ns] .>= 0)
end


"""
    isCoBongartzCompletionAt(mx::TauTilting, ns::Vector{Int64})

Checks if a τ-tilting pair is the coBongartz completion of its summands *outside* the list of indices `ns`. That is, if `(M,P)=(M',P')⊕(M", P")` where the summands of `(M', P')` are indexed over `ns`, then the function returns `true` if `(M,P)` is the coBongartz completion of `(M", P")`.
"""
function isCoBongartzCompletionAt(mx::TauTilting, ns::Vector{Int64})
    if !all(ns .>= 1) || !all(ns .<= length(mx))
        throw(DomainError("Invalid list of indices!"))
    end

    cmatrix = inv(gmatrix(mx))'
    @views return all(cmatrix[:, ns] .<= 0)
end
