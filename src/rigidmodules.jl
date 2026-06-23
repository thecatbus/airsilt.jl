# src/rigidmodules.jl

using FromFile, Oscar
@from "utils.jl" import presentgvec

"""
    TauRigid(obj::GAP.GapObj, gvec::Vector{Int64}, cvec::Vector{Int64})
    TauRigid(obj::GAP.GapObj; validate=false)

A wrapper for a non-zero indecomposable τ-rigid module over a quotient of a path algebra. Can be initialised by simply giving the underlying module with optional validation of indecomposability and τ-rigidity.

# Fields 
- `obj::GAP.GapObj`: The underlying module as a `GAP`-object.
- `cvec::Vector{Int64}`: The dimension vector of the module as a representation of the quiver.
- `gvec::Vector{Int64}`: The g-vector of the module.
"""
struct TauRigid
    obj::GAP.GapObj
    cvec::Vector{Int64}
    gvec::Vector{Int64}
end

function TauRigid(x::GAP.GapObj, validate=false)
    if validate && !GAP.Globals.IsIndecomposableModule(x) && !GAP.Globals.IsTauRigidModule(x)
        throw(ArgumentError("Supplied module is not indecomposable τ-rigid!"))
    end

    dim = GAP.gap_to_julia ∘ GAP.Globals.DimensionVector
    top = GAP.Globals.TopOfModule
    syz = GAP.evalstr("1stSyzygy")

    cvec = dim(x)
    gvec = dim(top(x)) - dim(top(syz(x)))
    return TauRigid(x, cvec, gvec)
end


"""
    ShiftedProjective(obj::GAP.GapObj, gvec::Vector{Int64}, cvec::Vector{Int64})
    ShiftedProjective(obj::GAP.GapObj; validate::Bool=false)
    ShiftedProjective(alg::GAP.GapObj, vertex::Int64)

A wrapper for a suspended indecomposable-projective module over a quotient of a path algebra, with optional validation of indecomposability and projectivity. Can be initialised by simply giving the unsuspended module object, or by giving the algebra and the index `n` of the vertex whose idempotent generates the underlying module.

# Fields 
- `obj::GAP.GapObj`: The underlying (unsuspended) module as a `GAP`-object.
- `cvec::Vector{Int64}`: Negation of the dimension vector of the underlying quiver-representation.
- `gvec::Vector{Int64}`: The g-vector of the underlying two-term complex.
"""
struct ShiftedProjective
    obj::GAP.GapObj
    cvec::Vector{Int64}
    gvec::Vector{Int64}
end


function ShiftedProjective(x::GAP.GapObj, validate::Bool=false)
    if validate && !GAP.Globals.IsIndecomposableModule(x) && !GAP.Globals.IsProjectiveModule(x)
        throw(ArgumentError("Supplied module is not indecomposable projective!"))
    end

    dim = GAP.gap_to_julia ∘ GAP.Globals.DimensionVector
    top = GAP.Globals.TopOfModule

    cvec = -dim(x)
    gvec = -dim(top(x))
    return ShiftedProjective(x, cvec, gvec)
end

function ShiftedProjective(alg::GAP.GapObj, vertex::Int64)
    x = GAP.Globals.IndecProjectiveModules(alg)[vertex]
    dim = GAP.gap_to_julia ∘ GAP.Globals.DimensionVector
    top = GAP.Globals.TopOfModule

    cvec = -dim(x)
    gvec = -dim(top(x))
    return ShiftedProjective(x, cvec, gvec)
end

Base.show(io::IO, x::Union{TauRigid,ShiftedProjective}) =
    print(io, ('⇾' in presentgvec(x.gvec) ? "Cok" : "") * presentgvec(x.gvec))
