module TauRigids
export IndecTauRigid, SuppTauTilting, rank, gvectors

using Oscar, UUIDs

const SUPERSCRIPT = Dict('0'=>'⁰','1'=>'¹','2'=>'²','3'=>'³','4'=>'⁴','5'=>'⁵','6'=>'⁶','7'=>'⁷','8'=>'⁸','9'=>'⁹')
const SUBSCRIPT = Dict('0'=>'₀','1'=>'₁','2'=>'₂','3'=>'₃','4'=>'₄','5'=>'₅','6'=>'₆','7'=>'₇','8'=>'₈','9'=>'₉')
PRINTSUB(n) = String([SUBSCRIPT[digit] for digit in string(n)])
PRINTSUP(n) = String([SUPERSCRIPT[digit] for digit in string(n)])

function PRESENTGVECTOR(gvec)
    degZero = []
    degOne = []
    for (i,a) in enumerate(gvec) 
        if a == 1 
            push!(degZero, "P$(PRINTSUB(i))")
        elseif a > 1 
            push!(degZero, "P$(PRINTSUB(i))$(PRINTSUP(a))")
        elseif a == -1 
            push!(degOne, "P$(PRINTSUB(i))")
        elseif a < -1
            push!(degOne, "P$(PRINTSUB(i))$(PRINTSUP(a))")
        end
    end

    if isempty(degOne) 
        return join(degZero, "⊕")
    elseif isempty(degZero) 
        return join(degOne, "[1]⊕") * "[1]"
    else
        return join(degOne, "⊕") * "⇾" * join(degZero, "⊕") 
    end
end

struct IndecTauRigid
    obj::GAP.GapObj
    dim::Vector{Int16}
    gvec::Vector{Int16}

    function IndecTauRigid(M::GAP.GapObj)
        dimvec = GAP.gap_to_julia ∘ GAP.Globals.DimensionVector
        topOf = GAP.Globals.TopOfModule
        syzygy = GAP.evalstr("1stSyzygy")

        dim = dimvec(M)
        gvec = dimvec(topOf(M)) - dimvec(topOf(syzygy(M)))

        return new(M, dim, gvec)
    end
end

Base.isequal(M1::IndecTauRigid, M2::IndecTauRigid) = M1.gvec == M2.gvec
Base.:(==)(M1::IndecTauRigid, M2::IndecTauRigid) = M1.gvec == M2.gvec
Base.hash(M::IndecTauRigid, h::UInt) = hash(M.gvec, hash(:IndecTauRigid, h))
function Base.show(stream::IO, M::IndecTauRigid)
    result = PRESENTGVECTOR(M.gvec) 
    if '⇾' in result 
        result = "Cok("*result*")"
    end
    print(stream, result)
end

mutable struct SuppTauTilting
    uid::UUID
    summands::Set{IndecTauRigid}
    complement::Set{Int16}
    mutations::Dict{Union{IndecTauRigid,Int16},UUID}

    function SuppTauTilting(smds; algebra_rank=-1)
        uid = uuid4()
        muts = Dict()
        cmpl = if isempty(smds)
            if algebra_rank < 0
                throw(ArgumentError("To initialise the support τ-tilting module 0, the number of vertices in the algebra must be provided as an argument."))
            end
            1:algebra_rank
        else
            findall(i -> all(d[i] == 0 for d in getproperty.(smds, :dim)),
                1:length(first(smds).gvec))
        end
        return new(uid, Set{IndecTauRigid}(smds), Set{Int16}(cmpl), muts)
    end
end

Base.isequal(M1::SuppTauTilting, M2::SuppTauTilting) = M1.uid == M2.uid
Base.:(==)(M1::SuppTauTilting, M2::SuppTauTilting) = M1.uid == M2.uid
Base.hash(M::SuppTauTilting, h::UInt) = hash(M.uid, hash(:SuppTauTilting, h))
function Base.show(stream::IO, M::SuppTauTilting)
    result = "(" 
    if isempty(M.summands) 
        result = result * "0" 
    else 
        result = result * join(map(smd -> sprint(show, smd), collect(M.summands)), " ⊕ ")
    end 
    result = result * " ; "
    if isempty(M.complement) 
        result = result * "0"
    else 
        result = result * join(map(i -> "P$(PRINTSUB(i))", collect(M.complement)), " ⊕ ")
    end 
    result = result * ")"
    print(stream, result)
end

rank(M::SuppTauTilting) = length(M.summands) + length(M.complement)
gvectors(M::SuppTauTilting) =
    Set([m.gvec for m in M.summands]) ∪
    Set([[j == i ? -1 : 0 for j in 1:rank(M)] for i in M.complement])

end
