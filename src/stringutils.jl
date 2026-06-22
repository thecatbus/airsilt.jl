module StringUtils
export PRINTSUP, PRINTSUB, PRESENTGVECTOR

const SUPERSCRIPT = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴', '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹', '-' => '⁻')
const SUBSCRIPT = Dict('0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄', '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉', '-' => '₋')
PRINTSUB(n) = String([SUBSCRIPT[digit] for digit in string(n)])
PRINTSUP(n) = String([SUPERSCRIPT[digit] for digit in string(n)])

function PRESENTGVECTOR(gvec)
    degZero = []
    degOne = []
    for (i, a) in enumerate(gvec)
        if a == 1
            push!(degZero, "P$(PRINTSUB(i))")
        elseif a > 1
            push!(degZero, "P$(PRINTSUB(i))$(PRINTSUP(a))")
        elseif a == -1
            push!(degOne, "P$(PRINTSUB(i))")
        elseif a < -1
            push!(degOne, "P$(PRINTSUB(i))$(PRINTSUP(-a))")
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


end
