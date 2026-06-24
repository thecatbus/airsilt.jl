# src/utils.jl

const SUPERSCRIPT = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴', '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹', '-' => '⁻')
const SUBSCRIPT = Dict('0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄', '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉', '-' => '₋')

subscr(n) = String([SUBSCRIPT[digit] for digit in string(n)])
supscr(n) = String([SUPERSCRIPT[digit] for digit in string(n)])

function presentgvec(gvec::AbstractVector{Int64})
    deg0, deg1 = [], []

    for (i, a) in enumerate(gvec)
        if a == 1
            push!(deg0, "P$(subscr(i))")
        elseif a > 1
            push!(deg0, "P$(subscr(i))$(supscr(a))")
        elseif a == -1
            push!(deg1, "P$(subscr(i))")
        elseif a < -1
            push!(deg1, "P$(subscr(i))$(supscr(-a))")
        end
    end

    if isempty(deg1)
        return join(deg0, "⊕")
    elseif isempty(deg0)
        return join(deg1, "[1]⊕") * "[1]"
    else
        return "(" * join(deg1, "⊕") * "⇾" * join(deg0, "⊕") * ")"
    end
end

function presentgmat(gmat::AbstractMatrix)
    summands = presentgvec.(eachcol(gmat))
    return join(summands, "⊕")
end

function presentgmat_tex(gmat::AbstractMatrix)
    summands = []

    for gvec in eachcol(gmat)
        deg0, deg1 = [], []

        for (i, a) in enumerate(gvec)
            if a == 1
                push!(deg0, "P_{$(i)}")
            elseif a > 1
                push!(deg0, "P_{$(i)}^{\\oplus $(a)}")
            elseif a == -1
                push!(deg1, "P_{$(i)}")
            elseif a < -1
                push!(deg1, "P_{$(i)}^{\\oplus $(-a)}")
            end
        end

        if isempty(deg1)
            push!(summands, join(deg0, "\\oplus "))
        elseif isempty(deg0)
            push!(summands, join(deg1, "[1]\\oplus ") * "[1]")
        else
            push!(summands, "(" * join(deg1, "\\oplus ") *
                            "\\to " * join(deg0, "\\oplus ") * ")")
        end
    end

    return "\$" * join(summands, "\\oplus ") * "\$"
end

function matrix_to_tex(mat::AbstractMatrix)
    rows = String[]
    for i in axes(mat, 1)
        push!(rows, join(string.(mat[i, :]), " & "))
    end
    alignment = String(['r' for _ in axes(mat,2)])
    body = join(rows, " \\\\ ")
    return "\$\\left[\\begin{array}{$(alignment)} $(body) \\end{array}\\right]\$"
end
