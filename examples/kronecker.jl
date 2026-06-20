### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 3dee82da-6b2c-11f1-9d3c-e508785fce07
begin
    using Pkg 
    Pkg.activate("/Users/parth/Projects/airsilt-2026/examples/")
    Pkg.add("TikzGraphs")
    using Oscar
    using Revise
    using AIRSilt
    using TikzGraphs
end

# ╔═╡ 1b0a18d8-380b-4674-9422-e1aee663e34c
recon = GAP.@gap """
    Q := Quiver(2, [[1,2,"x"],[1,2,"y"],[2,1,"xx"],[2,1,"yy"],[2,1,"xy"]]);
    kQ := PathAlgebra(GF(97), Q);
    AssignGeneratorVariables(kQ);
    relations := [x*yy-y*xy, yy*x-xy*y, y*xx-x*xy, xx*y-xy*x, x*xx, xx*x, y*yy, yy*y, y*xx, xx*y, x*yy, yy*x];
    kQ/relations
    """

# ╔═╡ b293c3b3-a308-475f-bb7b-4692b1bf504e
a1 = GAP.@gap """
    Q := Quiver(2, [[1,2,"x"]]);
    PathAlgebra(Rationals, Q)
    """

# ╔═╡ 9c41449a-7c0b-4e85-a017-ca6b105b7b2d
k2 = GAP.@gap """
    Q := Quiver(2, [[1,2,"x"],[1,2,"y"]]);
    PathAlgebra(Rationals, Q)
    """

# ╔═╡ cbf8397c-dc14-47dd-b803-6577154cbd11
a2 = GAP.@gap """
    Q := Quiver(3, [[1,2,"x"],[2,3,"y"]]);
    PathAlgebra(Rationals, Q)
    """

# ╔═╡ 6f48198a-ca14-4422-b0e4-2b8231241bd6
F = enumerateTauTilts(a2, timeout = 10)

# ╔═╡ 02a25052-3e94-427f-9ae2-6de2282848ff
quiv = LeftMutationQuiver(F)

# ╔═╡ cc8a9fa7-26ed-46b4-9734-1fccd9432b58
TikzGraphs.plot(quiv.quiver)

# ╔═╡ Cell order:
# ╠═3dee82da-6b2c-11f1-9d3c-e508785fce07
# ╟─1b0a18d8-380b-4674-9422-e1aee663e34c
# ╠═b293c3b3-a308-475f-bb7b-4692b1bf504e
# ╟─9c41449a-7c0b-4e85-a017-ca6b105b7b2d
# ╠═cbf8397c-dc14-47dd-b803-6577154cbd11
# ╠═6f48198a-ca14-4422-b0e4-2b8231241bd6
# ╠═02a25052-3e94-427f-9ae2-6de2282848ff
# ╠═cc8a9fa7-26ed-46b4-9734-1fccd9432b58
