### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 3dee82da-6b2c-11f1-9d3c-e508785fce07
begin
    using Pkg 
    Pkg.activate("/Users/parth/Projects/airsilt-2026/examples/")
    using Oscar
    using Revise
    using AIRSilt
end

# ╔═╡ 1b0a18d8-380b-4674-9422-e1aee663e34c
a1 = GAP.@gap """
    Q := Quiver(2, [[1,2,"x"],[1,2,"y"],[2,1,"xx"],[2,1,"yy"],[2,1,"xy"]]);
    kQ := PathAlgebra(GF(97), Q);
    AssignGeneratorVariables(kQ);
    relations := [x*yy-y*xy, yy*x-xy*y, y*xx-x*xy, xx*y-xy*x, x*xx, xx*x, y*yy, yy*y, y*xx, xx*y, x*yy, yy*x];
    kQ/relations
    """

# ╔═╡ e9dfc656-8592-4107-880d-b4ce92717a9c
F = enumerateTauTilts(a1, timeout=1)

# ╔═╡ 6f48198a-ca14-4422-b0e4-2b8231241bd6
[gvectors(m) for m in F]

# ╔═╡ 02a25052-3e94-427f-9ae2-6de2282848ff
length(F)

# ╔═╡ Cell order:
# ╠═3dee82da-6b2c-11f1-9d3c-e508785fce07
# ╠═1b0a18d8-380b-4674-9422-e1aee663e34c
# ╠═e9dfc656-8592-4107-880d-b4ce92717a9c
# ╠═6f48198a-ca14-4422-b0e4-2b8231241bd6
# ╠═02a25052-3e94-427f-9ae2-6de2282848ff
