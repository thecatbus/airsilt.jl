### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 3dee82da-6b2c-11f1-9d3c-e508785fce07
# ╠═╡ show_logs = false
begin
    using Pkg 
    Pkg.activate("/Users/parth/Projects/airsilt-2026/examples/")
    Pkg.add("TikzGraphs")
    using Oscar
    using Revise
    using AIRSilt
    using TikzGraphs
end

# ╔═╡ 4737cdd7-cf2a-4f71-8817-45131435c326
@bind quivercode html"""<textarea style="font-family: 'JuliaMono'; font-size: 14px; width: 100%; height: 100px; padding: 8px; border-radius: 4px; border: 1px solid #ccc;">quiver      := Quiver(2, [[1, 2, "x"], [1, 2, "y"], [2, 1, "z"]]);
basefield   := GF(97); # or Rationals, or GF(4)</textarea>"""

# ╔═╡ bfd34a90-dc75-45ca-ac52-45ac8fb774e4
@bind relscode html"""<textarea style="font-family: 'JuliaMono'; font-size: 14px; width: 100%; height: 50px; padding: 8px; border-radius: 4px; border: 1px solid #ccc;">relations := [x*z, z*y];</textarea>"""

# ╔═╡ b293c3b3-a308-475f-bb7b-4692b1bf504e
begin 
    algebra = GAP.evalstr("""
                $(quivercode)  
                kQ := PathAlgebra(basefield, quiver);
                AssignGeneratorVariables(kQ);
                $(relscode)
                kQ/relations
                """)
end

# ╔═╡ acf7c1a3-ba38-4a7f-82c5-7029103212fe
tauTilts = enumerateTauTilts(algebra, timeout = 10)

# ╔═╡ a6bfe6ad-42b2-4f60-b9a7-9f590eda1589
hasseQuiver = LeftMutationQuiver(tauTilts)

# ╔═╡ 183986d4-a294-457c-af59-1a459876bae0
plotQuiver(hasseQuiver; skeleton = true)

# ╔═╡ bd536da4-ee46-49c3-af61-48ae22a39e1d
# ╠═╡ disabled = true
#=╠═╡
begin
    hasseQuiver = LeftMutationQuiver(tauTilts)
    plotQuiver(hasseQuiver, skeleton = true)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─4737cdd7-cf2a-4f71-8817-45131435c326
# ╟─bfd34a90-dc75-45ca-ac52-45ac8fb774e4
# ╟─b293c3b3-a308-475f-bb7b-4692b1bf504e
# ╟─acf7c1a3-ba38-4a7f-82c5-7029103212fe
# ╟─a6bfe6ad-42b2-4f60-b9a7-9f590eda1589
# ╠═183986d4-a294-457c-af59-1a459876bae0
# ╟─bd536da4-ee46-49c3-af61-48ae22a39e1d
# ╟─3dee82da-6b2c-11f1-9d3c-e508785fce07
