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
    using Oscar
    using Revise
    using AIRSilt
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
    md"The $(GAP.Globals.Dimension(algebra))-dimensional algebra has been initialised."
end

# ╔═╡ 183986d4-a294-457c-af59-1a459876bae0
begin 
    tauTilts = enumerateTauTilts(algebra, timeout = 10)
    hasseQuiver = LeftMutationQuiver(tauTilts)
    plotQuiver(hasseQuiver)
end

# ╔═╡ 9bed4558-0d48-4347-a515-b584794a9216


# ╔═╡ ecf5011f-c049-49f5-bf1d-11d83f676cd6


# ╔═╡ fd577b65-9cce-4fea-917e-5aed90d9cc71
collect(tauTilts)[2].complement

# ╔═╡ cce160e3-aeb8-4b32-8860-6d028d58de5f
'⇾' in "⇾"

# ╔═╡ Cell order:
# ╟─4737cdd7-cf2a-4f71-8817-45131435c326
# ╟─bfd34a90-dc75-45ca-ac52-45ac8fb774e4
# ╟─b293c3b3-a308-475f-bb7b-4692b1bf504e
# ╠═183986d4-a294-457c-af59-1a459876bae0
# ╠═9bed4558-0d48-4347-a515-b584794a9216
# ╠═ecf5011f-c049-49f5-bf1d-11d83f676cd6
# ╠═fd577b65-9cce-4fea-917e-5aed90d9cc71
# ╠═cce160e3-aeb8-4b32-8860-6d028d58de5f
# ╟─3dee82da-6b2c-11f1-9d3c-e508785fce07
