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

# ╔═╡ 3cd0db02-6322-4601-bebf-44fba4f39565
using Pkg, Revise, PlutoUI

# ╔═╡ 3dee82da-6b2c-11f1-9d3c-e508785fce07
# ╠═╡ show_logs = false
Pkg.develop(path=joinpath(@__DIR__,".."))

# ╔═╡ 33008548-eb49-4695-bc62-158875005c59
using AIRSilt

# ╔═╡ 5ac5c666-fb44-402c-8fb6-9b786169e0e4
using Oscar

# ╔═╡ be6ea5fb-4eca-417a-beed-083ccc5e47be
md"# AIRSilt.jl/Playground
This reactive notebook showcases some of the methods availabel in `AIRSilt.jl`, a `Julia` tool to compute and visualise the posets of τ-tilting pairs defined by Adachi--Iyama--Reiten \[AIR\]. 

These posets arise under various other guises in nature -- 2-term silting complexes, functorially finite torsion pairs, and algebraic intermediate hearts to name a few. To play around with the resulting structures, knowing what each of these mean is not necessary, just pick one. Some learning resources are listed at the end.

To begin, load the package.
"

# ╔═╡ ecf061c5-1e9e-4dee-8a1f-02a21a65e95f

md"Documentation for functions can be viewed by clicking on \"Live docs\" in this window."

# ╔═╡ 8227c3fa-41bf-493b-a110-6a34ac62264e
md"## Initialising the algebra
The package builds on top of the `QPA` package for the `GAP` computer algebra system. To interface with `Julia`, we employ `Oscar.jl`. 
"

# ╔═╡ d9a90539-d361-4d78-af4a-5a58030bfa11
md"The simplest way to initialise an algebra is to write the `GAP` code for the algebra as a string and then call `GAP.evalstr()` on it."

# ╔═╡ 40ad318a-f42f-459c-a64e-05621bcf84f8
md"Try editing the `GAP` code to change the algebra! All methods in `AIRSilt.jl` have been implemented for finite dimensional algebras over finite fields (or the field of rational numbers). However the methods can sometimes be extended to infinite dimensional algebras via central reduction. 

This said, some other examples have been preloaded into the notebook and can be selected via the dropdown below.
"

# ╔═╡ 4e98e03a-a996-4601-9823-e073eb7004ed
@bind preset Select([
    (quiver="2, [[1,2,\"a1\"]]", 
     rels="") => "Linear A2",
    (quiver="3, [[1,2,\"a1\"],[2,3,\"a2\"]]", 
     rels="") => "Linear A3",
    (quiver="4, [[1,2,\"a1\"],[2,3,\"a2\"],[3,4,\"a3\"]]", 
     rels="") => "Linear A4",
    (quiver="5, [[1,2,\"a1\"],[2,3,\"a2\"],[3,4,\"a3\"],[4,5,\"a5\"]]", 
     rels="") => "Linear A5",
    
    (quiver="2, [[1,2,\"x\"],[1,2,\"y\"]]", 
     rels="") => "2-Kronecker",
    (quiver="2, [[1,2,\"x\"],[1,2,\"y\"],[1,2,\"z\"]]", 
     rels="") => "3-Kronecker",
    (quiver="2, [[1,2,\"x\"],[1,2,\"y\"],[1,2,\"z\"],[1,2,\"w\"]]", 
     rels="") => "4-Kronecker",

    (quiver="3, [[1,2,\"x1\"],[1,2,\"y1\"],[1,2,\"z1\"],
                          [1,2,\"x2\"],[1,2,\"y2\"],[1,2,\"z2\"]]", 
     rels="x1*y2-y1*x2, x1*z2-z1*x2, y1*z2-z1*y2") => "Beilinson(P2)",

])


# ╔═╡ 4737cdd7-cf2a-4f71-8817-45131435c326
@bind quivercode HTML("""<textarea style="font-family: 'JuliaMono'; font-size: 14px; width: 100%; height: 150px; padding: 8px; border-radius: 4px; border: 1px solid #ccc;">LoadPackage("qpa");
    quiver      := Quiver(""" * preset.quiver * """);
basefield   := GF(97); # or Rationals, or GF(4)
kQ          := PathAlgebra(basefield, quiver);
AssignGeneratorVariables(kQ);
relations   := [""" * preset.rels * """];
kQ/relations</textarea>""")

# ╔═╡ b293c3b3-a308-475f-bb7b-4692b1bf504e
algebra = GAP.evalstr("$(quivercode)")


# ╔═╡ 4a966c1e-2c09-4631-95b7-92f8c1368662
md"A $(GAP.Globals.Dimension(algebra))-dimensional algebra based on user input has been initialised above."


# ╔═╡ b17115b9-6f13-4190-8c69-5bd6b56ffb6a
md"## τ-tilting pairs
The poset of τ-tilting pairs is stored as a graph object. It can be computed by iterative left and right mutations from a starting seed (the free module of rank 1). The method below populates this poset, with an optional `timeout` that can be specified -- this is the number of seconds after which a partially computed result is returned.
"

# ╔═╡ ee829e86-f563-466e-9514-2df295b13b15
poset = tauposet(algebra, timeout=15)

# ╔═╡ b7c9bdbd-8159-4350-a360-5cc026b1540f
begin
    using MetaGraphsNext
    M_gmat = collect(labels(poset))[4]
    M = poset[M_gmat]
end

# ╔═╡ bd390c00-a283-48f5-b1ec-21067f0aa25b
begin 
    using PlutoPlotly # To display the plot object in this notebook
    PlutoPlot(gfanplot(poset, normalize=true, height=500, width=500))
end

# ╔═╡ 0e24a5c6-e4df-4722-8ea5-70035712c14b
md"The above function initialised the Hasse quiver of the τ-poset, or atleast a partial computation thereof with $(MetaGraphsNext.nv(poset)) vertices and $(MetaGraphsNext.ne(poset)) edges. 

The poset can be visualised as a directed graph, with an optional argument to display the vertices as 2-term silting complexes (default), or as the matrix of g-vectors or as a simple `\bullet`."

# ╔═╡ ca374d5f-6efa-4324-ac9a-3438ce1af9cf
tikzplot(poset, vertexlabel=:complex, vxfont="\\large", options="scale=2")

# ╔═╡ 8bfef074-4b9b-4df2-bfe2-986d2432f208
md"To access properties of the poset we use the `MetaGraphsNext.jl` package. Each vertex is labelled (and can be accessed) by its matrix of g-vectors, but it carries the data of the complete underlying module that is exposed to all `GAP` functions."

# ╔═╡ 1a0fd458-18b4-4d02-9612-949637f24b5c
md"Above we extracted a τ-tilting pair which is stored as a list of its summands. We exhibit some methods below." 

# ╔═╡ 974f24bb-9567-4d15-bcaf-71158b826e10
begin
    @info("The 2nd summand of M has K0-class $(M[2].cvec) and g-vector $(M[2].gvec).")
    @info("The result of checking tau-rigidity of this summand is: $(GAP.Globals.IsTauRigidModule(M[2].obj)).")
    @info("Is M the Bongartz completion of its 2nd summand? Answer: $(isBongartzCompletionAt(M,[1]))")

    # We can also manually compute mutation 
    _, N_gmat = mutate_in_tauposet!(poset,gmatrix(M),2)
    N = poset[N_gmat]
end

# ╔═╡ 02e1052f-f0c1-4e0d-8f1f-422841dc5338
begin
    using Graphs
    ngbd_codes = Graphs.neighborhood(poset.graph, code_for(poset,N_gmat), 2)
    map(n->poset[label_for(poset,n)], ngbd_codes)
end

# ╔═╡ a33efe3c-fa22-4e1c-ba99-1647e0a4ce4c
md"The τ-poset is also exposed to all `Graphs.jl` functions. For example below we find all elements in the poset that can be reached in at most two left-mutations from N (initialised above)."

# ╔═╡ 515402f7-230f-49a6-9f4f-6d62825bd6fc
md"## Plotting the g-fan
Associating each τ-tilting pair to the cone generated by its g-vectors collects the entire τ-poset into a simplicial fan. This is also the full-dimensional part of the heart fan. We implement a method for plotting two- and three-dimensional g-fans, with the option to normalize cone-slices to the unit sphere. Try hovering your pointer over cones and marked points!"

# ╔═╡ be147f29-0ee8-462c-babf-53fe5f6c5b39
md"## References
1. Adachi, Iyama, & Reiten. _τ-tilting theory_ (2013). [[arXiv]](http://arxiv.org/abs/1210.1036). 
3. Adachi, Iyama, & Reiten. _On τ-tilting theory_ (2024). [[doi]](https://doi.org/10.48550/arXiv.2410.15842)
4. Broomhead, Pauksztello, Ploog & Woolf. _The heart fan of an abelian category_ (2024). [[doi]](https://doi.org/10.48550/arXiv.2310.02844)
5. Cao, Gyoda & Yurikusa. _Bongartz completion via c-vectors_ (2022). [[doi]](https://doi.org/10.48550/arXiv.2106.11668)
6. Eisele, Janssens, & Raedschelders. _A reduction theorem for τ-rigid modules_ (2018). [[arXiv]](http://arxiv.org/abs/1603.04293)
7. Treffinger. _τ-tilting theory -- An introduction_ (2022). [[arxiv]](http://arxiv.org/abs/2106.00426)
"

# ╔═╡ 41eeb8dd-3813-4702-a082-9ee46c964172
md"The package `AIRSilt.jl` is written and maintained by [Parth Shimpi](https://pas201.user.srcf.net). It grew out of `GAP` code written by Theo Raedschelders and Michael Wemyss."

# ╔═╡ Cell order:
# ╟─be6ea5fb-4eca-417a-beed-083ccc5e47be
# ╠═33008548-eb49-4695-bc62-158875005c59
# ╟─ecf061c5-1e9e-4dee-8a1f-02a21a65e95f
# ╟─8227c3fa-41bf-493b-a110-6a34ac62264e
# ╠═5ac5c666-fb44-402c-8fb6-9b786169e0e4
# ╟─d9a90539-d361-4d78-af4a-5a58030bfa11
# ╟─4737cdd7-cf2a-4f71-8817-45131435c326
# ╟─b293c3b3-a308-475f-bb7b-4692b1bf504e
# ╟─4a966c1e-2c09-4631-95b7-92f8c1368662
# ╟─40ad318a-f42f-459c-a64e-05621bcf84f8
# ╟─4e98e03a-a996-4601-9823-e073eb7004ed
# ╟─b17115b9-6f13-4190-8c69-5bd6b56ffb6a
# ╠═ee829e86-f563-466e-9514-2df295b13b15
# ╟─0e24a5c6-e4df-4722-8ea5-70035712c14b
# ╠═ca374d5f-6efa-4324-ac9a-3438ce1af9cf
# ╟─8bfef074-4b9b-4df2-bfe2-986d2432f208
# ╠═b7c9bdbd-8159-4350-a360-5cc026b1540f
# ╟─1a0fd458-18b4-4d02-9612-949637f24b5c
# ╠═974f24bb-9567-4d15-bcaf-71158b826e10
# ╟─a33efe3c-fa22-4e1c-ba99-1647e0a4ce4c
# ╠═02e1052f-f0c1-4e0d-8f1f-422841dc5338
# ╟─515402f7-230f-49a6-9f4f-6d62825bd6fc
# ╠═bd390c00-a283-48f5-b1ec-21067f0aa25b
# ╟─be147f29-0ee8-462c-babf-53fe5f6c5b39
# ╟─41eeb8dd-3813-4702-a082-9ee46c964172
# ╟─3cd0db02-6322-4601-bebf-44fba4f39565
# ╟─3dee82da-6b2c-11f1-9d3c-e508785fce07
