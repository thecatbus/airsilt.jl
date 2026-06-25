# AIRSilt.jl

Tools for Adachi--Iyama--Reiten's tau-tilting theory, built upon the [OSCAR Computer Algebra System](https://www.oscar-system.org). Allows for computation and visualisation of mutation posets. 

## Installation

Because `AIRSilt.jl` is not currently registered on the official General Registry, it must be installed directly from this GitHub repository. 

Open your Julia REPL, press `]` to enter package mode, and run:

```julia
pkg> add https://github.com/thecatbus/airsilt.jl
```

Alternatively, you can add it programmatically within standard Julia scripts:

```julia
using Pkg
Pkg.add(url="https://github.com/thecatbus/airsilt.jl")
```
