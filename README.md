[//]: Logo
<img
    src="./docs/src/assets/logo.svg"
    width=1000px
    >

[//]: Badges
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://imperialcollegelondon.github.io/Interesso.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://imperialcollegelondon.github.io/Interesso.jl/dev)
[![Build Status](https://github.com/ImperialCollegeLondon/Interesso.jl/workflows/CI/badge.svg)](https://github.com/ImperialCollegeLondon/Interesso.jl/actions)
[![Coverage](https://codecov.io/gh/ImperialCollegeLondon/Interesso.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ImperialCollegeLondon/Interesso.jl)
 
[//]: Description
Solves Differential Equations using Integrated Residuals Methods

## Installation

```julia
using Pkg
Pkg.add("Interesso") #WIP
```

## Usage

WIP

## ODE Example

```julia
using Interesso

function f(u,p,t)
    return -1.0*u
end
u0 = 1.0 
tspan = (0.0, 10.0)

prob = ODEProblem(f,u0,tspan)
sol = solve(prob, interesso())
```
For documentation on `ODEProblem`, see [here](https://diffeq.sciml.ai/stable/types/ode_types/).

## DAE Example

```julia
using Interesso

function f(du,u,p,t)
    return du-1.0*u
end
du0 = -1.0
u0 = 1.0 
tspan = (0.0, 10.0)

prob = DAEProblem(f,du0,u0,tspan)
sol = solve(prob, interesso())
```
For documentation on `DAEProblem`, see [here](https://diffeq.sciml.ai/stable/types/dae_types/).