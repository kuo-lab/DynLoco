# DynLoco Documentation

`DynLoco` is a Julia package for simulating and optimizing simple dynamic locomotion models. The package includes a set of `struct`s that describe models, functions to manipulate models, and optimization routines to perform dynamic optimization on them. The optimization is based on `JuMP.jl`.

There is presently one main locomotion model, `WalkRW2l`.

# Getting started
The following demonstrates the generation of a simple walking model,
continous-time simulation of a single step, and then simulation of multiple
discrete steps. 
```julia
using DynLoco, Plots
w = WalkRW2l() # a default walking model
(t, theta, angvel) = simulatestep(w)
plot(t, theta, xlabel="time", ylabel="theta")
nsteps = multistep(w, Ps = fill(w.P, (10,1)))
plotvees(nsteps, xlabel="step #", ylabel="speed")
```



```@docs
WalkRW2l
```

