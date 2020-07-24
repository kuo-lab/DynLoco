module DynLoco

# Locomotion models

using Parameters #, DifferentialEquations
using JuMP, Ipopt, Plots, Setfield
using StructArrays
using Dierckx # spline package

export Walk, AbstractWalkRW2, WalkRW2l, Parms, onestep
export findgait

abstract type Walk end
abstract type AbstractWalkRW2 <: Walk end

"""
    w = WalkRW2l(vm, P, α, γ, L, M, g, parms, limitcycle)

Rimless wheel in 2D, linearized. Struct containing model parameters. Parameter values can be
overridden with `WalkRW2l(w, α=0.35)` which produces a new model like `w` except with a new
value for `alpha`. Similarly any parameter may be replaced using named
keyword arguments. State: mid-stance speed `vm`. Control: Push-off impulse `P`.

`parms` is a list of numerical model parameters.
`limitcycle` is a named tuple with f, a function whose root is a fixed point.
`safety` (logical) prevents model from falling backwards, rather it always moves forward,
albeit very slowly if given insufficient forward speed.
"""
@with_kw struct WalkRW2l <: AbstractWalkRW2
    vm = 0.35 # initial mid-stance velocity (stance leg upright)
    P = 0.1   # push-off impulse (mass-normalized, units of Δv)
    α = 0.3 # angle between legs
    γ = 0.  # downward slope
    L = 1.  # leg length
    M = 1.  # body mass
    g = 1.  # gravitational acceleration
    parms = (:α, :γ, :L, :M, :g) # a list of the model parameters
    limitcycle = (parms = (:vm,), f = w -> onestep(w).vm - w.vm)
    safety = false # model falls backward if not enough momentum
end

export StepResults, MultiStepResults

"""
    steps = StepResults(vm, θnew, tf, P, C, Pwork, Cwork, speed, steplength,
        tf1, tf2, Ωminus, Ωplus)

Struct containing step outcomes, including `vm` = middle stance velocity following
step-to-step transition, `θnew` = angle of new stance leg (+ccw about vertical),
`tf` = step time (mid-stance to mid-stance), `P`/`C` = push-off/collision,
`Pwork`/`Cwork` = work, `tf1`/`tf2` = times until s2s transition, `Ωminus`/`Ωplus` =
angular velocities of model before and after s2s.

`StepResults` is a StructArray, and can be referred to by step index `steps[i]`,
or by field name `steps.tf` (an array of step times).
"""
struct StepResults
    vm    # mid-stance velocity after s2s transition
    θnew  # angle of new stance leg after s2s transition (+ccw from vertical)
    tf    # step time (mid-stance to mid-stance)
    P
    C
    Pwork
    Cwork
    speed
    steplength
    stepfrequency
    tf1
    tf2
    Ωminus
    Ωplus
    vm0    # middle-stance velocity at start of step
    δ      # slope wrt horizontal
end

function Base.show(io::IO, w::W) where W <: Walk
    print(io, W, ": ")
    for (i,field) in enumerate(fieldnames(W))
        print(io, "$field = $(getfield(w, field))")
        if i < length(fieldnames(W))
            print(io, ", ")
        else
            println(io)
        end
    end
end


function Base.show(io::IO, s::StepResults)
    print(io, "StepResults: ")
    for (i,field) in enumerate(fieldnames(StepResults))
        print(io, "$field = $(getfield(s, field))")
        if i < length(fieldnames(StepResults))
            print(io, ", ")
        else
            println(io)
        end
    end
end


"""
    msr = MultiStepResults(steps, totalcost, totaltime, vm0, δangles, boundaryvels)

Struct containing outputs of a multistep. Can also be fed into multistepplot.
The steps field is a StepResults struct, which can be referred by index, e.g.
`msr.steps[1]` and also by fields, e.g. `msr.steps.tf`.
"""
struct MultiStepResults
    steps::StructArray{StepResults}
    totalcost
    totaltime
    vm0             # initial speed
    δangles         # angles
    boundaryvels::Tuple
end

export MultiStepResults


function Base.show(io::IO, s::MultiStepResults)
    print(io, "MultiStepResults: ")
    for (i,field) in enumerate(fieldnames(MultiStepResults))
        print(io, "$field = $(getfield(s, field))")
        if i < length(fieldnames(MultiStepResults))
            print(io, ", ")
        else
            println(io)
        end
    end
end

getfields(w::Walk, fields) = map(f->getfield(w,f), fields)

"""
    onestep(walk [,parameters, safety])

Take one step with a `walk` struct, from mid-stance to mid-stance. Parameter
values can be supplied to override the existing field values. Notably,
includes vm, P, δangle.
Output is a `NamedTuple` of outcomes such as `vm` (at end), `tf` (final time),
`work`, `speed`, `steplength`, etc.
Set `safety=true` to prevent step from blowing up if stance leg has too little
momentum to reach middle stance. A step will be returned with very long
step time, and very slow mid-stance velocity.
"""
function onestep(w::WalkRW2l; vm=w.vm, P=w.P, δangle = 0.,
    α=w.α, γ=w.γ,g=w.g,L=w.L,safety=w.safety)
    mylog = safety ? logshave : log # logshave doesn't blow up on negative numbers
    # Start at mid-stance leg vertical, and find the time tf1
    # to heelstrike, and angular velocity Ωminus
    # δangle is the upward change in slope from nominal γ
    # Phase 1: From mid-stance to just before impact
    # θ is angle ccw from surface normal (which may be at slope)
    θf = -α + δangle # angle of stance leg just before impact, negative means cw
    # angle wrt vertical is -α-γ+δangle
    #Ωminus = -√(2*g*(cos(0)-cos(θf))+L*vm^2)/√L  # pre-impact speed
    Ωminus = -√(2*g/L*(cos(-γ)-cos(θf-γ))+(vm/L)^2)  # pre-impact speed
    tf1 = mylog((L*(γ-θf)+√(vm^2-2γ*θf*L^2+θf^2*L^2))/(vm+L*γ)) # time until impact, phase 1
    # Step-to-step transition: Push-off and collision
    Pwork = 1/2 * P^2
    v0 = √(Ωminus^2 + P^2) # intermediate velocity after push-off
    Ωplus = cos(2α)*Ωminus - P*sin(2α) # initial ang vel of next leg
    Cwork = 1/2*(Ωplus^2 - v0^2) # negative collision work
    C = √(-2Cwork)
    # Phase 2: From just after impact to next mid-stance
    θnew = α + δangle # new stance leg's initial angle
    #    tf2 = mylog((γ + √(2γ*θnew-θnew^2+Ωplus^2))/(γ-θnew-Ωplus))
    gto = γ - θnew - Ωplus # needs to be positive for pendulum to reach next mid-stance
    if safety && gto <= 0
        tf2 = 1e3 - gto # more negative gto is worse, increasing a tf2 time
    else # safety off OR gto positive
        # time to mid-stance, phase 2
        inroot = (2γ*θnew-θnew^2+Ωplus^2)
        if inroot >= 0
#            tf2 = mylog((γ + √(2γ*θnew-θnew^2+Ωplus^2))/(gto)) # time to mid-stance, phase 2
            tf2 = mylog((γ + √inroot)/(gto)) # time to mid-stance, phase 2
        else # inroot negative,
            tf2 = 1e4 - inroot # more negative inroot extends time
        end

    end
    twicenextenergy = (2g*L*cos(θnew-γ)+L^2*Ωplus^2-2g*L*cos(-γ)) # to find next vm
    if twicenextenergy >= 0.
        vmnew = √twicenextenergy
    elseif safety # not enough energy
        vmnew = (1e-3)*exp(twicenextenergy)
    else # no safety, not enough energy
        vmnew = √twicenextenergy # this should fail
    end

    # Step metrics
    steplength = 2L*sin(α) # rimless wheel step length
    tf = tf1 + tf2 # total time mid-stance to mid-stance
    speed = steplength / tf # average speed mid-stance to mid-stance
    return (vm=vmnew, θnew=θnew, tf=tf, P=P, C=C, Pwork=Pwork,Cwork=Cwork,
        speed=speed, steplength=steplength, stepfrequency=speed/steplength,tf1=tf1, tf2=tf2,
        Ωminus=Ωminus,Ωplus=Ωplus,
        vm0=vm,δ=δangle)
end



"""
    wnew = findgait(walk [,target=target] [,varying=varying] [,parameters])

Find a limit cycle for a `walk` struct, using optional new `parameters` specified.
Optional `target` is a desired output in symbol form accompanied by associated
`varying` parameter, with 1-to-1 `target` and `varying`. Returns a new walk struct
that satisfies limit cycle and desired target.

# Examples
        findgait(walk, target=(:speed=>0.4), varying=(:P), P=0.15)
Finds a limit cycle with walk starting at P=0.15, but varying P to yield
a desired speed of 0.4.

Note that `target` and `varying` may be singletons of tuples, with each
referring to a parameter symbol. Other parameters are set as equalities.
"""
function findgait(w::W; target=(), varying=(), walkparms...) where W <: Walk
    if target isa Tuple{Vararg{Pair}} && varying isa Tuple
        return findgait(W(w, walkparms), target, varying)
    elseif target isa Pair && varying isa Symbol
        return findgait(W(w, walkparms), (target,), (varying,))
    else
        error("findgait unable to parse target and varying")
    end
end
# Can be called with target and varying tuples, e.g.
# target=(:speed=>0.4, :steptime=0.6), varying=(:P, :Kp)

# Non-keyword form, with ordered arguments: w, target, varying
function findgait(w::W, target::Tuple{Vararg{Pair}}, varying::Tuple) where W <: Walk
    # target = (:speed=>0.4, ) or even nothing--there's a function related to this
    # varying = (:P, :α) a tuple of symbols representing model parameters
    varying = (varying..., w.limitcycle.parms...)   # we always add vm as a limit cycle variable
    target = (target..., # and a target, where 0 is dummy
        map(Pair,w.limitcycle.parms,zeros(length(w.limitcycle.parms)))...)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
    nvars = length(varying)
    @assert nvars == length(target) "number of target should equal number of varying"
    startvalues = getfields(w, varying)  # initial guess taken from Walk parameters
    @variable(model, var[i=1:nvars], start=startvalues[i]) # JuMP variables var[i]

    # Set up target
    # For nonlinear constraints, produce functions to send to JuMP as
    # fv1(var[1], var[2],...) == 0. In Julia these are an array of
    # anonymous functions that use JuMP "vars" as input and yield
    # appropriate target fields, equivalent to onestep().speed - speedtargert.
    fv=Array{Function}(undef,nvars) # array of functions with vars as input
    for i in 1:nvars # turn each :speed=>0.4 into onestep().speed - speedtarget
        fv[i] = function(vars...) # fv = function of vars
            # convert vars[...] into parameter assignments for onestep
            assignedparms = (; zip(varying, vars)...) # P = var[1], etc.
            # Ex: given target=:speed =>0.4, do onestep(w; P=var[1]).speed - 0.4
            return getfield(
                onestep(w; assignedparms...), target[i].first) - target[i].second
        end # function
        # each fv[i] is a registered JuMP function fvi with automatic differentiation.
        register(model, Symbol("fv", i), nvars, fv[i], autodiff=true)
        if i < nvars # each constraint function should be matched to its target value
            JuMP.add_NL_constraint(model, :($(Symbol("fv",i))($(var...)) == 0.))
        else # Final constraint is always limit cycle onestep.vm == vm
            # this parses to fvi[var[1],var[2]...] == var[nvars] (which is vm)
            JuMP.add_NL_constraint(model, :($(Symbol("fv",i))($(var...)) == $(var[nvars])))
        end
    end # loop over targets
    optimize!(model)
    if termination_status(model) == MOI.LOCALLY_SOLVED || termination_status(model) == MOI.OPTIMAL
        optimal_solution = value.(var[i] for i in 1:nvars) # e.g. P, vm, etc. (JuMP.value returns optima)
    else
        error("The model was not solved correctly.")
        println(termination_status(model))
    end
    # produce a new Walk struct
    assignedoptimalparms = (; zip(varying, optimal_solution)...) # P = 0.17, etc.
    return W(w; assignedoptimalparms...)
end

"""
    multistep(walk; Ps=Ps)

Take multiple steps with a `walk` struct, from mid-stance to mid-stance, with as many
steps as in push-off array `Ps`. Returns a `MultiStepResults` struct containing individual steps.

Optional named keyword version can be called with any order of arguments after
`walk`:
    multistep(walk; Ps = [0.15, 0.15], δangles = [-0.1, 0.1], vm0=walk.vm, boundaryvels=(),
        extracost = 0, walkparms...)

where slope angles `δangles` (default level ground), `vm0` initial mid-stance velocity,
`boundaryvels` boundary velocities, `extracost` added onto push-off work, and `walkparms`
can be used to apply model parameters such as `α`, `γ`, `M`, etc.
"""
function multistep(w::W; Ps=w.P*ones(5), δangles=zeros(length(Ps)), vm0 = w.vm, extracost = 0,
    boundaryvels=(), walkparms...) where W <: Walk
    return multistep(W(w; walkparms...), Ps, δangles, boundaryvels=boundaryvels, extracost = extracost)
end

function multistep(w::W, Ps::AbstractArray, δangles=zeros(length(Ps)); vm0 = w.vm,
    boundaryvels=(), extracost = 0) where W <: Walk
    steps = StructArray{StepResults}(undef, length(Ps))
    for i in 1:length(Ps)
        steps[i] = StepResults(onestep(w, P=Ps[i], δangle=δangles[i])...)
        w = W(w, vm=steps[i].vm)
    end
    totalcost = sum(steps[i].Pwork for i in 1:length(Ps)) + extracost
    # 1/2 (v[1]^2-boundaryvels[1]^2)
    totaltime = sum(getfield.(steps,:tf))
    finalvm = w.vm
    return MultiStepResults(steps, totalcost, totaltime, vm0, δangles, boundaryvels)
    # boundaryvels is just passed forward to plots
end

export multistep

totalwork(steps::StructArray) = sum(steps.Pwork)
totaltime(steps::StructArray) = sum(steps.tf)
export totalwork, totaltime

"""
    multistepplot(::MultiStepResults)

A Plots.jl recipe for plotting output of `multistep` as 3x1 subplots with mid-stance velocity,
push-off impulse (or work), and terrain profile. To plot push-off work, use keyword argument
`plotwork=true`.
"""
multistepplot

@userplot MultiStepPlot

@recipe function f(h::MultiStepPlot; plotwork=false, boundarywork=true)

    markershape --> :circle

    msr = h.args[1]
    v = [msr.vm0; msr.steps.vm] # all velocities
    P = msr.steps.P
    δ = msr.δangles
    boundaryvels = msr.boundaryvels

    #noslope = all(δ .== 0.)
    doslope = true

    n = length(msr.steps) # how many steps

    # set up the subplots: v, P, slopes
    #legend := false
    #link := :both
    grid := false
    if !doslope
        layout := @layout [vplot
                           Pplot]
    else
        layout := @layout [vplot
                           Pplot
                           δplot]
    end

    # vplot
    @series begin
        subplot := 1
        ylabel := "V midstance"
        ylims := (0, Inf)
        if boundaryvels == nothing || isempty(boundaryvels)
            0:n, v
        else
            [0; 0:n; n], [boundaryvels[1]; v; boundaryvels[2]]
        end
    end

    # Pplot or workplot, interleaved with step numbers, with
    @series begin
        subplot := 2
        ylabel := plotwork ? "Work" : "P"
        ylims := (0, Inf)
        #xlims := (0.5, n+0.5)
        plotwork ? ([0.5; 1:n; n+0.5], [1/2*(v[1]^2-boundaryvels[1]^2); 1/2 .* P.^2; NaN]) :
                   ([0.5;1:n;n+0.5], [v[1]-boundaryvels[1]; P; NaN]) # just plot the push-offs
    end

    if doslope
        @series begin
            subplot := 3
            ylabel := "Terrain height"
            xlabel := "Step number"
            0:n, cumsum([0.; δ])
        end
    end
end

using JuMP, Ipopt
export optwalk, optwalkslope

"""
    optwalk(w::Walk, numsteps=5)

Optimizes push-offs to walk `numsteps` steps. Returns a `MultiStepResults`
struct. Allows slopes to be specified as keyword `δs` array of the slope of each successive
step.

Other keyword arguments: `boundaryvels = (vm,vm)` can specify a tuple of initial and
final speeds, default nominal middle-stance speed `vm`. To start and end at rest, use `(0,0)`.
`boundarywork = true` whether cost includes the work needed to
start and end from `boundaryvels`. `totaltime` is the total time to take the steps, by default
the same time expected of the nominal `w` on level ground.

See also `optwalkslope`
"""
function optwalk(w::W, numsteps=5; boundaryvels::Union{Tuple,Nothing} = nothing,
    boundarywork = true, totaltime=numsteps*onestep(w).tf,
    δ = zeros(numsteps)) where W <: Walk # default to taking the time of regular steady walking

    optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
    @variable(optsteps, P[1:numsteps]>=0, start=w.P) # JuMP variables P
    @variable(optsteps, v[1:numsteps+1]>=0, start=w.vm) # mid-stance speeds

    if boundaryvels == nothing || isempty(boundaryvels)
        boundaryvels = (w.vm, w.vm) # default to given gait if nothing specified
    end

    if !boundarywork # no hip work at beginning or end; apply boundary velocity constraints
        @constraint(optsteps, v[1] == boundaryvels[1])
        @constraint(optsteps, v[numsteps+1] == boundaryvels[2])
    end

    # Constraints
    # produce separate functions for speeds and step times
    register(optsteps, :onestepv, 3, # velocity after a step
        (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).vm, autodiff=true) # output vm
    register(optsteps, :onestept, 3, # time after a step
        (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).tf, autodiff=true)
    @NLexpression(optsteps, summedtime, # add up time of all steps
        sum(onestept(v[i],P[i],δ[i]) for i = 1:numsteps))
    for i = 1:numsteps  # step dynamics
        @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i],δ[i]))
    end
    @NLconstraint(optsteps, summedtime == totaltime) # total time

    if boundarywork
        @objective(optsteps, Min, 1/2*(sum((P[i]^2 for i=1:numsteps))+v[1]^2-boundaryvels[1]^2)) # minimum pos work
    else
        @objective(optsteps, Min, 1/2*sum((P[i]^2 for i=1:numsteps))) # minimum pos work
    end
    optimize!(optsteps)
    if termination_status(optsteps) == MOI.LOCALLY_SOLVED || termination_status(optsteps) == MOI.OPTIMAL
        optimal_solution = (vms=value.(v), Ps=value.(P))
    else
        error("The model was not solved correctly.")
        println(termination_status(optsteps))
    end
    return multistep(W(w,vm=value(v[1])), value.(P), δ, vm0=value(v[1]), boundaryvels=boundaryvels,
        extracost = boundarywork ? 1/2*(value(v[1])^2 - boundaryvels[1]^2) : 0) #, optimal_solution
end

"""
    optwalkslope(w<:Walk, numsteps=5)

Optimizes push-offs and terrain profile to walk `numsteps` steps. Returns a `MultiStepResults`
struct. Optional keyword arguments: `boundaryvels = (vm,vm)` can specify a tuple of initial and
final speeds, default nominal middle-stance speed `vm`. To start and end at rest, use `(0,0)`.
`boundarywork = true` whether cost includes the work needed to
start and end from `boundaryvels`. `symmetric = false` whether to enforce fore-aft symmetry in
the slope profile. `totaltime` is the total time to take the steps, by default the same
time expected of the nominal `w` on level ground.

See also `optwalk`
"""
function optwalkslope(w::W, numsteps=5; boundaryvels::Union{Tuple,Nothing} = nothing,
    boundarywork = true, symmetric=false,
    totaltime=numsteps*onestep(w).tf) where W <: Walk # default to taking the time of regular steady walking
    optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
    @variable(optsteps, P[1:numsteps]>=0, start=w.P) # JuMP variables P
    @variable(optsteps, δ[1:numsteps], start=0.)         # delta slope
    @constraint(optsteps, sum(δ[i] for i =1:numsteps) == 0.)  # zero height gain
    @variable(optsteps, v[1:numsteps+1]>=0, start=w.vm) # mid-stance speeds

    if boundaryvels == nothing
        boundaryvels = (w.vm, w.vm) # default to given gait if nothing specified
    end

    if !boundarywork # no hip work at beginning or end; apply boundary velocity constraints
        @constraint(optsteps, v[1] == boundaryvels[1])
        @constraint(optsteps, v[numsteps+1] == boundaryvels[2])
    end

    # Constraints
    # produce separate functions for speeds and step times
    register(optsteps, :onestepv, 3, # velocity after a step
        (v,P,δ)->onestep(w,P=P,vm=v,δangle=δ).vm, autodiff=true) # output vm
    register(optsteps, :onestept, 3, # time after a step
        (v,P,δ)->onestep(w,P=P,vm=v,δangle=δ).tf, autodiff=true)
    @NLexpression(optsteps, summedtime, # add up time of all steps
        sum(onestept(v[i],P[i],δ[i]) for i = 1:numsteps))
    for i = 1:numsteps  # step dynamics
        @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i],δ[i]))
    end
    @NLconstraint(optsteps, summedtime == totaltime) # total time

    if symmetric
        for j in 1:numsteps÷2
            @constraint(optsteps, δ[numsteps-j+1] == -δ[j])
        end
        if isodd(numsteps)
            @constraint(optsteps, δ[numsteps÷2+1] == 0)
        end
    end

    if boundarywork
        @objective(optsteps, Min, sum(P[i]^2 for i=1:numsteps)+v[1]^2-boundaryvels[1]^2) # minimum pos work
    else
        @objective(optsteps, Min, sum(P[i]^2 for i=1:numsteps)) # minimum pos work
    end

    optimize!(optsteps)
    if termination_status(optsteps) == MOI.LOCALLY_SOLVED || termination_status(optsteps) == MOI.OPTIMAL
        optimal_solution = (vms=value.(v), Ps=value.(P), δs=value.(δ))
    else
        error("The model was not solved correctly.")
        println(termination_status(optsteps))
    end
    return multistep(W(w,vm=value(v[1])), value.(P), value.(δ), vm0=value(v[1]),
        boundaryvels = boundaryvels,
        extracost = boundarywork ? 1/2*(value(v[1])^2-boundaryvels[1]^2) : 0)
end


function logshave(x, xmin=1e-10)
    x >= xmin ? log(x) : log(xmin) + (x - xmin)
end # logshave

export plotvees, plotvees!

"""
    plotvees(results::MultiStepResults [, tchange = 3, boundaryvels = (0.,0.))

Plots a series of discrete speeds for multiple steps, along with a spline
to connect the discrete points. Returns a `Plot` struct. See also plotvees!(p, ...).
"""
plotvees(results::MultiStepResults; veeparms...) = plotvees!(plot(), results; veeparms...)

"plotvees!([p,] results::MultiStepResults) adds to an existing plot."
plotvees!(results::MultiStepResults; veeparms...) = plotvees!(Plots.CURRENT_PLOT.nullableplot, results; veeparms...)

function plotvees!(p::Union{Plots.Plot,Plots.Subplot}, msr::MultiStepResults; tchange = 3, boundaryvels = (0.,0.),
    color = :auto, usespline=true, rampuporder = 2, plotoptions...)
    v = [msr.vm0; msr.steps.vm] # vm0 is the speed at beginning of first step, vm is the mid-stance speed of first step
    n = length(msr.steps)
    times = cumsum([tchange; msr.steps.tf]) # add up step times, starting from ramp-up
    # make smooth ramp-up in speed
    t0 = range(0, tchange, length=10)
    vstart = v[1]*(t0/tchange).^rampuporder      # ramp-up in speed, using a quadratic for now
    vend = v[n+1]*(1 .- t0/tchange).^rampuporder # ramp-down in speed
    if usespline     # make a smooth spline from discrete velocities
        k = 2 # spline order
        if length(v) > 2 # enough points to make splines from v alone
            spline = Spline1D(times, v; k=k)
        else # only say 1 point, so let's pad v with the ramp-up ramp-down
            spline = Spline1D([t0[1]; times; t0[end]+times[end]],
                [vstart[1]; v; vend[end]], k=k)
        end
        tspline = range(times[1], times[end], length=20)
        plot!(p,[t0; tspline; t0 .+ times[end]], [vstart; spline.(tspline); vend], color=color;
            plotoptions...)
    end
    plot!(p, times, v, seriestype=:scatter, legend=:none, color=color, # dots for discrete velocities
        xlabel="time", ylabel="speed"; plotoptions...)
end


using JuMP, Ipopt
export optwalktime
# ctime is the cost of time
"""
    optwalktime(w::WalkrW2l, numsteps=5)

Optimizes push-offs to walk `numsteps` steps in minimum work and time. Returns a `MultiStepResults`
struct. Optional keyword arguments: `δ` array of slopes for uneven terrain. `ctime` is the relative
cost of time in units of work/time.

Other keywords: `boundaryvels = (vm,vm)`
can specify a tuple of initial and
final speeds, default nominal middle-stance speed `vm`. To start and end at rest, use `(0,0)`.
`boundarywork = true` whether cost includes the work needed to
start and end from `boundaryvels`. `symmetric = false` whether to enforce fore-aft symmetry in
the slope profile. `totaltime` is the total time to take the steps, by default the same
time expected of the nominal `w` on level ground. `startv` initial guess at speeds.

See also `optwalk` and `optwalkslope`
"""
function optwalktime(w::W, nsteps=5; boundaryvels::Union{Tuple,Nothing} = (0.,0.), safety=true,
    ctime = 0.05, tchange = 3., boundarywork = true, δs = zeros(nsteps), startv = w.vm, walkparms...) where W <: Walk
    w = W(w; walkparms...)
    optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>1))
    @variable(optsteps, P[1:nsteps]>=0, start=w.P) # JuMP variables P
    # constraints: start
    if length(startv) == 1
        @variable(optsteps, v[1:nsteps+1]>=0, start=startv)
    else # starting guess for v is an array
        @variable(optsteps, v[i=1:nsteps+1]>=0, start=startv[i])
    end
    register(optsteps, :onestepv, 3, (v,P,δ)->onestep(w,P=P,vm=v,δangle=δ,safety=safety).vm, autodiff=true) # input P, output vm
    register(optsteps, :onestept, 3, (v,P,δ)->onestep(w,P=P,vm=v,δangle=δ,safety=safety).tf, autodiff=true)
    @NLexpression(optsteps, steptime[i=1:nsteps], onestept(v[i],P[i],δs[i]))
    @NLexpression(optsteps, totaltime, sum(onestept(v[i],P[i],δs[i]) for i = 1:nsteps))
    for i = 1:nsteps # collocation points
        @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i],δs[i]))
    end
    if boundarywork
        @NLobjective(optsteps, Min, 1/2*(sum(P[i]^2 for i=1:nsteps) + v[1]^2 - boundaryvels[1]^2) +
            ctime*totaltime)
    else
        @NLobjective(optsteps, Min, 1/2*(sum(P[i]^2 for i=1:nsteps)) +
            ctime*totaltime)
    end
    optimize!(optsteps)
    result = multistep(W(w,vm=value(v[1]),safety=safety), value.(P), δs, vm0=value(v[1]),
        boundaryvels=boundaryvels, extracost = ctime*value(totaltime) +
        (boundarywork ? 1/2*(value(v[1])^2-boundaryvels[1]^2) : 0))
    return result
end

export islimitcycle
"""
    islimitcycle(w<:Walk)

Checks whether `w` is a limit cycle, i.e. taking one step
returns nearly the same gait conditions for the next step.
"""
islimitcycle(w::Walk) = onestep(w).vm ≈ w.vm


# Optimize walking to match a velocity profile
# function optwalkvel(w::Walk, vels; boundaryvels::Union{Tuple,Nothing} = nothing,
#     boundarywork = true,
#     δ = zeros(numsteps))
#
#     optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
#     @variable(optsteps, P[1:numsteps]>=0, start=w.P) # JuMP variables P
#
#     if boundaryvels == nothing || isempty(boundaryvels)
#         boundaryvels = (w.vm, w.vm) # default to given gait if nothing specified
#     end
#
#     if !boundarywork # no hip work at beginning or end; apply boundary velocity constraints
#         @constraint(optsteps, v[1] == boundaryvels[1])
#         @constraint(optsteps, v[numsteps+1] == boundaryvels[2])
#     end
#
#     # Constraints
#     # produce separate functions for speeds and step times
#     register(optsteps, :onestepv, 3, # velocity after a step
#         (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).vm, autodiff=true) # output vm
#     register(optsteps, :onestept, 3, # time after a step
#         (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).tf, autodiff=true)
#     @NLexpression(optsteps, summedtime, # add up time of all steps
#         sum(onestept(v[i],P[i],δ[i]) for i = 1:numsteps))
#     for i = 1:numsteps  # step dynamics
#         @NLconstraint(optsteps, vels[i]==onestepv(v[i],P[i],δ[i]))
#     end
#     @NLconstraint(optsteps, summedtime == totaltime) # total time
#
#     if boundarywork
#         @objective(optsteps, Min, 1/2*(sum((P[i]^2 for i=1:numsteps))+v[1]^2-boundaryvels[1]^2)) # minimum pos work
#     else
#         @objective(optsteps, Min, 1/2*sum((P[i]^2 for i=1:numsteps))) # minimum pos work
#     end
#     optimize!(optsteps)
#     if termination_status(optsteps) == MOI.LOCALLY_SOLVED || termination_status(optsteps) == MOI.OPTIMAL
#         optimal_solution = (vms=value.(v), Ps=value.(P))
#     else
#         error("The model was not solved correctly.")
#         println(termination_status(optsteps))
#     end
#     return multistep(Walk(w,vm=value(v[1])), value.(P), δ, value(v[1]), boundaryvels,
#         extracost = boundarywork ? 1/2*(value(v[1])^2 - boundaryvels[1]^2) : 0) #, optimal_solution
# end



end # Module
