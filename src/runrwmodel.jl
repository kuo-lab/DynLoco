
using DynLocoT
using Plots

wrw = Walk(P=0.15,vm=0.38)
wrwg = Walk(P=0.05,vm=0.38,γ=0.02)
wrw2 = Walk(wrw, α=0.35)
wstar = findgait(wrw, target=:speed=>0.4, varying=:P)
onestep(wrw)
onestep(wrw, α=0.25, P=0.2)

findgait(wrw)
findgait(wrw, P=0.2, target=:speed=>0.4, varying=:P)
findgait(wrw; vm=0.35, P=0.15, α=0.35, :γ=>0.15)
findgait(wrw, target=:speed=>0.45, varying=:γ, P=0) # gravity only, no push-off
findgait(wrw, target=:vm=>0.4, varying=:P) # use vm as target
findgait(wrw, α=0.32, target=(:speed=>0.4,), varying=(:P,), vm=0.25, P=0.2) # use tuple of targets

# test safe step with too littl momentum
onestep(wrw, P=0., vm=0.1, safety=true) # should fail if safety=false

mystuff = multistep(wstar, [0.15 0.17 0.18])
multistepplot(mystuff)

# probably deprecate
wnew = findlimitcycle(wrw, target=:speed=>0.4, varying=:P)
onestep(wnew)
onestep(wnew).speed
onestep(wrwg)
onestep(wrwg,P=0.1)
wstar = findlimitcycle(wrw, target=:speed=>0.4, varying=:P)
findlimitcycle(wrw, target=:speed=>0.4, varying=:P, vm=0.2)
findlimitcycle(wrw, vm=0.2, target=:speed=>0.4, varying=:P)

# walk downhill
onestep(wstar, δangle=-0.05)

myout = multistep(wstar, 1.5*(wstar.P)*ones(15))
plot(getfield.(myout.steps,:vm))
multistepplot(myout)
## Optimize multiple steps for a short walk
# 5 steps of level ground, starting and ending from steady speed
using JuMP, Ipopt
nsteps = 5
optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
@variable(optsteps, P[1:nsteps], start=wstar.P) # JuMP variables P
#@variable(optsteps, δ[1:nsteps], start=0.)
δs = zeros(nsteps) # set bumps to zero
# constraints: start
@variable(optsteps, v[1:nsteps], start=wstar.vm)
@constraint(optsteps, v[1] == 0.3)
@constraint(optsteps, v[nsteps] == 0.3)
@objective(optsteps, Min, sum((P[i]^2 for i=1:nsteps)))
register(optsteps, :onestepv, 2, (v,P)->onestep(wstar,P=P,vm=v).vm, autodiff=true) # input P, output vm
for i = 1:nsteps-1
    @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i]))
end
optimize!(optsteps)
plot(value.(v))
plot(value.(P)) # push-off impulses

## Short walk in a certain amount of time
using JuMP, Ipopt
nsteps = 5
optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
@variable(optsteps, P[1:nsteps]>=0, start=wstar.P) # JuMP variables P
#@variable(optsteps, δ[1:nsteps], start=0.)
δs = zeros(nsteps) # set bumps to zero
# constraints: start
@variable(optsteps, v[1:nsteps+1], start=wstar.vm)
@constraint(optsteps, v[1] == 0.25)
@constraint(optsteps, v[nsteps+1] == 0.25)
@objective(optsteps, Min, sum((P[i]^2 for i=1:nsteps)))
register(optsteps, :onestepv, 2, (v,P)->onestep(wstar,P=P,vm=v).vm, autodiff=true) # input P, output vm
register(optsteps, :onestept, 2, (v,P)->onestep(wstar,P=P,vm=v).tf, autodiff=true)
@NLexpression(optsteps, steptime[i=1:nsteps], onestept(v[i],P[i]))
@NLexpression(optsteps, totaltime, sum(onestept(v[i],P[i]) for i = 1:nsteps))
#@NLexpression(optsteps, collisionimpulse, [use onestep to compute collision including negative collisions])
#@NLconstraint(optsteps, collisionimpulses[i=1:nsteps]>=0)
for i = 1:nsteps
    @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i]))
end
@NLconstraint(optsteps, totaltime == nsteps*onestep(wstar).tf)
optimize!(optsteps)
plot(value.(v); xlabel="step", ylabel="vm")
plot(value.(P), xlabel="step", ylabel="P") # push-off impulses
result = multistep(Walk(wstar,vm=0.25), value.(P), δs)
plot([0.25;getfield.(result.steps,:vm)])
## Short walk in a certain amount of time and slopes, plus some initial work
using JuMP, Ipopt
nsteps = 5
optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
@variable(optsteps, P[1:nsteps]>=0, start=wstar.P) # JuMP variables P
@variable(optsteps, δ[1:nsteps], start=0.)
#δs = zeros(nsteps) # set bumps to zero
# constraints: start
@variable(optsteps, v[1:nsteps+1]>=0, start=wstar.vm)
#@constraint(optsteps, v[1] == 0.05)
#@constraint(optsteps, v[nsteps+1] == 0.2)
#for i = 1:nsteps
#@constraint(optsteps, δ[i] == 0.)
#end
@constraint(optsteps, sum(δ[i] for i =1:nsteps) == 0.)  # zero height gain
@objective(optsteps, Min, sum((P[i]^2 for i=1:nsteps))+v[1]^2) # minimum work
register(optsteps, :onestepv, 3, (v,P,δ)->onestep(wstar,P=P,vm=v,δangle=δ).vm, autodiff=true) # input P, output vm
register(optsteps, :onestept, 3, (v,P,δ)->onestep(wstar,P=P,vm=v,δangle=δ).tf, autodiff=true)
#@NLexpression(optsteps, steptime[i=1:nsteps], onestept(v[i],P[i],δ[i]))
@NLexpression(optsteps, totaltime, sum(onestept(v[i],P[i],δ[i]) for i = 1:nsteps))
for i = 1:nsteps
    @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i],δ[i]))
end
@NLconstraint(optsteps, totaltime == nsteps*onestep(wstar).tf)
optimize!(optsteps)
plot(value.(v))
plot(value.(P)) # push-off impulses
plot([0.; cumsum(value.(δ))])
myout = optwalk(wstar)

# level walking
myout = optwalk(wstar, 6, boundaryvels = (0.1,0.1))
multistepplot(myout)
myout = optwalk(wstar, 6, boundaryvels = (0.,0.), δ=[-0.1,0.,0.,0.,0.,0.1])
multistepplot(myout)
myout1 = optwalkslope(wstar, 5, boundaryvels = (0., 0.), symmetric = true)
multistepplot(myout1)
myout2 = optwalk(wstar, 5, boundaryvels=(0.,0.), totaltime=myout1.totaltime,
    δ=myout1.δangles) # shoot doesn't work
multistepplot(myout2)


## Short walk with time as an objective not constraints
using JuMP, Ipopt
nsteps = 10
wstar4 = findgait(wrw, α=0.4, target=:speed=>0.4, varying=:P)
optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
@variable(optsteps, P[1:nsteps]>=0, start=wstar4.P) # JuMP variables P
#@variable(optsteps, δ[1:nsteps], start=0.)
δs = zeros(nsteps) # set bumps to zero
# constraints: start
@variable(optsteps, v[1:nsteps+1], start=wstar4.vm)
@constraint(optsteps, v[1] == 0.2)
@constraint(optsteps, v[nsteps+1] == 0.2)
register(optsteps, :onestepv, 2, (v,P)->onestep(wstar4,P=P,vm=v).vm, autodiff=true) # input P, output vm
register(optsteps, :onestept, 2, (v,P)->onestep(wstar4,P=P,vm=v).tf, autodiff=true)
@NLexpression(optsteps, steptime[i=1:nsteps], onestept(v[i],P[i]))
@NLexpression(optsteps, totaltime, sum(onestept(v[i],P[i]) for i = 1:nsteps))
#@NLexpression(optsteps, collisionimpulse, [use onestep to compute collision including negative collisions])
#@NLconstraint(optsteps, collisionimpulses[i=1:nsteps]>=0)
for i = 1:nsteps
    @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i]))
end
#@NLconstraint(optsteps, totaltime == nsteps*onestep(wstar4).tf)
ctime = 0.05
@NLobjective(optsteps, Min, sum((P[i]^2 for i=1:nsteps)) + ctime*totaltime)
optimize!(optsteps)
plot(value.(v); xlabel="step", ylabel="vm")
plot(value.(P), xlabel="step", ylabel="P") # push-off impulses
result = multistep(Walk(wstar4,vm=0.25), value.(P), δs)
plot([0.25;getfield.(result.steps,:vm)])
multistepplot(result)
# with coefficient of 0.05, model seems to prefer taking 5 steps to speed up
# over 50 steps
# I set steps to be longer 0.4, and it takes about 3 steps to
# speed up
# try plotting vs. time

## Short walks with time as an objective not constraints
using JuMP, Ipopt
wstar4 = findgait(wrw, α=0.4, target=:speed=>0.4, varying=:P)
#p = plot(layout=(3,1))
p = plot()
boundaryvels = (0.,0.)
for nsteps in [1, 2, 3, 4, 5, 7, 10, 15, 20]
    optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>1))
    @variable(optsteps, P[1:nsteps]>=0, start=wstar4.P) # JuMP variables P
    δs = zeros(nsteps) # set bumps to zero
    # constraints: start
    @variable(optsteps, v[1:nsteps+1], start=wstar4.vm)
    #@constraint(optsteps, v[1] == boundaryvel)
    #@constraint(optsteps, v[nsteps+1] == boundaryvel)
    register(optsteps, :onestepv, 2, (v,P)->onestep(wstar4,P=P,vm=v,safety=true).vm, autodiff=true) # input P, output vm
    register(optsteps, :onestept, 2, (v,P)->onestep(wstar4,P=P,vm=v,safety=true).tf, autodiff=true)
    @NLexpression(optsteps, steptime[i=1:nsteps], onestept(v[i],P[i]))
    @NLexpression(optsteps, totaltime, sum(onestept(v[i],P[i]) for i = 1:nsteps))
    for i = 1:nsteps # collocation points
        @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i]))
    end
    ctime = 0.04 # 0.07 looks pretty good
    @NLobjective(optsteps, Min, sum((P[i]^2 for i=1:nsteps)) + v[1]^2 - boundaryvels[1]^2 +
        ctime*totaltime)
    optimize!(optsteps)
    global result = multistep(Walk(wstar4,vm=value(v[1])), value.(P), δs, value(v[1]), boundaryvels)
    #Plots.display(multistepplot!(p,result))
    #Plots.display(plot!(cumsum([0;result.steps.tf]), [boundaryvels[1];result.steps.vm]))
    plotsmoothvees(result, tchange=4.)
    #plot!(cumsum([0; result.steps.tf]),[0; result.steps.speed])
end
Plots.display(p)
# a low time cost (min work) causes almost constant (slow) speeds
# and the more time costs, the more the top speed
# TODO: My boundary vels are a bit weird

# Plot speed vs. time, accumulating time
plot(cumsum([0;result.steps.tf]), [boundaryvels[1];result.steps.vm])

# Make a plot where it takes time to get up to speed
function plotsmoothvees(msr::MultiStepResults; tchange = 2)
    v = [msr.vm0; msr.steps.vm]
    #v = [msr.vm0; msr.steps.speed]
    n = length(msr.steps)
    times = cumsum([0; result.steps.tf])
    t0 = range(0, tchange, length=10)
    vstart = v[1]*(t0/tchange).^2
    vend = v[n+1]*(1 .- t0/tchange).^2
    alltimes = [t0; times .+ t0[end]; t0 .+ t0[end] .+ times[end]]
    allvees = [vstart; v; vend]
    plot!(p,alltimes, allvees)
end

# Another way to do it is with exponentials

# filter it
using DSP
using DSP.Filters


# With logshave we're able to use a boundaryvel as low as 0.16, but no lower

#
x = range(-0.001,0.1, length=100)
plot(x,logshave.(x))
