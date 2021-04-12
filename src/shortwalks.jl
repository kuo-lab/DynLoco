using DynLoco
using Plots
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

## Short walks of different distances
# Take walks of varying distances, and show how the optimal trajectory is to have a bell-shaped
# velocity profile, with peak speed that increases with distance up to about 12 steps.
# The cost function is total work, plus a linear cost of time with coefficient ctime.
wstar4 = findgait(WalkRW2l(α=0.35), target=:speed=>0.3, varying=:P)
ctime = 0.012 # cost of time, to encourage hurrying
tchange = 2
p = plot()
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
results = Array{MultiStepResults,1}(undef,0) # store each optimization result here
for (i,nsteps) in enumerate(walksteps)
    result = optwalktime(wstar4, nsteps, ctime=ctime) # optimize with a cost of time
    plotvees!(result, tchange=tchange, usespline=true, color=i, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(results, result) # add this optimization to results array
end
Plots.display(p) # instantaneous speed vs. distance profiles
savefig(p, "shortwalks.svg")
savefig(p, "shortwalks.pdf")

## Temp stuff, a changing step length; seems to work!
wstar4s = findgait(WalkRW2ls(α=0.35,safety=true), target=:speed=>0.3, varying=:P)
ctime = 0.012 # cost of time, to encourage hurrying
tchange = 2
p = plot()
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
resultss = Array{MultiStepResults,1}(undef,0) # store each optimization result here
for (i,nsteps) in enumerate(walksteps)
    result = optwalktime(wstar4s, nsteps, ctime=ctime) # optimize with a cost of time
    plotvees!(result, tchange=tchange, usespline=true, color=i, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(resultss, result) # add this optimization to results array
end
Plots.display(p) # instantaneous speed vs. distance profiles


## Try plotting stride speeds

function stridespeeds(steps)
    steptimes = [steps.tf; tchange]
    stepdistances = [steps.steplength; 0]
    return cumsum([0;steptimes],dims=1), [0; stepdistances./steptimes]
end

# splinify
using Dierckx
function splinify(t, v; ramporder = 1, thresholds = (0.05,0.95))
    if length(v) < 2 # enough points to make splines from v alone
        error("not enough points")
    end
    spline = Spline1D(t, v; k=ramporder)
    twhole = range(t[1], t[end], length=50)
    vwhole = spline.(twhole)
    tsteady = twhole[vwhole .>= thresholds[2]*maximum(vwhole)]
    taccel = twhole[thresholds[1]*maximum(vwhole) .<= vwhole .< thresholds[2]*maximum(vwhole)]
    tspeedup = taccel[taccel .< 0.5*twhole[end]]
    tslowdown = taccel[taccel .>= 0.5*twhole[end]]
    println(tspeedup[end]-tspeedup[1], " ", taccel[end]-taccel[1], " ", tslowdown[end]-tslowdown[1])
    return tspeedup[end]-tspeedup[1], taccel[end]-taccel[1], tslowdown[end]-tslowdown[1]
end

# plot stride speeeds vs distance
p = plot()
Tdata = Array{Float64}(undef, length(results), 3)
for i in 1:length(results)
    plot!(stridespeeds(results[i].steps)..., ylims=(0,Inf))
    Tspeedup, Tsteady, Tslowdown = splinify(stridespeeds(results[i].steps)...)
    Tdata[i,:] .= splinify(stridespeeds(results[i].steps)...; thresholds=(0.03,0.98))
end
Plots.display(p)
distances = [sum(result.steps.steplength) for result in results]
plot(distances, Tdata./[r.totaltime+2*tchange for r in results])




## Short walks: Peak speed vs. distances
peakspeeds = [maximum(result.steps.vm) for result in results]     # mid-stance speeds
peakspeeds = [maximum(stridespeeds(r.steps)[2]) for r in results] # stride speeds
distances = [sum(result.steps.steplength) for result in results]
p1 = plot(distances, peakspeeds, xlabel="Distance", ylabel="Peak speed", xlims=(0,Inf), ylims=(0,Inf))
p2 = plot(walksteps, peakspeeds, xlabel="# of steps", ylabel="Peak speed")
plot(p1, p2, layout = (1,2), legend=false)
#savefig("peakshortwalks.svg")
savefig("peakshortwalks.pdf")

## Short walks: Time to walk a distance
# A fairly linear increase in time to walk a distance, but with a slight curved toe-in
timetowalk = [result.totaltime+tchange for result in results]
plot(distances, timetowalk, xlims=(0,Inf), ylims=(0,Inf),
    xguide="Distance", yguide="Time", title="Time to walk a distance", label=nothing)
savefig("durationdistance.pdf")

## Short walks: Varying ctimes to demonstrate self-similarity
wstar4 = findgait(WalkRW2l(α=0.35), target=:speed=>0.3, varying=:P)
ctimes = range(0.006, 0.06, length=6)
tchange = 2
layout = @layout[ a{0.85w} grid(6,1)]
p = plot(;layout)
peaks = zeros(length(walksteps),length(ctimes))
durations = similar(peaks)
walksteps = [2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
results = Array{MultiStepResults,2}(undef,(length(walksteps),length(ctimes))) # store each optimization result here
for (j,ctime) in enumerate(ctimes)
    for (i,nsteps) in enumerate(walksteps)
        result = optwalktime(wstar4, nsteps, ctime=ctime) # optimize with a cost of time
        #peaks[i,j] = maximum(result.steps.vm)
        peaks[i,j] = maximum(stridespeeds(result.steps)[2])
        durations[i,j] = result.totaltime
        results[i,j] = result
    end
end
# after the fact, let's plot them all on top of each other
# using the ctime=0.02 result as the basis
tbase = durations[end,2]
vbase = peaks[end,2]
for (j, ctime) in enumerate(ctimes)
    for (i,nsteps) in enumerate(walksteps)
        result = results[i,j]

        plotvees!(result, tchange=tchange, color=i, usespline=:false, markersize=2, subplot=j+1,
            xticks = [20,40], yticks=[0.2,0.4,0.6],xguide="",yguide="",tickfontsize=4,
            xlims=(0,maximum(durations)+2tchange), ylims=(0,maximum(peaks))) # plot instantaneous speed vs. time
        plotvees!(result, tchange=tchange, color=i, usespline=:false, markersize=2, tscale = tbase/(durations[end,j]), 
            vscale = vbase/peaks[end,j],subplot=1)
    end
end
Plots.display(p) 
println("Durations of a factor of ", (durations[end,1]+2tchange)/(durations[end,end]+2tchange))
println("Peak speeds over a range of ", peaks[end,end]/peaks[end,1])
println("  about ", peaks[end,1]*sqrt(9.81)," to ", peaks[end,end]*sqrt(9.81), "m/s")
savefig("selfsimilarity.pdf")
savefig("selfsimilarity.svg")
# GR no fonts, doesn't do eps
# plotlyJS did export fonts, not necessarily the right one
# pyplot doesn't preserve fonts, but does export eps





## Short walks: Up and down slopes
# Compare walking uphill, downhill, and level, for a fixed number of steps, and including
# optimized time. The results show uphill walking is skewed to fast speed-up at beginning,
# slow coasting toward end. Downhill is skewed for slow speed-up aided by gravity, and
# abrupt end.
# `optwalktime` is an optimization for a given number of steps, with a `ctime` cost.

# Make plots comparing up, level, down for various numbers of steps
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.4, varying=:P)
myslope = 0.08
p = plot(layout=(3,1))
for (i,nsteps) = enumerate([5, 10, 15])
    resultlevel = optwalktime(wstar4s, nsteps, ctime = ctime)
    plotvees!(p[1],resultlevel, tchange=tchange, title="Level", rampuporder=1) # special function to include ramp-up in speed

    # walk up a 10% slope
    resultup = optwalktime(wstar4s, nsteps, ctime = ctime, δs=fill(myslope, nsteps))
    plotvees!(p[2],resultup, tchange=tchange, title="Up", rampuporder=1)

    # Walk down a slope
    resultdown = optwalktime(wstar4s, nsteps, ctime = ctime, δs=fill(-myslope, nsteps))
    plotvees!(p[3],resultdown, tchange=1, title="Down", rampuporder=1)
end
Plots.display(p)

## Short walks: Varying slope angles
# The steeper, the more skewed the speed profile

# Make plots comparing up, level, down for various numbers of steps
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.4, varying=:P)
myslopes = 0:0.02:0.08
p = plot(layout=(2,1))
nsteps = 6
startv = wstar4.vm
for slope in myslopes
    # walk up a slope
    resultup = optwalktime(wstar4s, nsteps, ctime = ctime, δs=fill(slope, nsteps))
    plotvees!(p[1],resultup, tchange=1, title="Up", rampuporder=1)
    startv = [resultup.vm0;resultup.steps.vm]

    # Walk down a slope
    resultdown = optwalktime(wstar4s, nsteps, ctime = ctime, δs=fill(-slope, nsteps))
    plotvees!(p[2],resultdown, tchange=1, title="Down", rampuporder=1)
end
Plots.display(p)

## Short walks: Varying step lengths
# Shorter steps will yield a more plateaued speed profile.
# Longer steps will reach more of a rounded profile with bigger range of speeds.
# (This optimization allows for faster speeds at short steps, due to a lack of
# swing leg cost. Humans will probably not walk faster with short steps.)
stepfreq = onestep(wstar4).stepfrequency
αs = [0.2, 0.3, 0.4, 0.5, 0.6]
nominalsteps = 6
totaldistance = 6*onestep(wstar4).steplength # walk a similar distance for all
steplengths = 2*wstar4.L .* sin.(αs)
p = plot()
results = Array{MultiStepResults,1}(undef,0) # store each optimization result here
for (i,α) in enumerate(αs)
    nsteps = Int(round(totaldistance/steplengths[i]))
    w = findgait(WalkRW2l(wstar4,α=α,safety=true), target=:stepfrequency=>stepfreq, varying=:P)
    result = optwalktime(w, nsteps, ctime=ctime) # optimize with a cost of time
    plotvees!(result, tchange=1, color=i, rampuporder=1) # plot instantaneous speed vs. time
    push!(results, result) # add this optimization to results array
end
Plots.display(p) # instantaneous speed vs. distance profiles
# savefig("steplengthsshortwalks.svg")


## Short walks: Compare trapezoid cruising against short walk
# Another way to walk a short distance with cruise speed, "trapezoid" velocity
# Start from a boundaryvel, do boundary work to get a certain speed
# then continue with a ramp up in speed, then stay at a cruise speed
# and then ramp down and hit another boundaryvel
# v[1] is the mid-stance velocity where you have say Naccel steps to get to cruise
#  speed up    cruise cruise cruise cruise   slow down
# 0 v[1] v[2]  v[Naccel+1] to v[Naccel + N]  v[end-2] v[end-1] v[end]
using JuMP, Ipopt
boundaryvels = (0.,0.) # actually the code below all assumes 0 boundary speed

Ncruise = 5 # how many steps of cruising
Naccel = 0  # how many steps of start-up (and again for end)
Nsteps = Naccel*2 + Ncruise
vcruise = 0.45 # cruising speed, always fixed
# linear increase deltavel = vcruise/(Naccel+1)
# constraint v[i=1..Naccel] = deltavel*i
deltavel = vcruise / (Naccel+1) # assuming starting from zero speed
velstart = [deltavel*i for i in 1:Naccel] # ramp up in speed
velcruise = [vcruise for i in 1:Ncruise]
velend = [deltavel*(Naccel+1-i) for i in 1:Naccel]
vels = [velstart; velcruise; velend]
# constraint v[Naccel+1 ... Naccel+N] == vcruise
# constraint v[Naccel+N+(1..Naccel)] == deltavel*(Naccel+1-i)
# solve for P that produces it, time will be an outcome
w = WalkRW2l(wstar4, safety=true)
optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
@variable(optsteps, P[1:Nsteps]>=0, start=w.P) # JuMP variables P

# Constraints
register(optsteps, :onestepv, 3, # velocity after a step
    (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).vm, autodiff=true) # output vm
for i = 1:Nsteps-1  # step dynamics
    @NLconstraint(optsteps, vels[i+1]==onestepv(vels[i],P[i],0.)) # put delta here
end

# leave out the objective, because we are fully prescribing everything with constraints
optimize!(optsteps)
if termination_status(optsteps) == MOI.LOCALLY_SOLVED || termination_status(optsteps) == MOI.OPTIMAL
    optimal_solution = Ps=value.(P)
else
    error("The model was not solved correctly.")
    println(termination_status(optsteps))
end
trapezoidresults=multistep(WalkRW2l(w,vm=vels[1]), Ps=optimal_solution, boundaryvels=(0.,0.),
    extracost = 1/2*(vels[1]^2 - boundaryvels[1]^2) )

# verify with multistep
multistepplot(trapezoidresults,plotwork=true, label="square")
savefig("trapezoidshortwalks.svg")

# If you just want to do square wave in speed, it costs a lot of initial push-off
# so let's compare with walking the same number of steps and same amount of time
optresults=optwalk(w, Ncruise, boundaryvels=(0,0),totaltime=trapezoidresults.totaltime  )
println("trapezoid cost = ", trapezoidresults.totalcost, "   optimal cost = ", optresults.totalcost)
# It's definitely more expensive to use the square wave
multistepplot!(optresults,plotwork=true,label="optimal")
savefig("trapezoidcomparisonshortwalks.svg")

# nice speed profile comparison
plotvees(trapezoidresults, tchange=1, rampuporder=1)
plotvees!(optresults, tchange=1, rampuporder=1)

# nice work comparison

plot([0.5; 1:Nsteps; Nsteps+0.5], [1/2*optresults.vm0^2; optresults.steps.Pwork; NaN],markershape=:circle)
plot!([0.5; 1:Nsteps; Nsteps+0.5], [1/2*trapezoidresults.vm0^2; trapezoidresults.steps.Pwork; NaN],markershape=:circle)


## Trajectory for short walks


## Brachistokuo Ramp
# Optimal slope and walk with ramp
# Compare walking ramp and flat in same amount of time, for three different times
wstar = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.4, varying=:P)
N = 6
walktime = N * onestep(wstar).tf *0.82 # meant to be a brisk walk
walkdistance = N * onestep(wstar).steplength
rampresult = optwalkslope(wstar, N, boundaryvels = (0., 0.), symmetric = true,
    totaltime = walktime)
p = multistepplot(rampresult; plotwork=true, label="ramp")
println("ramp total cost = ", rampresult.totalcost)
flatresult = optwalk(wstar, N, boundaryvels = (0., 0.),
    totaltime = rampresult.totaltime, δ = zeros(6))
println("flat total cost = ", flatresult.totalcost)
multistepplot!(flatresult; plotwork=true, label="flat")
# optionally, try a reversed ramp and see if it's higher cost still
#concaveresult = optwalk(wstar, 6, boundaryvels = (0.,0.), boundarywork=true,
#    totaltime = rampresult.totaltime, δ = -rampresult.δangles)
#multistepplot!(concaveresult; plotwork=true)

## Brachistokuo ramp: Compute the cost for different speeds
walktimes = (0.8:0.05:1.2) * N*onestep(wstar).tf
rampresults = Array{MultiStepResults,1}(undef, length(walktimes))
flatresults = Array{MultiStepResults,1}(undef, length(walktimes))
for (i,walktime) in enumerate(walktimes)
    rampresults[i] = optwalk(wstar, N, boundaryvels = (0.,0.),
        totaltime = walktime, δ = rampresult.δangles)
    flatresults[i] = optwalk(wstar, N, boundaryvels = (0.,0.),
        totaltime = walktime, δ = zeros(N))
end
# plot totalcost vs average speed
plot(walkdistance ./ walktimes, [getfield.(rampresults, :totalcost), getfield.(flatresults, :totalcost)],
    xlabel="Average speed", ylabel="Total Work", labels=["Ramp" "Flat"])
# savefig("rampvsflat.pdf")

## Brachistokuo ramp: Plot the ramp to scale
plot(onestep(wstar).steplength .* (0:6),cumsum(tan.([0;rampresult.δangles]).*onestep(wstar).steplength),
    aspect_ratio=1)
sl = onestep(wstar).steplength
plot(sl .* [0; cumsum(cos.(rampresult.δangles))],[0; cumsum(sin.(rampresult.δangles))].*sl,
    aspect_ratio=1)
# savefig("ramptoscale.pdf")


includet("drawingtrial.jl")
using .DynLocoGraphics



drawmodel(wstar,rampresult.steps[1])


p = plot(onestep(wstar).steplength .* (0:6),cumsum(tan.([0;rampresult.δangles]).*onestep(wstar).steplength),
    aspect_ratio=1)


# Use this with our model
sl = onestep(wstar).steplength
floorx = sl .* [0; cumsum(cos.(rampresult.δangles))]
floory = sl .* [0; cumsum(sin.(rampresult.δangles))]
p = plot(floorx, floory, aspect_ratio=1, showaxis=false, label=false)
for i = 1:length(rampresult.steps)
    drawmodel!(p, wstar, rampresult.steps[i],(floorx[i],floory[i]),scalev=0.6)
end
display(p)

# actuallly i'm not sure if we get a passive gait going downhill
findgait(WalkRW2l(wstar4,P=0), target=:speed=>0.3,varying=:γ)


## A downhill gait
wd=WalkRW2l(wstar4,P=0,γ=0.06) # this seems to work
wdstar = findgait(wd, target=:speed=>0.2, varying=:γ)
onestep(wd) # γ=0.0538776, vm = 0.14129694390484976
msr=multistep(WalkRW2l(wdstar,P=0,γ=0), P=zeros(14),δangles=-0.0538776800329528*ones(14))
# yeah, something is wrong with going downhill, this doesn't match with
# the wdstar gait, although it may have to do with how mid-stance is defined


## MPC attempt. Do a short walk, and do a finite-horizon MPC after each step
# Optimize a short walk on level ground in resultlevel, with a time objective
# and then do another optimization for a fixed (but same) amount of time (resultsettime)
# and then do a finite (fixed) horizon MPC to re-do that optimization
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
tchange = 2
myslope = 0.08
ctime = 0.02
nsteps = 10
resultlevel = optwalktime(wstar4s, nsteps, boundarywork = true, boundaryvels=(0,0), ctime = ctime,safety=false)
plotvees(resultlevel, tchange=tchange, title="Level", usespline=false,rampuporder=1, speedtype=:midstance) # special function to include ramp-up in speed

# check whether you can optimize for same steps but with time constrained
resultsettime=optwalk(wstar4s, nsteps, boundarywork=true, boundaryvels=(0,0),totaltime=resultlevel.totaltime  )
plotvees!(resultsettime, tchange=tchange, usespline=false, speedtype=:midstance)


# Start MPC
remainingtime = resultsettime.totaltime
remainingsteps = nsteps
currentstep = resultsettime.steps[1]
currentvm0 = resultsettime.steps[1].vm0
bcwork = 1/2*(currentstep.vm0^2 - 0^2) # applied boundary impulse, now ready to step
mysteps = Vector{StepResults}(undef,nsteps)
for i in 1:nsteps # finite horizon after each i'th step
    println("i = $i")
    # optimize starting from the most recent vm0
    nextmsr = optwalk(wstar4s, remainingsteps, boundarywork=(false,true), boundaryvels=(currentvm0,0), totaltime=remainingtime)
    # take a step (optional, should match what was optimized)
    nextstep = onestep(wstar4s, vm=currentvm0, P=currentstep.P, δangle=currentstep.δ)
    nextstep = nextmsr.steps[1] # or use next step: StepResults(nextstep...)
    mysteps[i] = nextstep 
    println("nextstepvm = ", nextstep.vm, "  result.vm = ", resultsettime.steps[i].vm)
    @assert isapprox(nextstep.vm, resultsettime.steps[i].vm, atol=1e-4) # check whether the steps agree
    plot!(cumsum([tchange+resultsettime.totaltime-remainingtime;nextmsr.steps.tf]),
        [nextstep.vm0;nextmsr.steps.vm],show=true)
    remainingsteps = remainingsteps - 1
    remainingtime = remainingtime - nextstep.tf
    # set up for the next one
    currentvm0 = nextstep.vm
end
# okay, this works; each step is a finite horizon MPC with a shortening horizon

## Finite horizon MPC walking over a bump
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
nsteps = 15
δs = zeros(nsteps); δs[Int((nsteps+1)/2)] = 0.05
nominalmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δs)
plotvees(nominalmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)

remainingtime = nominalmsr.totaltime
remainingsteps = nsteps
currentstep = nominalmsr.steps[1]
currentvm0 = nominalmsr.steps[1].vm0
mpcsteps = Vector{StepResults}(undef,nsteps)
for i in 1:nsteps-1 # finite horizon after each i'th step; don't optimize the last step
    println("i = $i")
    # optimize starting from the most recent vm0
    nextmsr = optwalk(wstar4s, remainingsteps, boundarywork=(false,false), boundaryvels=(currentvm0,wstar4s.vm), totaltime=remainingtime, δs=δs[i:nsteps])
    nextstep = nextmsr.steps[1] 
    mpcsteps[i] = nextstep 
    println("nextstepvm = ", nextstep.vm, "  result.vm = ", nominalmsr.steps[i].vm)
    @assert isapprox(nextstep.vm, nominalmsr.steps[i].vm, atol=1e-4) # check whether the steps agree
    plot!(cumsum([nominalmsr.totaltime-remainingtime;nextmsr.steps.tf]),
        [nextstep.vm0;nextmsr.steps.vm],show=true)
    remainingsteps = remainingsteps - 1
    remainingtime = remainingtime - nextstep.tf
    # set up for the next one
    currentvm0 = nextstep.vm
end
# last step can't be solved easily, because we have one push-off to satisfy both
# the last time and the last vm; so we don't bother with it

## Receding horizon MPC walking over a bump
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
nsteps = 15
δs = zeros(nsteps); δs[Int((nsteps+1)/2)] = 0.05
nominalmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δs)
plotvees(nominalmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)
nhorizon = 12

elapsedtime = 0
modeltime = 0
currentstep = nominalmsr.steps[1]
currentvm0 = nominalmsr.steps[1].vm0
mpcsteps = Vector{StepResults}(undef,nsteps)
upcomingδs = zeros(nhorizon)
for i in 1:nsteps-1 # receding horzion; don't optimize the last step
#i=1
    println("i = $i")
    # optimize starting from the most recent vm0
    errortime = elapsedtime - modeltime # how far you're ahead of nominal model
    remainingtime = tfstar * nhorizon - errortime # make up for lost time
    if i+nhorizon-1 > nsteps # need to pad zeros
        # let's say your horizon is 6, and nsteps = 10, but
        # but we are on step i = 6, so we have to pad with zeros
        extrazerosteps = i+nhorizon-1-nsteps
        upcomingδs[1:nhorizon-extrazerosteps] = δs[i:nsteps]
        upcomingδs[nhorizon-extrazerosteps+1:nhorizon] .= 0
    else # plenty of steps left
        upcomingδs[1:nhorizon] .= δs[i:i+nhorizon-1] # unless i+nhorizon-1 exceeds nsteps, then we need to shrink horizon
    end # setting upcomingδs
    nextmsr = optwalk(wstar4s, nhorizon, boundarywork=(false,false), boundaryvels=(currentvm0,wstar4s.vm), totaltime=remainingtime, δs=upcomingδs)
    nextstep = nextmsr.steps[1] 
    mpcsteps[i] = nextstep 
    println("nextstepvm = ", nextstep.vm, "  result.vm = ", nominalmsr.steps[i].vm)
    #@assert isapprox(nextstep.vm, nominalmsr.steps[i].vm, atol=1e-4) # check whether the steps agree
    plot!(cumsum([elapsedtime;nextmsr.steps[1:end-1].tf]),
        [nextstep.vm0;nextmsr.steps.vm[1:end-1]],show=true)
    remainingtime = remainingtime - nextstep.tf
    # set up for the next one
    currentvm0 = nextstep.vm
    elapsedtime = elapsedtime + nextstep.tf
    modeltime = modeltime + tfstar
end