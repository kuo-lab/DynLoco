using DynLoco
using Plots

## Short walks of different distances
wstar4 = findgait(Walk(α=0.4), target=:speed=>0.4, varying=:P)
p = plot()
results = Array{MultiStepResults,1}(undef,0)
for (i,nsteps) in enumerate([1, 2, 3, 4, 5, 7, 10, 15, 20])
    result = optwalktime(wstar4, nsteps, ctime=0.05)
    plotvees!(result, tchange=3, color=i)
    push!(results, result)
end
Plots.display(p)

## Peak speed vs. distances
peakspeeds = [maximum(result.steps.vm) for result in results]
distances = [sum(result.steps.steplength) for result in results]
plot(distances, peakspeeds, xlabel="Distance", ylabel="Peak speed", xlims=(0,Inf), ylims=(0,Inf))

## Time to walk a distance
# A fairly linear increase in time to walk a distance, but with a slight curved toe-in
timetowalk = [result.totaltime+4 for result in results]
plot(distances, timetowalk, xlims=(0,Inf), ylims=(0,Inf),
    xlabel="Distance", ylabel="Time")

## Short walks on different slopes
# Compare walking uphill, downhill, and level, for a fixed number of steps, and including
# optimized time. The results show uphill walking is skewed to fast speed-up at beginning,
# slow coasting toward end. Downhill is skewed for slow speed-up aided by gravity, and
# abrupt end.

myslope = 0.08
wstar4 = findgait(Walk(α=0.4), target=:speed=>0.4, varying=:P)
p = plot()
for (i,nsteps) = enumerate([5, 10, 15])
result = optwalktime(wstar4, nsteps, ctime = 0.025)
plotvees!(result, tchange=3)

# walk up a 10% slope
result = optwalktime(wstar4, nsteps, ctime = 0.025, δs=fill(myslope, nsteps))
plotvees!(result, tchange=3)

# Walk down a slope
result = optwalktime(wstar4, nsteps, ctime = 0.025, δs=fill(-myslope, nsteps))
plotvees!(result, tchange=3)
end
Plots.display(p)

## Optimal slope and walk
# Compare walking ramp and flat in same amount of time, for three different times
wstar = findgait(Walk(wrw,safety=true), target=:speed=>0.4, varying=:P)
N = 6
walktime = N * onestep(wstar).tf
walkdistance = N * onestep(wstar).steplength
rampresult = optwalkslope(wstar, N, boundaryvels = (0., 0.), boundarywork= true, symmetric = true,
    totaltime = walktime)
multistepplot(rampresult; plotwork=true)
flatresult = optwalk(wstar, N, boundaryvels = (0., 0.), boundarywork=true,
    totaltime = rampresult.totaltime, δ = zeros(6))
multistepplot!(flatresult; plotwork=true)
# optionally, try a reversed ramp and see if it's higher cost still
#concaveresult = optwalk(wstar, 6, boundaryvels = (0.,0.), boundarywork=true,
#    totaltime = rampresult.totaltime, δ = -rampresult.δangles)
#multistepplot!(concaveresult; plotwork=true)

## Compute the cost for different speeds
walktimes = (1:0.05:1.5) * N*onestep(wstar).tf
rampresults = Array{MultiStepResults,1}(undef, length(walktimes))
flatresults = Array{MultiStepResults,1}(undef, length(walktimes))
for (i,walktime) in enumerate(walktimes)
    rampresults[i] = optwalk(wstar, N, boundaryvels = (0.,0.), boundarywork = true,
        totaltime = walktime, δ = rampresult.δangles)
    flatresults[i] = optwalk(wstar, N, boundaryvels = (0.,0.), boundarywork = true,
        totaltime = walktime, δ = zeros(N))
end
plot(walkdistance ./ walktimes, [getfield.(rampresults, :totalcost), getfield.(flatresults, :totalcost)],
    xlabel="Speed", ylabel="Work")


## Another way to walk a short distance with cruise speed
# start from a boundaryvel, do boundary work to get a certain speed
# then continue with a ramp up in speed, then stay at a cruise speed
# and then ramp down and hit another boundaryvel
# v[1] is the mid-stance velocity where you have say Naccel steps to get to cruise
#  speed up    cruise cruise cruise cruise   slow down
# 0 v[1] v[2]  v[Naccel+1] to v[Naccel + N]  v[end-2] v[end-1] v[end]

boundaryvels = (0.,0.)
boundarywork = true


Ncruise = 3
Naccel = 0
Nsteps = Naccel*2 + Ncruise
vcruise = 0.45 # cruising speed, always fixed
# linear increase deltavel = vcruise/(Naccel+1)
# constraint v[i=1..Naccel] = deltavel*i
deltavel = vcruise / (Naccel+1)
velstart = [deltavel*i for i in 1:Naccel]
velcruise = [vcruise for i in 1:Ncruise]
velend = [deltavel*(Naccel+1-i) for i in 1:Naccel]
vels = [velstart; velcruise; velend]
# constraint v[Naccel+1 ... Naccel+N] == vcruise
# constraint v[Naccel+N+(1..Naccel)] == deltavel*(Naccel+1-i)
# solve for P that produces it, time will be an outcome
w = wstar


    optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
    @variable(optsteps, P[1:Nsteps]>=0, start=w.P) # JuMP variables P

    if boundaryvels == nothing || isempty(boundaryvels)
        boundaryvels = (w.vm, w.vm) # default to given gait if nothing specified
    end

    if !boundarywork # no hip work at beginning or end; apply boundary velocity constraints
        @constraint(optsteps, v[1] == boundaryvels[1])
        @constraint(optsteps, v[Nsteps+1] == boundaryvels[2])
    end

    # Constraints
    # produce separate functions for speeds and step times
    register(optsteps, :onestepv, 3, # velocity after a step
        (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).vm, autodiff=true) # output vm
    #register(optsteps, :onestept, 3, # time after a step
    #    (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).tf, autodiff=true)
    #@NLconstraint(optsteps, vels[])
    for i = 1:Nsteps-1  # step dynamics
        @NLconstraint(optsteps, vels[i+1]==onestepv(vels[i],P[i],0.)) # put delta here
    end

# leave out the objective

    optimize!(optsteps)
    if termination_status(optsteps) == MOI.LOCALLY_SOLVED || termination_status(optsteps) == MOI.OPTIMAL
        optimal_solution = Ps=value.(P)
    else
        error("The model was not solved correctly.")
        println(termination_status(optsteps))
    end
    #return multistep(Walk(w,vm=value(v[1])), value.(P), δ, value(v[1]), boundaryvels,
    #    extracost = boundarywork ? 1/2*(value(v[1])^2 - boundaryvels[1]^2) : 0) #, optimal_solution
    squarewaveresults=multistep(Walk(w,vm=vels[1]), Ps=optimal_solution, boundaryvels=(0.,0.),
        extracost = boundarywork ? 1/2*(vels[1]^2 - boundaryvels[1]^2) : 0)

# verify with multistep
multistepplot(squarewaveresults,plotwork=true)

# If you just want to do square wave in speed, it costs a lot of initial push-off
# so let's compare with walking the same number of steps and same amount of time
optresults=optwalk(w, 3, boundaryvels=(0,0),boundarywork=true, totaltime=results.totaltime  )
multistepplot!(optresults,plotwork=true)
println("square wave cost = ", squarewaveresults.totalcost, "   optimal cost = ", optresults.totalcost)
# It's definitely more expensive to use the square wave
