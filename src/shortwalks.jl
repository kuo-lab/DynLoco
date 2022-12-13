using DynLoco
using Plots, Statistics
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

## Short walks of different distances
# Take walks of varying distances, and show how the optimal trajectory is to have a bell-shaped
# velocity profile, with peak speed that increases with distance up to about 12 steps.
# The cost function is total work, plus a linear cost of time with coefficient ctime.
wstar4 = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.3, varying=:P)
ctime = 0.015 # cost of time, to encourage hurrying
tchange = 1.75 # boundary condition time to get up to speed (arbitrary, excluded from optimization) 
p = plot()
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
results = Array{MultiStepResults,1}(undef,0) # store each optimization result here
for (i,nsteps) in enumerate(walksteps)
    result = optwalktime(wstar4, nsteps, ctime=ctime)#,negworkcost=0.2) # optimize with a cost of time
    plotvees!(result, tchange=tchange, usespline=false, color=i, speedtype=:shortwalks, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(results, result) # add this optimization to results array
end
Plots.display(p) # instantaneous speed vs. distance profiles
#savefig(p, "shortwalks.svg")
#savefig(p, "shortwalks.pdf")
# plus, currently adding variable step lengths below into a single plot


## Short walks of shorter and longer step lengths
wstar4 = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.3, varying=:P)
wstar43 = findgait(WalkRW2l(α=0.3,safety=true), target=:speed=>0.3, varying=:P)
wstar44 = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.3, varying=:P)
wsteplens = [wstar43, wstar4, wstar44]
ctimes = [0.015, 0.015, 0.015]
ctime = 0.015 # cost of time, to encourage hurrying
tchange = 1.75
p = plot()
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
results43 = Array{MultiStepResults,1}(undef,0) # store each optimization result here
results = Array{MultiStepResults,1}(undef,0) # store each optimization result here
results44 = Array{MultiStepResults,1}(undef,0) # store each optimization result here
peakspds = zeros(length(walksteps))
peakspds43 = zeros(length(walksteps))
peakspds44 = zeros(length(walksteps))
durations = zeros(length(walksteps))
durations43 = zeros(length(walksteps))
durations44 = zeros(length(walksteps))
for (i,nsteps) in enumerate(walksteps)
    result = optwalktime(wsteplens[1], nsteps, ctime=ctimes[1])#,negworkcost=0.2) # optimize with a cost of time
    plotvees!(result, tchange=tchange, usespline=false, color=i, speedtype=:shortwalks, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(results43, result) # add this optimization to results array
    peakspds43[i] = maximum(stepspeeds(result.steps)[2])
    durations43[i] = result.totaltime
    result = optwalktime(wsteplens[2], nsteps, ctime=ctimes[1])#,negworkcost=0.2) # optimize with a cost of time
    plotvees!(result, tchange=tchange, usespline=false, color=i, speedtype=:shortwalks, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(results, result) # add this optimization to results array
    peakspds[i] = maximum(stepspeeds(result.steps)[2])
    durations[i] = result.totaltime
    result = optwalktime(wsteplens[3], nsteps, ctime=ctimes[3])#,negworkcost=0.2) # optimize with a cost of time
    plotvees!(result, tchange=tchange, usespline=false, color=i, speedtype=:shortwalks, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(results44, result) # add this optimization to results array
    peakspds44[i] = maximum(stepspeeds(result.steps)[2])
    durations44[i] = result.totaltime
end # longer steps took longer and resulted in almost same peak speed but of course traveled farther
distances43 = [sum(result.steps.steplength) for result in results43]
distances = [sum(result.steps.steplength) for result in results]
distances44 = [sum(result.steps.steplength) for result in results44]
Plots.display(p) # instantaneous speed vs. distance profiles
# Longer steps are more costly because of collisions, but doesn't change peak speed much
# and does increase total time, and slightly increase average speed

# Variable step lengths, using the v^0.42 curve (this actually looks pretty decent)
# with variable steps, there's a steeper drop-off in speed
wstar4vs = findgait(WalkRW2lvs(α=0.35,safety=true), target=:speed=>0.3, varying=:P)
ctime = 0.05 # cost of time, to encourage hurrying
tchange = 1.75
#pv = plot()
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
resultvss = Array{MultiStepResults,1}(undef,0) # store each optimization result here
tees = zeros(length(walksteps),3)
peakspdvss = zeros(length(walksteps))
durationvss = zeros(length(walksteps))
for (i,nsteps) in enumerate(walksteps)
    result = optwalktime(wstar4vs, nsteps, ctime=ctime) # optimize with a cost of time
    plotvees!(result, tchange=tchange, usespline=false, speedtype=:shortwalks, color=i, rampuporder=1, markersize=2, linestyle=:dot) # plot instantaneous speed vs. time
    push!(resultvss, result)
    peakspdvss[i] = maximum(stepspeeds(result.steps)[2])
    durationvss[i] = result.totaltime
     # add this optimization to results array
end
Plots.display(p) # instantaneous speed vs. distance profiles
#savefig(p, "shortwalks.svg")
#savefig(p, "shortwalks.pdf")

## Short walks with speed-dependent step lengths, linearized step length (not currently using)
wstar4s = findgait(WalkRW2ls(α=0.35,safety=true), target=:speed=>0.3, varying=:P)
wstar4s = WalkRW2ls(wstar4s, vmstar=wstar4s.vm, cstep = 0.2)
ctime = 0.03 # cost of time, to encourage hurrying
tchange = 2
p = plot()
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
results = Array{MultiStepResults,1}(undef,0) # store each optimization result here
for (i,nsteps) in enumerate(walksteps)
    result = optwalktime(wstar4s, nsteps, ctime=ctime) # optimize with a cost of time
    plotvees!(result, tchange=tchange, usespline=false, color=i, speedtype=:shortwalks, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(results, result) # add this optimization to results array
end
Plots.display(p) # instantaneous speed vs. distance profiles



## Peak speed vs distance
#peakspeeds = [maximum(result.steps.vm) for result in results]     # mid-stance speeds
peakspeeds = [maximum(stepspeeds(r.steps)[2]) for r in results] # step speeds
distances = [sum(result.steps.steplength) for result in results]
p1 = plot(distances, peakspeeds, xlabel="Distance", ylabel="Peak speed", xlims=(0,Inf), ylims=(0,Inf))
p2 = plot(walksteps, peakspeeds, xlabel="# of steps", ylabel="Peak speed")
peakspeeds43 = [maximum(stepspeeds(r.steps)[2]) for r in results43] # step speeds
distances43 = [sum(result.steps.steplength) for result in results43]
peakspeeds44 = [maximum(stepspeeds(r.steps)[2]) for r in results44] # step speeds
distances44 = [sum(result.steps.steplength) for result in results44]
peakspeedvss = [maximum(stepspeeds(r.steps)[2]) for r in resultvss] # step speeds
distancevss = [sum(result.steps.steplength) for result in resultvss]
plot!(p1, distances43, peakspeeds43)
plot!(p1, distances44, peakspeeds44)
plot!(p1, distancevss, peakspeedvss)
plot!(p2, walksteps, peakspeeds43, xlabel="# of steps", ylabel="Peak speed")
plot!(p2, walksteps, peakspeeds44, xlabel="# of steps", ylabel="Peak speed")
plot!(p2, walksteps, peakspeedvss, xlabel="# of steps", ylabel="Peak speed")
plot(p1, p2, layout = (1,2), legend=false)
#savefig("peakshortwalks.svg")
#savefig("peakshortwalks.pdf")

## Short walks: Time to walk a distance
# A fairly linear increase in time to walk a distance, but with a slight curved toe-in
timetowalk = [result.totaltime+tchange for result in results]
plot(distances, timetowalk, xlims=(0,Inf), ylims=(0,Inf),
    xguide="Distance", yguide="Time", title="Time to walk a distance", label=nothing)
# savefig("durationdistance.pdf")


## Short walks: Varying ctimes to demonstrate self-similarity THE BIG PLOT
wstar4 = findgait(WalkRW2l(α=0.35), target=:speed=>0.3, varying=:P)
ctimes = (0.006, 0.015, 0.0276, 0.0384, 0.0492, 0.06)
tchange = 1.75
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
peaks = zeros(length(walksteps),length(ctimes))
durations = similar(peaks)
results = Array{MultiStepResults,2}(undef,(length(walksteps),length(ctimes))) # store each optimization result here
for (j,ctime) in enumerate(ctimes)
    for (i,nsteps) in enumerate(walksteps)
        result = optwalktime(wstar4, nsteps, ctime=ctime) # optimize with a cost of time
        #peaks[i,j] = maximum(result.steps.vm)
        peaks[i,j] = maximum(stepspeeds(result.steps)[2])
        durations[i,j] = result.totaltime
        results[i,j] = result
    end
end
# after the fact, let's plot them all on top of each other (scaling time and velocity)
# using the ctime=0.02 result as the basis
tbase = durations[end,2]
vbase = peaks[end,2]
#layout = @layout[ a{0.85w} grid(6,1)]
pleft = plot(; )#@layout [grid(1,3); a{0.86h}])
pright = plot(; layout=grid(6,1))
ptop = plot(; layout=grid(1,3))
#p = plot(;layout)
#Plots.display(p) 
for (j, ctime) in enumerate(ctimes)
    for (i,nsteps) in enumerate(walksteps)
        result = results[i,j]
        plotvees!(pright,result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=0, subplot=j,
            xticks = [20,40], yticks=[0.2,0.4,0.6],xguide="",yguide="",tickfontsize=4,
            xlims=(0,maximum(durations)+3tchange), ylims=(0,maximum(peaks)), linewidth=0.5) # subplot instantaneous speed vs. time
        plotvees!(pleft,result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=2, tscale = tbase/(durations[end,j]), 
            vscale = vbase/peaks[end,j],subplot=1) # main scaled speed vs time
    end
end
for (i,result) in enumerate(resultvss) # add in the variable step length results computed above in resultvss (dash-dot lines)
    plotvees!(pleft,result, tchange=tchange, usespline=false, speedtype=:shortwalks, color=i, markersize=2, linestyle=:dashdot,subplot=1,
    tscale = tbase/(durationvss[end]),vscale = vbase/peakspdvss[end]) # plot instantaneous speed vs. time
    plotvees!(ptop,result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=0, 
        xticks = [20,40], yticks=[0.2,0.4,0.6],subplot=2,xguide="",yguide="",tickfontsize=4,
        xlims=(0,maximum(durations)+3tchange), ylims=(0,maximum(peaks)), linewidth=0.5)
end
for (i,result) in enumerate(results43) # add in shorter steps (dashed lines)
    plotvees!(pleft,result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=2, tscale = tbase/(durations43[end]), 
        vscale = vbase/peakspds43[end],subplot=1, linestyle=:dash)
    plotvees!(ptop,result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=0, 
        xticks = [20,40], yticks=[0.2,0.4,0.6],subplot=1,xguide="",yguide="",tickfontsize=4,
        xlims=(0,maximum(durations)+3tchange), ylims=(0,maximum(peaks)), linewidth=0.5)
end
for (i,result) in enumerate(results44) # add in longer steps (dotted lines)
    plotvees!(pleft,result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=2, tscale = tbase/(durations44[end]), 
        vscale = vbase/peakspds44[end],subplot=1, linestyle=:dot)
    plotvees!(ptop,result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=0,  
        xticks = [20,40], yticks=[0.2,0.4,0.6],subplot=3, xguide="",yguide="",tickfontsize=4,
        xlims=(0,maximum(durations)+3tchange), ylims=(0,maximum(peaks)), linewidth=0.5)
end
plot(ptop, pleft, pright, layout= @layout [ [a{0.15h}; b] c{0.15w}])
#plot(pleft, pright, layout=grid(1,2,widths=(0.85,0.15)))
Plots.display(p)
println("Durations of a factor of ", (durations[end,1]+2tchange)/(durations[end,end]+2tchange))
println("Peak speeds over a range of ", peaks[end,end]/peaks[end,1])
println("  about ", peaks[end,1]*sqrt(9.81)," to ", peaks[end,end]*sqrt(9.81), "m/s")
peaks[end,:]*sqrt(9.81)
#savefig("shortwalks.pdf")
#savefig("shortwalks.svg")
#0.9904474300428026
#1.3233122108084778
#1.6132937663732232
#1.797862141973894
#1.950737949278336
#2.0827899599175965
# c_t = 0.0126 yields 1.25 m/s, 0.018 yields 1.4 m/s which corresponds to c_t*Mg^1.5L^0.5
#  27 W and 38.7 W; typical work is P = 0.0122*MgL = 8 J push-off, which is less than human 20 J.
# If you make it relative to human work per step, you get close to 100W.
#savefig("selfsimilarity.pdf")
#savefig("selfsimilarity.svg")
# GR no fonts, doesn't do eps
# plotlyJS did export fonts, not necessarily the right one
# pyplot doesn't preserve fonts, but does export eps
#Durations of a factor of 1.9312552294263987
#Peak speeds over a range of 2.1028778476688945
#  about 0.9904474300428026 to 2.0827899599175965m/s

## Short walks: Vary c and # steps more tightly to get a plot of speed vs C
wstar4 = findgait(WalkRW2l(α=0.35), target=:speed=>0.3, varying=:P)
ctimes = (0.006, 0.015, 0.0276, 0.0384, 0.0492, 0.06)
morectimes = range(ctimes[begin], ctimes[end], length=16)
walksteps = [1,2,3,4,5,6, 7, 10, 15, 20] # take walks of this # of steps

results = [optwalktime(wstar4, nsteps, ctime=ctime) for nsteps in walksteps, ctime in morectimes]
peaks = [maximum(stepspeeds(r.steps)[2]) for r in results]
durations = [r.totaltime for r in results]
sec = sqrt(1/9.81); mps = sqrt(9.81)
p1=plot(morectimes, mps.*peaks',legend=false,xlabel="c_t",ylabel="peak speed",linewidth=0.2)
plot!(morectimes, mps.*(peaks .* middle(peaks[:,end])./ maximum(peaks,dims=2))',linewidth=1,legend=false, xlabel="c_t",ylabel="peak speed")
plot!([0.04,0.04], [1,2].*0.2*mps,lims=(0,Inf)) # 0.2*sqrt(gL)
# right now the cost of time is integral (Wdot+c)*dt
# or c_t = [J/s] in real units, or [(J/MgL)/(time/sqrt(L/g)] 
# so 0.06*MgL/sqrt(L/g) = Mg^1.5*L^0.5 = 0.03*2147 = 64 J/s seems high

# Plot peak speed vs distance, both regular and normalized to each other
distances = [sum(result.steps.steplength) for result in results]
p1=plot(distances, peaks, xlabel="Distance", ylabel="Peak speed", xlims=(0,maximum(distances)), ylims=(0,Inf),legend=false)
plot!(p1, distances43, peakspeeds43, linestyle=:dash)
plot!(p1, distances44, peakspeeds44, linestyle=:dot)
plot!(p1, distancevss, peakspeedvss, linestyle=:dashdot)
p2=plot(distances, (peaks .* middle(peaks[end,:]) ./ maximum(peaks,dims=1)), xlabel="Distance", ylabel="Peak speed",legend=false,xlims=(0,maximum(distances)) )
plot!(p2, [4 0; 4 0.5],[0.1 1.25 ;1.1 1.25]./mps) # 1 m/s
plot!(p2, distances43, (peakspeedvss .* middle(peaks[end,:]) ./ maximum(peakspeedvss,dims=1)), linestyle=:dash)
plot!(p2, distances43, (peakspeeds43 .* middle(peaks[end,:]) ./ maximum(peakspeeds43,dims=1)), linestyle=:dot)
plot!(p2, distances43, (peakspeeds44 .* middle(peaks[end,:]) ./ maximum(peakspeeds44,dims=1)), linestyle=:dashdot)
plot(p1,p2,layout=(1,2),link=:y, legend=false)
#savefig("peakspeedvdist.pdf")
#savefig("peakspeedvdist.svg")

# peak speed increases with distance, and of course c shifts it upward
# but basically the idea is that c and distance are the only two parameters, and everything collapses
# onto a single line
# And do duration vs distance
p1=plot(distances, durations, xlabel="Distance", ylabel="Duration",legend=false,linewidth=0.2)
plot!(distances, durations .* middle(durations[end,:])./ maximum(durations,dims=1), linewidth=0.5,xlabel="Distance", ylabel="Duration",legend=false)
plot!([5,5],[1,2]./sec,xlims=(0,maximum(distances))) # 1 sec
plot!(distances43, durations43 .* middle(durations[end,:])./ maximum(durations43,dims=1), linewidth=0.5,xlabel="Distance", ylabel="Duration",legend=false,linestyle=:dash)
plot!(distances44, durations44 .* middle(durations[end,:])./ maximum(durations44,dims=1), linewidth=0.5,xlabel="Distance", ylabel="Duration",legend=false,linestyle=:dash)
plot!(distancevss, durationvss .* middle(durations[end,:])./ maximum(durationvss,dims=1), linewidth=0.5,xlabel="Distance", ylabel="Duration",legend=false,linestyle=:dash)
#savefig("durationdistance.pdf")
#savefig("durationdistance.svg")



## Short walks: Walk a certain number of steps, with minimum energy-time vs steady minCOT
# where it's possible to stay almost exactly at minCOT, except for start-up 
# this is like trapezoidal, but using an objective to make every step have same speed
wstar4s = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.5, varying=:P)
wstar4n = findgait(WalkRW2l(α=0.35, safety=true), target=:speed=>0.4, varying=:P)
nsteps = 10
ctime = 0.0195
tchange = 1.75
nominalmsr=optwalktime(wstar4s, nsteps, ctime = ctime, boundarywork=true) # to compare with our usual solution
minvarmsr=optwalkvar(wstar4n, nsteps, boundarywork=true)
p = plot(layout=(2,1))
plotvees!(p[1],nominalmsr, tchange=tchange, rampuporder=1, usespline = false, speedtype=:shortwalks)
plotvees!(p[1],minvarmsr, tchange=tchange, rampuporder=1, usespline = false, speedtype=:shortwalks)
plot!(p[2],[0:nsteps+1], [1/2*nominalmsr.vm0^2; nominalmsr.steps.Pwork; NaN],markershape=:circle)
plot!(p[2],[0:nsteps+1], [NaN; nominalmsr.steps.Cwork; -1/2*nominalmsr.steps[end].vm0^2], markershape=:circle)
plot!(p[2],[0:nsteps+1], [1/2*minvarmsr.vm0^2; minvarmsr.steps.Pwork; NaN],markershape=:circle,xticks=0:nsteps+1,yticks=(-0.1,0.1))
plot!(p[2],[0:nsteps+1], [NaN; minvarmsr.steps.Cwork; -1/2*minvarmsr.steps[end].vm0^2], markershape=:circle)
plot!(p[2],xlabel="step", ylabel="push-off work", legend=false)
Plots.display(p)
#savefig("twohypotheses.pdf")
println("energy-time work = ", 1/2*minvarmsr.vm0^2 + sum(minvarmsr.steps.Pwork))
println("max steady = ", 1/2*nominalmsr.vm0^2 + sum(nominalmsr.steps.Pwork))
println("ratio = ",  (1/2*minvarmsr.vm0^2 + sum(minvarmsr.steps.Pwork))/(1/2*nominalmsr.vm0^2 + sum(nominalmsr.steps.Pwork)) )
println("ratio = ",  (sum(minvarmsr.steps.Pwork))/(sum(nominalmsr.steps.Pwork)) )
# or about 12.7% more costly to do steady gait

# use the following to see all the plots
#multistepplot(nominalmsr)
#multistepplot!(minvarmsr)

## Short walks: Trapezoid comparison for different numbers of steps
wstar4 = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.3, varying=:P)
ctime = 0.015 # cost of time, to encourage hurrying
tchange = 1.75
layout = @layout[ grid(2,1) b{0.40w}]
p = plot(;layout,legend=false)
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
results = Array{MultiStepResults,1}(undef,0) # store each optimization result here
for (i,nsteps) in enumerate(walksteps)
    result = optwalktime(wstar4, nsteps, ctime=ctime)#,negworkcost=0.2) # optimize with a cost of time
    plotvees!(p[2],result, tchange=tchange, usespline=false, color=i, speedtype=:shortwalks, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(results, result) # add this optimization to results array
end
tresults = Array{MultiStepResults,1}(undef,0) # store each optimization result here
for (i,nsteps) in enumerate(walksteps)
    result = optwalkvar(wstar4, nsteps, boundarywork=true) # optimize with a cost of variance
    plotvees!(p[1],result, tchange=tchange, usespline=false, color=i, speedtype=:shortwalks, rampuporder=1, markersize=2) # plot instantaneous speed vs. time
    push!(tresults, result) # add this optimization to results array
end


Plots.display(p) # instantaneous speed vs. distance profiles


# TODO: Make this comparable in duration


## Short walks: A ctimes parameter study with both fixed and varying step lengths in subplots
# Not currently using this, because we can already show fixed & varying in left plot
wstar4 = findgait(WalkRW2l(α=0.35), target=:speed=>0.3, varying=:P)
wstar4vs = findgait(WalkRW2lvs(α=0.35, safety=true), target=:speed=>0.3, varying=:P)
ctimes = range(0.006, 0.06, length=6)*2
tchange = 2
layout = @layout[ a{0.85w} grid(6,1)]
p = plot(;layout)
peaks = zeros(length(walksteps),length(ctimes),2)
durations = similar(peaks)
walksteps = [2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
results = Array{MultiStepResults,3}(undef,(length(walksteps),length(ctimes),2)) # store each optimization result here
for (k,w) in enumerate((wstar4, wstar4vs))
#k = 1; w = wstar4vs
for (j,ctime) in enumerate(ctimes)
    for (i,nsteps) in enumerate(walksteps)
        result = optwalktime(w, nsteps, ctime=ctime) # optimize with a cost of time
        #peaks[i,j] = maximum(result.steps.vm)
        peaks[i,j,k] = maximum(stepspeeds(result.steps)[2])
        durations[i,j,k] = result.totaltime
        results[i,j,k] = result
    end
end
end
# after the fact, let's plot them all on top of each other
# using the ctime=0.02 result as the basis
tbase = durations[end,2,1]
vbase = peaks[end,2,1]
for (k,w) in enumerate((wstar4, wstar4vs))
#k = 1; w = wstar4vs
for (j, ctime) in enumerate(ctimes)
    for (i,nsteps) in enumerate(walksteps)
        result = results[i,j,k]

        plotvees!(result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=2, subplot=j+1,
            xticks = [20,40], yticks=[0.2,0.4,0.6],xguide="",yguide="",tickfontsize=4,
            xlims=(0,maximum(durations)+2tchange), ylims=(0,maximum(peaks))) # plot instantaneous speed vs. time
        plotvees!(result, tchange=tchange, color=i, usespline=:false, speedtype=:shortwalks,markersize=2, tscale = tbase/(durations[end,j,1]), 
            vscale = vbase/peaks[end,j,1],subplot=1)
    end
end
end
Plots.display(p) 
#println("Durations of a factor of ", (durations[end,1]+2tchange)/(durations[end,end]+2tchange))
#println("Peak speeds over a range of ", peaks[end,end]/peaks[end,1])
#println("  about ", peaks[end,1]*sqrt(9.81)," to ", peaks[end,end]*sqrt(9.81), "m/s")




## Short walks: Up and down slopes
# Compare walking uphill, downhill, and level, for a fixed number of steps, and including
# optimized time. The results show uphill walking is skewed to fast speed-up at beginning,
# slow coasting toward end. Downhill is skewed for slow speed-up aided by gravity, and
# abrupt end.
# `optwalktime` is an optimization for a given number of steps, with a `ctime` cost.

# Make plots comparing up, level, down for various numbers of steps
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.4, varying=:P)
myslope = 0.08
ctime = 0.02
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

## Short walks: Varying fixed step lengths
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

## Short walks: Real trapezoid, using tchange to get up to speed
# this gives same results as the optwalkvar optimization, except here we specify all the
# push-offs and just try to get the velocity to be square
wstar4 = findgait(WalkRW2l(α=0.35), target=:speed=>0.45, varying=:P)
tchange = 1.75
nsteps = 5
trapezoidresult=multistep(wstar4, Ps=fill(wstar4.P, nsteps), boundaryvels=(0.,0.),extracost=1/2*wstar4.vm^2)
plotvees(trapezoidresult, tchange=tchange, speedtype=:shortwalks, rampuporder=1, usespline=false)

    
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
#savefig("trapezoidshortwalks.svg")

# If you just want to do square wave in speed, it costs a lot of initial push-off
# so let's compare with walking the same number of steps and same amount of time
optresults=optwalk(w, Ncruise, boundaryvels=(0,0),totaltime=trapezoidresults.totaltime  )
println("trapezoid cost = ", trapezoidresults.totalcost, "   optimal cost = ", optresults.totalcost)
# It's definitely more expensive to use the square wave
multistepplot!(optresults,plotwork=true,label="optimal")
#savefig("trapezoidcomparisonshortwalks.svg")

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
nsteps = 21
tfstar = onestep(wstar4s).tf
δs = zeros(nsteps); δs[floor(Int,(nsteps+1)/2)] = 0.05
nominalmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δs)
plotvees(nominalmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)
nhorizon = 13 # so far we shouldn't reduce the horizon much
# because we're assuming you want to match a twin brother at horizon
# but we could also try to match our original plan

elapsedtime = 0
modeltime = 0
currentvm0 = nominalmsr.steps[1].vm0
mpcsteps = Vector{StepResults}(undef,nsteps)
upcomingδs = zeros(nhorizon)
for i in 1:nsteps-1 # receding horzion; don't optimize the last step
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
    #plot!(cumsum([elapsedtime;nextmsr.steps[1:end-1].tf]),
    #    [nextstep.vm0;nextmsr.steps.vm[1:end-1]],show=true)
    plot!(cumsum([elapsedtime;nextmsr.steps[1].tf]),
        [nextstep.vm0;nextmsr.steps.vm[1]],show=true)
    remainingtime = remainingtime - nextstep.tf
    # set up for the next one
    currentvm0 = nextstep.vm
    elapsedtime = elapsedtime + nextstep.tf
    modeltime = modeltime + tfstar
end





## Walk over a single bump with fixed and with varying step lengths
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
nsteps = 15
δs = zeros(nsteps); δs[Int((nsteps+1)/2)] = 0.05 # one bump
nominalmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δs)

# WalkRW2ls has varying step lengths according to preferred human
wstar4ls = findgait(WalkRW2ls(α=0.4,safety=true), target=:speed=>0.45, varying=:P, cstep=0.35, vmstar=wstar4s.vm)
wstar4lvs = findgait(WalkRW2lvs(α=0.4,safety=true), target=:speed=>0.45, varying=:P, c=1., vmstar=wstar4s.vm)
varyingmsr = optwalk(wstar4ls, nsteps, boundarywork=false,δs=δs)
varyingmsrv = optwalk(wstar4lvs, nsteps, boundarywork=false,δs=δs)
plotvees(nominalmsr,boundaryvels=nominalmsr.boundaryvels, speedtype=:shortwalks)
plotvees!(varyingmsr,boundaryvels=(nothing,nothing), speedtype=:shortwalks)
plotvees!(varyingmsrv,boundaryvels=(nothing,nothing), speedtype=:shortwalks)

plot(cumsum(nominalmsr.steps.tf), nominalmsr.steps.vm,label="normal")
plot!(cumsum(varyingmsr.steps.tf), varyingmsr.steps.vm, label="varying", )

# step timings, per step, regular and varying step lengths
plot(cumsum(nominalmsr.steps.tf),nominalmsr.steps.tf)
plot!(cumsum(varyingmsr.steps.tf),varyingmsr.steps.tf)


## Walk over a single bump to minimize variance
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
nsteps = 15
δs = zeros(nsteps); δs[Int((nsteps+1)/2)] = 0.05
varmsr=optwalkvar(wstar4s, nsteps, boundarywork=false, δs=δs)
plotvees(nominalmsr,boundaryvels=nominalmsr.boundaryvels, speedtype=:midstance)
plotvees!(varmsr, boundaryvels=nominalmsr.boundaryvels, speedtype=:midstance)
multistepplot(varmsr, label=["tight speed" "" ""])
multistepplot!(nominalmsr,label=["min work" "" ""])
savefig("tightspeed1.pdf")
varmsr.totalcost, varmsr.totaltime
nominalmsr.totalcost, nominalmsr.totaltime
onestep(wstar4s).tf*nsteps
# plot cumulative time gain
levelsteptimes = fill(onestep(wstar4s).tf, nsteps)
bumpsteptimes = nominalmsr.steps.tf
varmsrsteptimes = varmsr.steps.tf

plot(-cumsum(bumpsteptimes .-levelsteptimes),label="min work")
plot!(-cumsum(varmsrsteptimes .- levelsteptimes),label="tight regulation", xlabel="steps", ylabel="cum time gain")
savefig("tightspeed2.pdf")

multistep
controlcost = nsteps*(0.5*wstar4s.P^2)

## Walk over a simple bump with no compensation
nocompmsr = multistep(wstar4s, Ps=fill(wstar4s.P,nsteps),δangles=δs,boundaryvels=(wstar4s.vm,wstar4s.vm))
multistepplot!(nocompmsr)
nocompmsr.totalcost

## Walk over a single bump with a reactive compensation
nbump = Int((nsteps+1)/2)
reactmsr1 = multistep(wstar4s, Ps=fill(wstar4s.P,nbump),δangles=δs[1:nbump],boundaryvels=(wstar4s.vm,wstar4s.vm))
reactmsr2 = optwalk(wstar4s, nsteps-nbump, totaltime = nominalmsr.totaltime - reactmsr1.totaltime - 2,boundaryvels=(reactmsr1.steps[end].vm,wstar4s.vm), boundarywork=(false,false))
multistepplot(reactmsr2)



## Triangle walk, based on min var walk
# Short walks: Walk a certain number of steps, with minimum energy-time vs steady minCOT
# where it's possible to stay almost exactly at minCOT, except for start-up 
# this is like trapezoidal, but using an objective to make every step have same speed
wstar4s = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.5, varying=:P)
wstar4n = findgait(WalkRW2l(α=0.35, safety=true), target=:speed=>0.4, varying=:P)
nsteps = 10
ctime = 0.0195
tchange = 1.75
nominalmsr=optwalktime(wstar4s, nsteps, ctime = ctime, boundarywork=true) # to compare with our usual solution
minvarmsr=optwalkvar(wstar4n, nsteps, boundarywork=true)
A = 1.9*wstar4s.vm/(nsteps*onestep(wstar4s).tf)
v0 = 0.11#0.8*A*tchange#0.12
mintrimsr=optwalktriangle(wstar4n, nsteps, A = A, boundarywork=false,boundaryvels=(v0,v0))
p = plot(layout=(1,2))
plotvees!(p[1],nominalmsr, tchange=tchange, rampuporder=1, usespline = false, markershape=:circle,speedtype=:shortwalks)
plotvees!(p[1],minvarmsr, tchange=tchange, rampuporder=1, usespline = false,markershape=:circle, speedtype=:shortwalks)
plotvees!(p[1],mintrimsr, tchange=tchange, rampuporder=1, usespline = false,markershape=:circle, speedtype=:shortwalks)
plot!(p[2],[0:nsteps+1], [1/2*nominalmsr.vm0^2; nominalmsr.steps.Pwork; NaN],markershape=:circle)
#plot!(p[2],[0:nsteps+1], [NaN; nominalmsr.steps.Cwork; -1/2*nominalmsr.steps[end].vm0^2], markershape=:circle)
plot!(p[2],[0:nsteps+1], [1/2*minvarmsr.vm0^2; minvarmsr.steps.Pwork; NaN],markershape=:circle,xticks=0:nsteps+1)
#plot!(p[2],[0:nsteps+1], [NaN; minvarmsr.steps.Cwork; -1/2*minvarmsr.steps[end].vm0^2], markershape=:circle)
plot!(p[2],[0:nsteps+1], [1/2*mintrimsr.vm0^2; mintrimsr.steps.Pwork; NaN],markershape=:circle,xticks=0:nsteps+1)
#plot!(p[2],[0:nsteps+1], [NaN; mintrimsr.steps.Cwork; -1/2*minvarmsr.steps[end].vm0^2], markershape=:circle)
plot!(p[2],xlabel="step", ylabel="push-off work", legend=false)
Plots.display(p)

#savefig("threehypotheses.pdf")
println("energy-time work = ", 1/2*nominalmsr.vm0^2 + sum(nominalmsr.steps.Pwork))
println("max steady = ", 1/2*minvarmsr.vm0^2 + sum(minvarmsr.steps.Pwork))
println("triangle   = ", 1/2*mintrimsr.vm0^2 + sum(mintrimsr.steps.Pwork))
println("ratio = ",  (1/2*minvarmsr.vm0^2 + sum(minvarmsr.steps.Pwork))/(1/2*nominalmsr.vm0^2 + sum(nominalmsr.steps.Pwork)) )
println("ratio = ",  (1/2*mintrimsr.vm0^2 + sum(mintrimsr.steps.Pwork))/(1/2*nominalmsr.vm0^2 + sum(nominalmsr.steps.Pwork)) )
#println("ratio = ",  (sum(minvarmsr.steps.Pwork))/(sum(nominalmsr.steps.Pwork)) )

threecosts = [1/2*nominalmsr.vm0^2 + sum(nominalmsr.steps.Pwork), 1/2*minvarmsr.vm0^2 + sum(minvarmsr.steps.Pwork), 1/2*mintrimsr.vm0^2 + sum(mintrimsr.steps.Pwork)]
bar(threecosts,xticks=((1,2,3),("Energy-Time", "Steady min-COT", "Steady accel")))
savefig("threebars.pdf")