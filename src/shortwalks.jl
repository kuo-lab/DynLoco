using DynLoco
using Plots

## Short walks of different distances
wstar4 = findgait(Walk(α=0.4), target=:speed=>0.4, varying=:P)
p = plot()
results = Array{MultiStepResults,1}(undef,0)
for (i,nsteps) in enumerate([2, 3, 4, 5, 7, 10, 15, 20])
    result = optwalktime(wstar4, nsteps, ctime=0.05)
    plotvees!(result, tchange=3, color=i)
    push!(results, result)
end
Plots.display(p)


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
# Need to fix the total cost calculation
wstar = findgait(wrw, target=:speed=>0.4, varying=:P)
myout1 = optwalkslope(wstar, 6, boundaryvels = (0., 0.), boundarywork= true, symmetric = true)
multistepplot(myout1; plotwork=true)
myout2 = optwalk(wstar, 6, boundaryvels=(0.,0.), totaltime=myout1.totaltime,
    δ=myout1.δangles) # Yes we get the same result
multistepplot!(myout2, plotwork=true)
# level ground
myout3 = optwalk(wstar, 6, boundaryvels=(0.,0.), totaltime=myout1.totaltime, δ=zeros(6)) # level
multistepplot!(myout3, plotwork=true)
# upward bulge
myout4 = optwalk(wstar, 6, boundaryvels=(0.,0.), totaltime=myout1.totaltime, δ=-myout1.δangles) # level
multistepplot!(myout4, plotwork=true)
Plots.display(Plots.CURRENT_PLOT.nullableplot)
println("work = ", myout1.totalcost, " ", myout2.totalcost, " ", myout3.totalcost, " ", myout4.totalcost)

# Try a few different speeds
myout1 = optwalkslope(wstar, 6, boundaryvels = (0., 0.), boundarywork= true, symmetric = true)
bestramp = myout1.δangles
myout2 = optwalk(wstar, 6, boundaryvels=(0.,0.), totaltime=myout1.totaltime,
    δ=bestramp) # Yes we get the same result
multistepplot(myout2, plotwork=true)
myout5 = optwalk(wstar, 6, boundaryvels=(0.,0.), totaltime=1.1*myout1.totaltime,
    δ=bestramp) # a bit slower
multistepplot!(myout5, plotwork=true)
myout6 = optwalk(wstar, 6, boundaryvels=(0.,0.), totaltime=1.2*myout1.totaltime,
    δ=bestramp) # a bit slower
multistepplot!(myout6, plotwork=true)
println("work = ", myout1.totalcost, " ", myout2.totalcost, " ", myout5.totalcost, " ", myout6.totalcost)

# and compare against level ground at three speeds
myout7 = optwalk(wstar, 6, boundaryvels=(0.,0.), totaltime=1.1*myout1.totaltime)
myout8 = optwalk(wstar, 6, boundaryvels=(0.,0.), totaltime=1.2*myout1.totaltime)
multistepplot!(myout7, plotwork=true)
multistepplot!(myout8, plotwork=true)
println("work = ", myout1.totalcost, " ", myout5.totalcost, " ", myout6.totalcost) # ramp slower

println("work = ", myout5.totalcost, " ", myout7.totalcost, " ", myout8.totalcost) # level slower
