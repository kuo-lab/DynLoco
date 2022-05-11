
using DynLoco, JuMP, Ipopt
using Plots, Statistics
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

# catch someone walking at a certain speed and distance
function optwalktag(w::W, nsteps=50; boundaryvels::Union{Tuple,Nothing} = (0.,0.), safety=true,
    ctime = 0.05, tchange = 3., boundarywork = true, δs = zeros(nsteps), startv = w.vm, negworkcost = 0., targetspeed = 0.4, targetdist = 1., cerror=0., walkparms...) where W <: Walk
    #println("walkparms = ", walkparms)
    w = W(w; walkparms...)
    optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>1))
    @variable(optsteps, P[1:nsteps]>=0, start=w.P) # JuMP variables P
    # constraints: starting guess for velocities
    
    if length(startv) == 1 # startv can be just a scalar, a 1-element vector, or n elements
        @variable(optsteps, v[1:nsteps+1]>=0, start=startv)
    else # starting guess for v is an array
        @variable(optsteps, v[i=1:nsteps+1]>=0, start=startv[i])
    end
    register(optsteps, :onestepv, 3, (v,P,δ)->onestep(w,P=P,vm=v,δangle=δ,safety=safety).vm, autodiff=true) # input P, output vm
    register(optsteps, :onestept, 3, (v,P,δ)->onestep(w,P=P,vm=v,δangle=δ,safety=safety).tf, autodiff=true)
    register(optsteps, :onestepc, 3, (v,P,δ)->onestep(w,P=P,vm=v,δangle=δ,safety=safety).C, autodiff=true)
    register(optsteps, :onestepcps, 3, (v,P,δ)->onestep(w,P=P,vm=v,δangle=δ,safety=safety).costperstep, autodiff=true)
    register(optsteps, :softrelu, 1, x->0.5*log(1+2*x), autodiff=true)
    register(optsteps, :sigmoid, 1, x->1/(1+exp(-8*x)), autodiff=true)
    @NLexpression(optsteps, steptime[i=1:nsteps], onestept(v[i],P[i],δs[i]))
    @NLexpression(optsteps, targetposn[i=1:nsteps], targetdist + targetspeed*sum(steptime[j] for j in 1:i))
    steplength = 2*sin(w.α)
    @NLexpression(optsteps, posn[i=1:nsteps], i*steplength)
    @NLexpression(optsteps, totaltime, sum(onestept(v[i],P[i],δs[i]) for i = 1:nsteps))
    @NLexpression(optsteps, targetahead[i=1:nsteps], sigmoid(targetposn[i]-posn[i]))
    for i = 1:nsteps # collocation points
        @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i],δs[i]))
    end
    if boundarywork
        # mean-square target error, plus work, plus time
        @NLobjective(optsteps, Min, 1/2*(sum(P[i]^2 for i=1:nsteps) + negworkcost*sum(onestepc(v[i],P[i],δs[i])^2 for i=1:nsteps) + v[1]^2 - boundaryvels[1]^2 ) +
            ctime*sum(targetahead[i]*steptime[i] for i=1:nsteps) + cerror*sum((targetposn[i]-posn[i])^2 for i=1:nsteps))
#        @NLobjective(optsteps, Min, 1/2*(sum(P[i]^2 for i=1:nsteps) + negworkcost*sum(onestepc(v[i],P[i],δs[i])^2 for i=1:nsteps) + v[1]^2 - boundaryvels[1]^2 ) +
#            ctime*sum(targetahead[i]*steptime[i] for i=1:nsteps))
    else
        @NLobjective(optsteps, Min, 1/2*(sum(P[i]^2 for i=1:nsteps) + negworkcost*sum(onestepc(v[i],P[i],δs[i])^2 for i=1:nsteps)) +
            ctime*totaltime)
    end
    optimize!(optsteps)
    if termination_status(optsteps) == MOI.LOCALLY_SOLVED || termination_status(optsteps) == MOI.OPTIMAL
        optimal_solution = (vms=value.(v), Ps=value.(P))
    else
        error("The model was not solved correctly.")
        println(termination_status(optsteps))
    end
    @show value.(targetahead)
    result = multistep(W(w,vm=value(v[1]),safety=safety), value.(P), δs, vm0=value(v[1]),
        boundaryvels=boundaryvels, extracost = ctime*value(totaltime) +
        (boundarywork ? 1/2*(value(v[1])^2-boundaryvels[1]^2) : 0))
    return result
end


wstar4 = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.3, varying=:P)
ctime = 0.015 # cost of time, to encourage hurrying
tchange = 1.75
p = plot(layout=(2,1))
targetspeed = 0.3
targetdist = 0.5
result = optwalktag(wstar4, 50, ctime=ctime, targetspeed=targetspeed, targetdist=targetdist)
targettime = range(0, result.totaltime, length=50)
targetpos(t) = targetdist .+ targetspeed.*t
chasertime = cumsum([0;result.steps.tf])
chaserposn = cumsum([0; result.steps.steplength])
catcher = findfirst(chaserposn .>= targetpos.(chasertime))
plot!(p[1], chasertime[1:catcher], result.steps.speed[1:catcher],ylims=(0,Inf))
plot!(p[2], chasertime[1:catcher], targetpos(chasertime[1:catcher]))
plot!(p[2], chasertime[1:catcher], chaserposn[1:catcher])


## Methods: 
# minimize mean-square error to target, plus work and time

# minimize time to reach target (cost of time) + work to get there
# (but need to figure out a cost afterwards)


# minimum gross cost - temporally discounted reward
#                    + cost of task time  (which yields more after-task time)
# some tasks have a fixed reward, and temporal discounting favors getting it sooner
# but we find that a linear cost of time actually does pretty well

# Mean square error
wstar4 = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.3, varying=:P)
ctime = 0.015 # cost of time, to encourage hurrying
cerror = 0.015
targetspeed = 0.3
targetdist = 4.5
result = optwalktag(wstar4, 20, ctime=0, cerror=cerror,targetspeed=targetspeed, targetdist=targetdist)
p = plot(layout=(2,1))
targettime = range(0, result.totaltime, length=50)
targetpos(t) = targetdist .+ targetspeed.*t
chasertime = cumsum([0;result.steps.tf])
chaserposn = cumsum([0; result.steps.steplength])
catcher = length(chaserposn)-1 #findfirst(chaserposn .>= targetpos.(chasertime))
plot!(p[1], chasertime[1:catcher], result.steps.speed[1:catcher],ylims=(0,Inf))
plot!(p[2], chasertime[1:catcher], targetpos(chasertime[1:catcher]))
plot!(p[2], chasertime[1:catcher], chaserposn[1:catcher])
# MSE you want to hit a peak speed quickly and then slow down toward the target.
# with no ctime

