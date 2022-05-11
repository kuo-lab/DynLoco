# 
using JuMP, Ipopt, DynLoco, Plots
# optimize only the horizontal or vertical component of work or force
function optwalkcomp(w::W, numsteps=5; boundaryvels::Union{Tuple,Nothing} = nothing,
    boundarywork::Union{Tuple{Bool,Bool},Bool} = (true,true), totaltime=numsteps*onestep(w).tf,
    δs = zeros(numsteps), horvert = 1, impwork = 2) where W <: Walk # default to taking the time of regular steady walking

    #cumdelta = cumsum(δs)

    # horvert = 1 means use horizontal, 2 = vertical, 3 = both
    if horvert == 1
        comp = sin.(w.α .- δs)
    elseif horvert == 2
        comp = cos.(w.α .- δs)
    elseif horvert == 3 # use both horizontal and vertical, i.e. push-off (set impwork = 2)
        comp = ones(length(δs))
    else
        error("need horvert = 1 or 2")
    end

    # impwork =  1 for impulse, = 2 for work
    if !(impwork == 1 || impwork == 2)
        error("need impwork = 1 or 2")
    end
    
    optsteps = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0))
    @variable(optsteps, P[1:numsteps]>=0, start=w.P) # JuMP variables P
    @variable(optsteps, v[1:numsteps+1]>=0, start=w.vm) # mid-stance speeds

    if boundaryvels === nothing || isempty(boundaryvels)
        boundaryvels = (w.vm, w.vm) # default to given gait if nothing specified
    end

    if typeof(boundarywork) <: Bool
        boundarywork = (boundarywork, boundarywork)
    end

    if !boundarywork[1] # no hip work at beginning or end; apply boundary velocity constraints
        @constraint(optsteps, v[1] == boundaryvels[1])
    end
    if !boundarywork[2]
        @constraint(optsteps, v[numsteps+1] == boundaryvels[2])
    end

    # Constraints
    # produce separate functions for speeds and step times
    register(optsteps, :onestepv, 3, # velocity after a step
        (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).vm, autodiff=true) # output vm
    register(optsteps, :onestept, 3, # time after a step
        (v,P,δ)->onestep(w,P=P,vm=v, δangle=δ).tf, autodiff=true)
    @NLexpression(optsteps, summedtime, # add up time of all steps
        sum(onestept(v[i],P[i],δs[i]) for i = 1:numsteps))
    for i = 1:numsteps  # step dynamics
        @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i],δs[i]))
    end
    @NLconstraint(optsteps, summedtime == totaltime) # total time

    if boundarywork[1] # including boundary work
        @objective(optsteps, Min, 1/2*(sum(((comp[i]*P[i])^impwork for i=1:numsteps))+v[1]^2-boundaryvels[1]^2)+0*(v[end]^2-boundaryvels[2]^2)) # minimum pos work
    else 
        @objective(optsteps, Min, 1/2*sum(((comp[i]*P[i])^impwork for i=1:numsteps))) # minimum pos work alone
    end
    optimize!(optsteps)
    if termination_status(optsteps) == MOI.LOCALLY_SOLVED || termination_status(optsteps) == MOI.OPTIMAL
        optimal_solution = (vms=value.(v), Ps=value.(P))
    else
        error("The model was not solved correctly.")
        println(termination_status(optsteps))
    end

    return multistep(W(w,vm=value(v[1])), value.(P), δs, vm0=value(v[1]), boundaryvels=boundaryvels,
        extracost = boundarywork[1] ? 1/2*(value(v[1])^2 - boundaryvels[1]^2)+0/2*(value(v[end])^2-boundaryvels[2]^2) : 0) #, optimal_solution
end

wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
nsteps = 15
δs = zeros(nsteps); δs[Int((nsteps+1)/2)] = 0.05 # one bump
nominalmsr=optwalkcomp(wstar4s, nsteps, boundaryvels=(), boundarywork=false, δs=δs, horvert=3, impwork=2)
horizwmsr = optwalkcomp(wstar4s, nsteps, boundaryvels = (), boundarywork=false, δs=δs, horvert=1, impwork=2)
horizimsr = optwalkcomp(wstar4s, nsteps, boundaryvels = (), boundarywork=false, δs=δs, horvert=1, impwork=1)
vertwmsr = optwalkcomp(wstar4s, nsteps, boundaryvels = (), boundarywork=false, δs=δs, horvert=2, impwork=2)
vertimsr = optwalkcomp(wstar4s, nsteps, boundaryvels = (), boundarywork=false, δs=δs, horvert=2, impwork=1)
plotvees(nominalmsr,boundaryvels=nominalmsr.boundaryvels, speedtype=:midstance, usespline=false, legend=true, label =["H & V Work" ;])
plotvees!(horizwmsr, boundaryvels=nominalmsr.boundaryvels, speedtype=:midstance, usespline=false, legend=true,label = ["H work"])
plotvees!(horizimsr, boundaryvels=nominalmsr.boundaryvels, speedtype=:midstance, usespline=false, legend=true,label = ["H imp"])
savefig("h&v.png")
## Horizontal gives clearly the wrong kind of response (propulsion)
# Horizontal work gives a nearly opposite response, with work 0.46 or 1.0047% greater than nominal
# Horizontal impulse gives a square but opposite response, with work 1.033% greater
# Vertical resembles normal push-off, since leg is mostly vertical anyway


# now do a down-step
δs = zeros(nsteps); δs[Int((nsteps+1)/2)] = -0.05 # one bump
nominalmsr=optwalkcomp(wstar4s, nsteps, boundaryvels=(), boundarywork=false, δs=δs, horvert=3, impwork=2)
horizwmsr = optwalkcomp(wstar4s, nsteps, boundaryvels = (), boundarywork=false, δs=δs, horvert=1, impwork=2)
horizimsr = optwalkcomp(wstar4s, nsteps, boundaryvels = (), boundarywork=false, δs=δs, horvert=1, impwork=1)
vertwmsr = optwalkcomp(wstar4s, nsteps, boundaryvels = (), boundarywork=false, δs=δs, horvert=2, impwork=2)
vertimsr = optwalkcomp(wstar4s, nsteps, boundaryvels = (), boundarywork=false, δs=δs, horvert=2, impwork=1)
plotvees(nominalmsr,boundaryvels=nominalmsr.boundaryvels, speedtype=:midstance, usespline=false, legend=true, label =["H & V Work" ;])
plotvees!(horizwmsr, boundaryvels=nominalmsr.boundaryvels, speedtype=:midstance, usespline=false, legend=true,label = ["H work"])
plotvees!(horizimsr, boundaryvels=nominalmsr.boundaryvels, speedtype=:midstance, usespline=false, legend=true,label = ["H imp"])
