## Travelator or moving walkway

# let's first optimize without time
# the goal is to start walking at nominal speed, and then 
# meet the travelator and walk at nominal speed on it

using DynLoco, JuMP, Ipopt
using Plots, Statistics, ToeplitzMatrices, DSP
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

"""
    optwalkwalkway(w::Walk, numsteps=5)

Optimizes push-offs to walk `numsteps` steps. Returns a `MultiStepResults`
struct. Allows slopes to be specified as keyword `δs` array of the slope of each successive
step.

Walkway version assumes a moving walkway at step (n+1)/2, with speed v0

Other keyword arguments: `boundaryvels = (vm,vm)` can specify a tuple of initial and
final speeds, default nominal middle-stance speed `vm`. To start and end at rest, use `(0,0)`.
`boundarywork = true` whether cost includes the work needed to
start and end from `boundaryvels`. `totaltime` is the total time to take the steps, by default
the same time expected of the nominal `w` on level ground.

"""
function optwalkwalkway(w::W, numsteps=5; vwalkway=0., boundaryvels::Union{Tuple,Nothing} = nothing,
    boundarywork::Union{Tuple{Bool,Bool},Bool} = (true,true), totaltime=numsteps*onestep(w).tf,
    δs = zeros(numsteps)) where W <: Walk # default to taking the time of regular steady walking

    nwalkway = (numsteps+1) ÷ 2
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
        (v,P,δ)->onestepwalkway(w,P=P,vm=v, δangle=δ,vwalkway=0.).vm, autodiff=true) # output vm
    register(optsteps, :onestept, 3, # time after a step
        (v,P,δ)->onestepwalkway(w,P=P,vm=v, δangle=δ,vwalkway=0.).tf, autodiff=true)
    register(optsteps, :onestepvw, 3, # velocity after a step
        (v,P,δ)->onestepwalkway(w,P=P,vm=v, δangle=δ,vwalkway=vwalkway).vm, autodiff=true) # output vm
    register(optsteps, :onesteptw, 3, # time after a step
        (v,P,δ)->onestepwalkway(w,P=P,vm=v, δangle=δ,vwalkway=vwalkway).tf, autodiff=true)
    @NLexpression(optsteps, summedtime, # add up time of all steps
        sum(onestept(v[i],P[i],δs[i]) for i = 1:nwalkway-1) + sum(onestept(v[i],P[i],δs[i]) for i = nwalkway+1:numsteps) + 
           onesteptw(v[nwalkway],P[nwalkway],δs[nwalkway]))
    for i = 1:numsteps  # step dynamics
        if i == nwalkway 
            @NLconstraint(optsteps, v[i+1] == onestepvw(v[i],P[i],δs[i]))
        else
            @NLconstraint(optsteps, v[i+1]==onestepv(v[i],P[i],δs[i]))
        end
    end
    @NLconstraint(optsteps, summedtime == totaltime) # total time

    if boundarywork[1]
        @objective(optsteps, Min, 1/2*(sum((P[i]^2 for i=1:numsteps))+v[1]^2-boundaryvels[1]^2)+0*(v[end]^2-boundaryvels[2]^2)) # minimum pos work
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

    return multistep(W(w,vm=value(v[1])), value.(P), δs, vm0=value(v[1]), boundaryvels=boundaryvels,
        extracost = boundarywork[1] ? 1/2*(value(v[1])^2 - boundaryvels[1]^2)+0/2*(value(v[end])^2-boundaryvels[2]^2) : 0) #, optimal_solution
end
# normally vplus = vminus*cos(2a) + P*sin(2a) , but in this case
# vplus + v0*cos(α) == vminus*cos(2α) + P*sin(2α)

"""
    onestepwalkway(walk [,parameters, safety])

Take one step with a `walk` struct, from mid-stance to mid-stance. Parameter
values can be supplied to override the existing field values. Notably,
includes vm, P, δangle.
Walkway version has a keyword arg `vwalkway=0` for stepping onto walkway at that speed.
Output is a `NamedTuple` of outcomes such as `vm` (at end), `tf` (final time),
`work`, `speed`, `steplength`, etc.
Set `safety=true` to prevent step from blowing up if stance leg has too little
momentum to reach middle stance. A step will be returned with very long
step time, and very slow mid-stance velocity.
"""
function onestepwalkway(w::WalkRW2l; vm=w.vm, P=w.P, vwalkway=0.,δangle = 0.,
    α=w.α, γ=w.γ,g=w.g,L=w.L,safety=w.safety)
    mylog = safety ? DynLoco.logshave : log # logshave doesn't blow up on negative numbers
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
    Ωplus = cos(2α)*Ωminus - P*sin(2α) + vwalkway*cos(α) # initial ang vel of next leg (Ω should be negative)
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


wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
tstar = t0 = onestep(wstar4s).tf # nominal step time
Pstar = wstar4s.P
vstar = wstar4s.vm
workstar = onestep(wstar4s).Pwork

nsteps = 15; 
δlevel = zeros(nsteps);
slowwalkway=optwalkwalkway(wstar4s, nsteps, boundarywork=false,vwalkway=0.1)
plotvees(slowwalkway,speedtype=:midstance,boundaryvels=(wstar4s.vm,wstar4s.vm))





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