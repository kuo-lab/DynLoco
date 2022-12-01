## Differentiating the rimless wheel

using DynLoco
using Plots, Statistics, ToeplitzMatrices, DSP
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

demeanandnormalize(response, normfactor, withmean=mean(response)) = (response .- withmean) ./ normfactor
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
tstar = t0 = onestep(wstar4s).tf # nominal step time
Pstar = wstar4s.P
vstar = wstar4s.vm
workstar = onestep(wstar4s).Pwork

nsteps = 15; halfwindow = (nsteps+1) ÷ 2
δou = zeros(nsteps); δou[Int((nsteps+1)/2)] = 0.05
oumsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δou)
h = demeanandnormalize(convert(Vector{Float64},oumsr.steps.vm),0.05) # normalized impulse response
hstar = copy(h)
plotvees(oumsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)


# I want to differentiate onestep, where we send it P, vm, and delta, and
# get back the next v and tf. I will want the derivative with respect to v?
# J = 1/2 P^2 + 1/2 tf^2
# dJ/dtheta = P * dP/dtheta + tf*dtf/dtheta where v is the desired speed
#           = P * dP/dv*dv/dtheta + tf * dtf/dv*dv/dtheta, where we want to treat v as the command
# dP/dv = sensitivity of P to a change in vnext
# dtf/dv = sensitivity of tf to a change in vnext

onestep(wstar4s, P = wstar4s.P*1.05)
onestepp(wstar4s)
#oofstep(P, )

onestepp(wstar4s, vm = wstar4s.vm, vnext = 0.4, δangle = 0.).P
using ForwardDiff

# take a derivative of P with respect to v
dPdv = ForwardDiff.derivative( v->onestepp(wstar4s, vnext=v).P, wstar4s.vm)# 1.001646
# which implies that if I change v a tiny bit, then P will change by amount
onestepp(wstar4s, vnext=wstar4s.vm+0.0001).P - wstar4s.P # 0.00010017059
# so it pretty much works

# take a derivative of tf with respect to v
# change v a little bit, then tf will change by amount
dtfdv = ForwardDiff.derivative( v->onestepp(wstar4s, vnext=v).tf, wstar4s.vm) # -1.727319
onestepp(wstar4s, vnext=wstar4s.vm+0.0001).tf - onestep(wstar4s).tf # -0.0001726995
# so it pretty much works

# now combine d(1/2P^2)/dtheta = sum P*dP/dv*dv/dtheta
h = hstar
P = onestepp(wstar4s, δangle = 0.).P
dPdv = ForwardDiff.derivative( v->onestepp(wstar4s, vnext=v).P, wstar4s.vm)# 1.001646
dvdtheta = h
dworkdtheta = P * dPdv .* dvdtheta # a gradient with respect to parameters
tf = onestepp(wstar4s, δangle = 0.).tf
dtfdv = ForwardDiff.derivative( v->onestepp(wstar4s, vnext=v).tf, wstar4s.vm) # -1.727319
dtimedtheta = (tf-tstar) * dtfdv .* dvdtheta



# I'm trying to compute the derivative by chain rule
vees = zeros(nterrain); vcheck2 = zeros(nterrain); Ps = zeros(nterrain); taus = zeros(nterrain)
vmprev = vstar; movingtau = zeros(nterrain); movingwork = zeros(nterrain); workcorr = zeros(nterrain,nsteps); taucorr = zeros(nterrain,nsteps)
dworkdtheta = zeros(nsteps); dtimedtheta = zeros(nsteps)
for i in (1:nsteps).+8
    predictedv = sum(reverse(padme(δs,i-halfwindow+1,i+halfwindow-1)) .* h)
#    predictedv = sum(reverse(h) .* padme(δs,i-halfwindow+1,i+halfwindow-1))
    osr = onestepp(wstar4s, vm=vmprev, vnext=predictedv+vstar, δangle=δs[i]) 
    osr2 = onestep(wstar4s, vm=vmprev, P=osr.P, δangle=δs[i])
    vees[i] = predictedv
    Ps[i] = osr.P
    taus[i] = osr.tf
    vmprev = predictedv + vstar
    vcheck2[i] = osr2.vm
    dPdv = ForwardDiff.derivative( v-> onestepp(wstar4s, vm=vmprev, vnext=v, δangle=δs[i]).P, predictedv+vstar)
    dtfdv = ForwardDiff.derivative( v-> onestepp(wstar4s, vm=vmprev, vnext=v, δangle=δs[i]).tf, predictedv+vstar)
    dvdtheta = reverse(padme(δs,i-halfwindow+1,i+halfwindow-1))
    dworkdtheta .+= osr.P * dPdv .* dvdtheta
    dtimedtheta .+= (osr.tf-tstar) * dtfdv .* dvdtheta

    movingtau[i] = 1/nsteps*sum(padme(taus,i-nsteps+1,i)) - tstar
    movingwork[i] = 1/nsteps*sum(abs2,0.5*padme(Ps,i-2+1,i)) - 0.5*Pstar^2
    if i >= nsteps 
        workcorr[i,:] = movingwork[i] * padme(δs,i-nsteps+1,i) # should complain about high work
        taucorr[i,:] = movingtau[i] * padme(δs, i-nsteps+1,i)
    end
end
plot(vs0); plot!(vees); plot!(vcheck2.-vstar) # good


# Meanwhile try to compute derivative by making a function
# that computes work and time for the terrain
function takesteps(w::WalkRW2l;δs=δs, h=hstar)
# return the work and time for taking a bunch of steps on terrain
# using policy h    
vees = zeros(nterrain); vcheck2 = zeros(nterrain); Ps = zeros(nterrain); taus = zeros(nterrain)
vmprev = vstar; movingtau = zeros(nterrain); movingwork = zeros(nterrain); workcorr = zeros(nterrain,nsteps); taucorr = zeros(nterrain,nsteps)
dworkdtheta = zeros(nsteps); dtimedtheta = zeros(nsteps)
for i in (1:nsteps).+8
    predictedv = sum(reverse(padme(δs,i-halfwindow+1,i+halfwindow-1)) .* h)
#    predictedv = sum(reverse(h) .* padme(δs,i-halfwindow+1,i+halfwindow-1))
    osr = onestepp(wstar4s, vm=vmprev, vnext=predictedv+vstar, δangle=δs[i]) 
    osr2 = onestep(wstar4s, vm=vmprev, P=osr.P, δangle=δs[i])
    vees[i] = predictedv
    Ps[i] = osr.P
    taus[i] = osr.tf
    vmprev = predictedv + vstar
    vcheck2[i] = osr2.vm
    dPdv = ForwardDiff.derivative( v-> onestepp(wstar4s, vm=vmprev, vnext=v, δangle=δs[i]).P, predictedv+vstar)
    dtfdv = ForwardDiff.derivative( v-> onestepp(wstar4s, vm=vmprev, vnext=v, δangle=δs[i]).tf, predictedv+vstar)
    dvdtheta = reverse(padme(δs,i-halfwindow+1,i+halfwindow-1))
    dworkdtheta .+= osr.P * dPdv .* dvdtheta
    dtimedtheta .+= (osr.tf-tstar) * dtfdv .* dvdtheta

    movingtau[i] = 1/nsteps*sum(padme(taus,i-nsteps+1,i)) - tstar
    movingwork[i] = 1/nsteps*sum(abs2,0.5*padme(Ps,i-2+1,i)) - 0.5*Pstar^2
    if i >= nsteps 
        workcorr[i,:] = movingwork[i] * padme(δs,i-nsteps+1,i) # should complain about high work
        taucorr[i,:] = movingtau[i] * padme(δs, i-nsteps+1,i)
    end
end
plot(vs0); plot!(vees); plot!(vcheck2.-vstar) # good











# onestepptf returns the P needed to get the desired tf
"""
  `onestepp(w::Walk; vm=0.4, vnext=0.45, δangle = 0.)`

  Takes a step given a desired mid-stance velocity `vnext`, computing
  the push-off `P` needed. Returns a named tuple of step info, including
  step time `tf` and other info. The step starts at mid-stance speed `vm`
  defined as time of surface normal.
"""
function onestepptf(w::WalkRW2l; vm=w.vm, vnext=w.vm, δangle = 0., destf=0., 
    α=w.α, γ=w.γ,g=w.g,L=w.L,safety=w.safety)
    # Start at mid-stance leg vertical, and find the time tf1
    # to heelstrike, and angular velocity Ωminus
    # δangle is the upward change in slope from nominal γ
    # Phase 1: From mid-stance to just before impact
    # θ is angle ccw from surface normal (which may be at slope)

    if destf == 0. # zero desired time, so fill in the default
        destf = onestepp(w, vm=vm, vnext=vnext,δangle=δangle).tf
    end

    θf = -α + δangle # angle of stance leg just before impact, negative means cw
    # angle wrt vertical is -α-γ+δangle
    Ωminus = -√(2*g/L*(cos(-γ)-cos(θf-γ))+(vm/L)^2)  # pre-impact speed
    tf1 = log((L*(γ-θf)+√(vm^2-2γ*θf*L^2+θf^2*L^2))/(vm+L*γ)) # time until impact, phase 1

    remainingtime = destf - tf1 # this is the desired time, need to find the P that produces this time


    # Phase 2: From just after impact to next mid-stance
    θnew = α + δangle # new stance leg's initial angle
    Ωplus = -sqrt(vnext^2 - 2g*L*cos(θnew-γ) + 2g*L*cos(-γ))/L

    # Step-to-step transition: Push-off and collision
    #Ωplus = cos(2α)*Ωminus - P*sin(2α) # initial ang vel of next leg
    #P*sin(2α) = cos(2α)*Ωminus - Ωplus 
    P = (cos(2α)*Ωminus - Ωplus)/sin(2α)
    Pwork = 1/2 * P^2
    v0 = √(Ωminus^2 + P^2) # intermediate velocity after push-off
    Cwork = 1/2*(Ωplus^2 - v0^2) # negative collision work
    C = √(-2Cwork)
    #    tf2 = mylog((γ + √(2γ*θnew-θnew^2+Ωplus^2))/(γ-θnew-Ωplus))
    gto = γ - θnew - Ωplus # needs to be positive for pendulum to reach next mid-stance
    if safety && gto <= 0
        tf2 = 1e3 - gto # more negative gto is worse, increasing a tf2 time
    else # safety off OR gto positive
        # time to mid-stance, phase 2
        inroot = (2γ*θnew-θnew^2+Ωplus^2)
        if inroot >= 0
#            tf2 = mylog((γ + √(2γ*θnew-θnew^2+Ωplus^2))/(gto)) # time to mid-stance, phase 2
            tf2 = log((γ + √inroot)/(gto)) # time to mid-stance, phase 2
        else # inroot negative,
            tf2 = 1e4 - inroot # more negative inroot extends time
        end

    end

# derivation of omegaplus, using old equations
    #twicenextenergy = (2g*L*cos(θnew-γ)+L^2*Ωplus^2-2g*L*cos(-γ)) # to find next vm
    #0.5*vmnew^2 = (g*L*cos(θnew-γ)+0.5*L^2*Ωplus^2-g*L*cos(-γ))
    #KEupright   =   PEplus + KEplus             -PEupright
    #KEplus = KEupright - PEplus + PEupright
    #0.5*L^2*Ωplus^2 = 0.5*vmnew^2 - g*L*cos(θnew-γ) + g*L*cos(-γ)
    #L^2*Ωplus^2 = vmnew^2 - 2g*L*cos(θnew-γ) + 2g*L*cos(-γ)
    #Ωplus^2 = (vmnew^2 - 2g*L*cos(θnew-γ) + 2g*L*cos(-γ))/L^2
    #Ωplus = sqrt(vmnew^2 - 2g*L*cos(θnew-γ) + 2g*L*cos(-γ))/L

    # Step metrics
    steplength = 2L*sin(α) # rimless wheel step length
    tf = tf1 + tf2 # total time mid-stance to mid-stance
    speed = steplength / tf # average speed mid-stance to mid-stance
    return (vm=vnext, θnew=θnew, tf=tf, P=P, C=C, Pwork=Pwork,Cwork=Cwork,
        speed=speed, steplength=steplength, stepfrequency=speed/steplength,tf1=tf1, tf2=tf2,
        Ωminus=Ωminus,Ωplus=Ωplus,
        vm0=vm,δ=δangle)
end


#v+ = v-*cos + P * sin(2a)

#pendulum isa
#theta(tf) = (theta+v+)*e^-t/tau + (theta-v+)B*e^t/tau = 0
# solve for v+

# if you ask for a different v, we will produce a different P, t
# if you ask for a different t, we will produce a different P, vnext
 1/2 P^2 + dP/dtheta + 1/2tf* dtf/dtheta