

using DynLoco, Plots, Statistics
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

wstar4 = findgait(WalkRW2l(α=0.35,safety=true), target=:speed=>0.3, varying=:P)

# take multiple steps with some deltas
nsteps = 4
δs = 0.02*rand(nsteps)
Ps = wstar4.P*ones(nsteps)
msr = multistep(wstar4, Ps = Ps, δangles=δs)


# This is taking multiple steps in a loop of onesteps, to get
# same result
vmprev = wstar4.vm # start with regular speed
vees = zeros(nsteps)
taus = zeros(nsteps)
for i in 1:nsteps
    osr = onestep(wstar4, vm=vmprev, P=Ps[i], δangle=δs[i])
    vees[i] = osr.vm
    taus[i] = osr.tf
    vmprev = osr.vm
    movingaveragev[i] = mean(vees[i-window:i])
    movingaveraget[i] = mean(taus[i-window:i])
end

# Include push-off selection
nh = 4
h = [-0.1, 0.2, 0., -0.1] # a sample impulse response
vees = zeros(nsteps)
taus = zeros(nsteps)
expectedaveragev = wstar4.vm
expectedaveraget = onestep(wstar4).tf
expectedsumwork = nh*wstar4.P^2
movingaveragev = zeros(nsteps)
movingaveraget = zeros(nsteps)
movingsumwork = zeros(nsteps)
errorv = zeros(nsteps)
errort = zeros(nsteps)
errorwork = zeros(nsteps)
for i in 1:nsteps
    vnext = conv(bumps[i:i+nh-1],h) # select the next speed
    # convert vnext to the appropriate push-off
    osr = onestepp(wstar4, vm=vmprev, vnext=vnext, δangle=δs[i])
    P = osr.P # the push-off needed to produce the desired step
    Phistory[i] = P
    movingaveragev[i] = mean(vees[i-window:i])
    movingaveraget[i] = mean(taus[i-window:i])
    movingsumwork[i] = sum(Phistory[i-window:i]^2)
    errorv[i] = movingaveragev[i] - expectedaveragev[i]
    errort[i] = movingaveraget[i] - expectedaveraget[i]
    errorwork[i] = movingsumwork[i] - expectedsumwork[i]

    vmprev = osr.vm
end



stepresult = onestep(wstar4, δangle=0.01)
onestep(wstar4, vm = stepresult.vm, δangle = stepresult.)

ctime = 0.015 # cost of time, to encourage hurrying
tchange = 1.75
p = plot()
walksteps = [1, 2, 3, 4, 5, 6, 7, 10, 15, 20] # take walks of this # of steps
results = Array{MultiStepResults,1}(undef,0) # store each optimization result here


for i in 1:nsteps
    



end


## onestepp computes the push-off and timing as result of speed
function onestepp(w::WalkRW2l; vm=w.vm, vnext=w.vm, δangle = 0.,
    α=w.α, γ=w.γ,g=w.g,L=w.L,safety=w.safety)
    # Start at mid-stance leg vertical, and find the time tf1
    # to heelstrike, and angular velocity Ωminus
    # δangle is the upward change in slope from nominal γ
    # Phase 1: From mid-stance to just before impact
    # θ is angle ccw from surface normal (which may be at slope)
    θf = -α + δangle # angle of stance leg just before impact, negative means cw
    # angle wrt vertical is -α-γ+δangle
    Ωminus = -√(2*g/L*(cos(-γ)-cos(θf-γ))+(vm/L)^2)  # pre-impact speed
    tf1 = log((L*(γ-θf)+√(vm^2-2γ*θf*L^2+θf^2*L^2))/(vm+L*γ)) # time until impact, phase 1

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

osr = onestep(wstar4)
ospr = onestepp(wstar4, vnext=wstar4.vm)
ospr.P

# let's check whether we can reproduce the previous steps
vmprev = wstar4.vm # start with regular speed
Ps .= Ps .+ randn(nsteps)*0.02
vees = zeros(nsteps)
taus = zeros(nsteps)
Pchecks = zeros(nsteps)
for i in 1:nsteps
    osr = onestep(wstar4, vm=vmprev, P=Ps[i], δangle=δs[i])
    ospr = onestepp(wstar4, vm=vmprev, vnext=osr.vm, δangle=δs[i])
    vees[i] = osr.vm
    taus[i] = osr.tf
    Pchecks[i] = ospr.P
    vmprev = osr.vm
end
Pchecks
Ps