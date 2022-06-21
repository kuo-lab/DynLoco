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
plotvees(oumsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)

# downstep
δod = zeros(nsteps); δod[Int((nsteps+1)/2)] = -0.05
odmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δod)
hd = demeanandnormalize(convert(Vector{Float64},odmsr.steps.vm),0.05) # normalized down impulse response
plotvees(odmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)

# up & down-step
δud = zeros(nsteps); δud[Int((nsteps+1)/2)] = 0.05; δud[Int((nsteps+1)/2+1)] = -0.05; 
udmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δud)
hud = demeanandnormalize(convert(Vector{Float64},udmsr.steps.vm),0.05)
plotvees(udmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)

## Simulate a stretch of uneven terrain
nterrain = 100
δs = rand(nterrain)./10
coefs = hcat(0:nterrain-1,ones(nterrain))\δs # slope and offset
δs .= δs .- hcat(0:nterrain-1,ones(nterrain))*coefs # detrend the bumps
nominalmsr=optwalk(wstar4s, nterrain, boundarywork=false, δs=δs)
pfig = plotvees(nominalmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)
display(pfig)
vs = convert(Vector{Float64}, nominalmsr.steps.vm)
vs0 = vs .- nominalmsr.vm0
verror = vs0
ts = convert(Vector{Float64}, nominalmsr.steps.tf)
ts0 = ts .- t0
terror = ts0
pushoffs = nominalmsr.steps.P

## Assemble bumps and impulse response
# Make a square Toeplitz matrix from the bumps
A = TriangularToeplitz(δs,:L)
pv = A[:,1:nsteps]*h # also known as conv(δs,h)
pv = pv[halfwindow:end] # remove first steps
plot(pv,label="predicted v")
plot!(vs0, label="demeaned v*") # so linearly predicted v resembles computed optimum
# So a convolution of bumps and impulse response
# does a good job of predicting the velocities

# so A*h = (v-mean(v)), and you can also go backward, yielding the impulse response
# we can invert the A (bump) matrix to get the impulse response
plot(h,label="optimal up",title="Optimal h and regressed h")
plot!(A[halfwindow:end,1:2*halfwindow-1] \ vs0[1:end-halfwindow+1],label="pinv response")



## Take a moving average of speeds
movingavekernel = ones(nsteps)./nsteps
movingavev = conv(movingavekernel, vs0)[halfwindow+0:end+1-halfwindow] # moving average speed error
plot(movingavev, title="moving average speed")
plot!(vs0)

## How well correlated is the moving average with the bumps?
plot(movingavev .* circshift(δs,20))
plot!(movingavev .* δs)


## we want to do gradient descent, where the
# error in speed, v-v0 is correlated with the bumps 

matrixofcorr = zeros(nterrain,nsteps)
newh = h.*0
mu = 0.01
pfig = plot(newh)
for i in halfwindow:nterrain-halfwindow+1
    matrixofcorr[i,:] .= vs0[i-1]*δs[i-halfwindow+1:i+halfwindow-1]
    newh .= newh .- mu*(vs0[i-1]*δs[i-halfwindow+1:i+halfwindow-1])
    plot!(pfig, newh)
end
display(pfig)
plot(newh*1000)
plot!(h)

## now take some steps, replicating the optimum
Ps = zeros(nterrain); steptimes = zeros(nterrain); vchecks = zeros(nterrain)
vmprev = vstar
for i in eachindex(δs)
    predictedv = vs[i]
    stepresult = onestepp(wstar4s; vm=vmprev, vnext=predictedv, δangle=δs[i])
    Ps[i] = stepresult.P
    steptimes[i] = stepresult.tf
    stepresult = onestep(wstar4s, vm=vmprev, P=Ps[i], δangle=δs[i])
    vchecks[i] = stepresult.vm
    vmprev = stepresult.vm
end
plot(vs)
plot!(vchecks,linestyle=:dash)
plot(Ps)
plot!(nominalmsr.steps.P,linestyle=:dash)

## another replication of optimum, but using convolution
Ps = zeros(nterrain); steptimes = zeros(nterrain); vchecks = zeros(nterrain)
vmprev = vstar
for i in 8:length(δs)-7
    predictedv = sum(reverse(δs[i-halfwindow+1:i+halfwindow-1]) .* h)
    vchecks[i] = predictedv
end
plot(vchecks)
plot!(vs0)


## Use convolution to predict the next v and P
Ps = zeros(nterrain); steptimes = zeros(nterrain); vchecks = zeros(nterrain); vchecks2 = zeros(nterrain)
terror = zeros(nterrain); workerror = zeros(nterrain)
vmprev = vstar
δpad = vcat(zeros(halfwindow-1), δs, zeros(halfwindow-1))
newh = h#(rand(nsteps) .- 0.5)*0.25#h.*0
muw = 0.5; mut = 0.5
for i in 8:length(δpad)-7
    predictedv = sum(reverse(δpad[i-halfwindow+1:i+halfwindow-1]) .* h)
    vchecks[i-7] = predictedv
    stepresult = onestepp(wstar4s; vm=vmprev, vnext=predictedv+vstar, δangle=δpad[i])
    P = stepresult.P # based on our current impulse response, this is the push-off to apply
    Ps[i-7] = P
    steptimes[i-7] = stepresult.tf
    stepresult = onestep(wstar4s, vm=vmprev, P=P, δangle=δpad[i])
    vchecks2[i-7] = stepresult.vm
    vmprev = stepresult.vm
    terror[i-7] = stepresult.tf - tstar
    workerror[i-7] = stepresult.Pwork - workstar
    newh .= newh .- muw*(workerror[i-7]*δpad[i-halfwindow+1:i+halfwindow-1]) .-
      mut*(terror[i-7]*δpad[i-halfwindow+1:i+halfwindow-1])
end
plot(vs0)
plot!(vchecks, linestyle=:dash) # this is okay
plot!(vchecks2.-vstar, linestyle=:dash)
plot(Ps)
plot!(nominalmsr.steps.P,linestyle=:dash)
plot(newh)



matrixofcorr = zeros(nterrain,nsteps)
newh = h.*0
muw = 0.01; mut = 0.01
vmprev = vstar
vs = zeros(nterrain); Ps = zeros(nterrain); steptimes = zeros(nterrain); terror = zeros(nterrain);
workerror = zeros(nterrain)
for i in eachindex(δs)
    window = max(1, i-halfwindow+1):min(length(δs),i+halfwindow-1)
    predictedv = conv(δs[window], newh)[nsteps-1] + wstar4s.vm
    stepresult = onestepp(wstar4s, vm=vmprev, vnext=predictedv, δangle=δs[i])
    P, tf = (stepresult.P, stepresult.tf)
    vs[i] = predictedv
    Ps[i] = P
    steptimes[i] = tf
    terror[i] = tf - tstar
    workerror[i] = stepresult.Pwork - workstar
    newh .= newh .- muw*(workerror[i-1]*δs[i-halfwindow+1:i+halfwindow-1]) .-
      mut*(terror[i-1]*δs[i-halfwindow+1:i+halfwindow-1])
end
plot(newh*1000)
plot!(h)

## Now do the same thing, except take a walk using the h0
vmprev = wstar4s.vm
for i in 1:nterrain
    # calculate the next v using the impulse response
    # transform the v into appropriate P, using onestepp
    stepresult = onestepp(wstar4s, vm=vmprev, vnext=predictedv, δangle=δs[i])
    # apply that step
    # evaluate the time, speed, energy cost
    P, tf, Pwork = (stepresult.P, stepresult.tf, stepresult.Pwork)
    # correlate time error with bump, timeerror should be a moving average
    terror * δs 
    # correlate energy error with bump, workerror should be a moving average
    workerror * δs
    # make adjustment to the impulse response with
    # stochastic gradient descent
    osr = onestep(wstar4, vm=vmprev, P=Ps[i], δangle=δs[i])
    vees[i] = osr.vm
    taus[i] = osr.tf
    vmprev = osr.vm
    movingaveragev[i] = mean(vees[i-window:i])
    movingaveraget[i] = mean(taus[i-window:i])
end

## 
newh = h * 0
halfwin = (nsteps-1) ÷ 2
vpred = conv(δs,h) # a prediction of the speeds
vmprev = wstar4s.vm
for i in halfwin+1:length(vpred)-halfwin
    osr = onestepp(wstar4s, vm=vmprev, vnext = vpred[i], δangle=δ[])
end

#for i in eachindex(δs)
i = 1
oneclumpv = conv(padme(δs, i-halfwin, i+halfwin), h)[1+halfwin:end-halfwin]


Ps = zeros(nterrain); taus = zeros(nterrain); vees = zeros(nterrain); vcheck2 = zeros(nterrain)
vmprev = vstar; halfwin = 7
j = 1 # where we are in δs
for clump in 1:6
    # take some steps 
    oneclumpv = conv(padme(δs, j-halfwin, j+halfwin), h)[1+halfwin:end-halfwin]
    vmprev = conv(padme(δs, j-halfwin, j+halfwin), h)[halfwin] + vstar
    for i in 1:nsteps
        osr = onestepp(wstar4s, vm=vmprev, vnext=oneclumpv[i]+vstar, δangle=δs[j]) # this doesn't work
        #osr = onestepp(wstar4s, vm=vmprev, vnext=vs[j], δangle=δs[j])
        vees[j] = oneclumpv[i] + vstar
        Ps[j] = osr.P
        osr2 = onestep(wstar4s, vm=vmprev, P=Ps[j], δangle=δs[j])
        vcheck2[j] = osr2.vm
        taus[j] = osr.tf
        vmprev = oneclumpv[i]+vstar
        j = j + 1
    end
end
plot(vs); plot!(vees[halfwin+1:end]); plot!(vcheck2[halfwin+1:end])
plot(pushoffs); plot!(Ps[halfwin+1:end])


plot(conv(padme(δs,i-halfwin,i+halfwin),h))
plot!(conv(padme(δs,i-halfwin,i+halfwin),h)[1+halfwin:end-halfwin])
#end

## take a bunch of single steps (simulation works)
vees = zeros(nterrain); vcheck2 = zeros(nterrain); Ps = zeros(nterrain); taus = zeros(nterrain)
movingtau = zeros(nterrain); movingwork = zeros(nterrain); workcorr = zeros(nterrain,nsteps); taucorr = zeros(nterrain,nsteps)
for i in 1:nterrain
    predictedv = sum(reverse(padme(δs,i-halfwindow+1,i+halfwindow-1)) .* h)
    osr = onestepp(wstar4s, vm=vmprev, vnext=predictedv+vstar, δangle=δs[i]) 
    osr2 = onestep(wstar4s, vm=vmprev, P=osr.P, δangle=δs[i])
    vees[i] = predictedv
    Ps[i] = osr.P
    taus[i] = osr.tf
    vmprev = predictedv + vstar
    vcheck2[i] = osr2.vm
    movingtau[i] = 1/nsteps*sum(padme(taus,i-nsteps+1,i)) - tstar
    movingwork[i] = 1/nsteps*sum(abs2,0.5*padme(Ps,i-2+1,i)) - 0.5*Pstar^2
    if i >= nsteps 
        workcorr[i,:] = movingwork[i] * padme(δs,i-nsteps+1,i) # should complain about high work
        taucorr[i,:] = movingtau[i] * padme(δs, i-nsteps+1,i)
    end
end
plot(vs0); plot!(vees); plot!(vcheck2.-vstar) # good
plot(pushoffs); plot!(Ps) # also good
plot(nominalmsr.steps.tf); plot!(taus); plot!(movingtau)

plot(plot(mean(workcorr,dims=1)',title="workcorr"),
plot(mean(taucorr,dims=1)',title="taucorr"))
# when workcorr is high, those are weights that should be reduced to
# reduce work. When taucorr is high, those are weights that should be
# reduced to reduce time
plot(mean(workcorr,dims=1)')
plot!(0.2*mean(taucorr,dims=1)')



msr = multistep(wstar4s, Ps=nominalmsr.steps.P, δangles=δs)
plot(vs); plot!(msr.steps.vm)
msr.steps[1]
onestep(wstar4s, vm=vstar, P=pushoffs[1], δangle=δs[1])

vtest = zeros(length(δs)); for i in eachindex(δs); osr = onestep(wstar4s, vm=vmprev, P=pushoffs[i], δangle=δs[i]); vtest[i] = osr.vm; vmprev=osr.vm; end
plot(vs); plot!(vtest)
Ptest = zeros(length(δs)); for i in eachindex(δs); osr = onestepp(wstar4s, vm=vmprev, vnext=vs[i], δangle=δs[i]); Ptest[i] = osr.P; vmprev=osr.vm; end
plot(pushoffs); plot!(Ptest)
plot([1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]*h)
h0 = h .- mean(h)
Aa = [ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 
      -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 -1 1 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 -1 1 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 -1 1 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 -1 1 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 -1 1];
plot(Aa*h0)
plot!(hud.-mean(hud))



# so toeplitz*h = v



# we think there is a sequence h * b = v that will
# produce the desired speed sequence. 
# convolution is sum h(i)*b(n-i) or sum b(i)*h(n-i)
# or [b1 0   0] * [h1]    [h1 h2 h3] * [b1]
#    [b2 b1  0]   [h2]
#    [b3 b2 b1]   [h3]

# v = b1*h + b2*h(t-1) + b3*h(t-2) + ...
# v = [b1 0  0 0 0 0 0 0] * [h]
#     [b2 b1 0 0 0 0 0 0] * [h]
#     [b3 b2 b1 0 0 0 0 0] * [h]
#     [b4 b3 b2 b1 0 0 0 0] * [h]
#     [b5 b4 b3 b2]
#     [b6]
#     [b7]
#     [b8]
#     [b9 b8 b7]

# integral from 0 to t of u(tau)*h(t-tau) dtau


#


# moving window average
function movingwindow(x,n)
    firsthalf = floor(Int, (n+1)/2)
    remainder = n - firsthalf
    y = zeros(length(x))
    for i in firsthalf:length(x)-remainder
        y[i] = mean(x[(i-firsthalf+1):(i+remainder)])
    end
    return y
end
movingwindow(1:10,2)
plot(movingwindow(1:100,20))

function movingsum(x,n)
    firsthalf = floor(Int, (n+1)/2)
    remainder = n - firsthalf
    y = zeros(length(x))
    for i in firsthalf:length(x)-remainder
        y[i] = sum(x[(i-firsthalf+1):(i+remainder)])
    end
    return y
end

using DSP
function moving_average(x::AbstractArray,M::Integer)
    M > 1 || error("window size must be larger that 1")
    conv(x,rect(M)/M)[M ÷ 2 + 1 : length(x) + M÷2]
end

plot(moving_average(1:100, 20))


# now let's see what the moving window average of 
# speeds looks like
plot(vs0)
plot!(movingwindow(vs0,15))
# so the moving window speed is pretty small

# and same with the cumulative timings
plot(movingsum(nominalmsr.steps.tf,15))
plot!(ones(length(δs))*onestep(wstar4s).tf*15)
# and timing for moving average is close to nominal



## I need to take multiple steps here
# where multi-step is 
#multistep(wstar4s, Ps=X, δangles=X,

#)

using ForwardDiff

onestep(wstar4s,P=0.1,δangle=0).vm
nextv = P -> onestep(wstar4s,P=P).vm
nextv(0.1)
ForwardDiff.derivative(nextv, 0.1)

multistep(wstar4s, Ps=[0.1,0.2]).totaltime
costandtime(Ps) = (msr->[msr.totalcost,msr.totaltime])(multistep(wstar4s,Ps=Ps))
ForwardDiff.gradient( Ps->multistep(wstar4s, Ps=Ps).totaltime, [0.1, 0.2])
ForwardDiff.jacobian(costandtime,Ps=[0.1,0.2])
# and I can evaluate a Jacobian
ForwardDiff.jacobian( Ps->multistep(wstar4s,Ps=Ps))


crap(v) = [v[1]^2,2*v[2]]
crap([1,2])
ForwardDiff.jacobian(crap,[1,2])
candt(Ps) = (msr->[msr.totalcost,msr.totaltime])(multistep(wstar4s,Ps=Ps))
candt([0.1,0.2])
ForwardDiff.jacobian(candt,[0.1,0.2])
ForwardDiff.jacobian(candt,[0.1,0.17,0.1])
multistep(wstar4s,Ps=[0.1,0.15,0.18]).totalcost


"""
  `onestepp(w::Walk; vm=0.4, vnext=0.45, δangle = 0.)`

  Takes a step given a desired mid-stance velocity `vnext`, computing
  the push-off `P` needed. Returns a named tuple of step info, including
  step time `tf` and other info. The step starts at mid-stance speed `vm`
  defined as time of surface normal.
"""
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


"""
  `padme(array, from, to)``

  zero-padded indexing of `array`, allowing equivalent of `array[from:to]`
  where padded by zeros if `from <= 1` or `to >= length(array)`.
  For example, `padme(1:5, -2, 7)` automatically pads by 2 on each end.
"""
function padme(x::AbstractArray, from, to) 
    N = length(x)
    middlelength = to - from + 1 # length of middle section
    numberbefore1 = max(0, 1-from)
    numberafterN = max(0, to-N)
    return vcat(zeros(numberbefore1), x[max(1,from):min(N,to)], zeros(numberafterN))
end


padme([1,2,3],1,3)


