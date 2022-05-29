using DynLoco
using Plots, Statistics, ToeplitzMatrices, DSP
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

demeanandnormalize(response, normfactor, withmean=mean(response)) = (response .- withmean) ./ normfactor
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
t0 = onestep(wstar4s).tf # nominal step time
nsteps = 15; window = nsteps; halfwindow = (window+1) ÷ 2
δou = zeros(nsteps); δou[Int((nsteps+1)/2)] = 0.05
oumsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δou)
h = demeanandnormalize(convert(Vector{Float64},oumsr.steps.vm),0.05) # normalized impulse response
plotvees(oumsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)

δod = zeros(nsteps); δod[Int((nsteps+1)/2)] = -0.05
odmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δod)
hd = demeanandnormalize(convert(Vector{Float64},odmsr.steps.vm),0.05) # normalized down impulse response
plotvees(odmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)

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
plotvees(nominalmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)
vs = convert(Vector{Float64}, nominalmsr.steps.vm)
vs0 = vs .- nominalmsr.vm0
ts = convert(Vector{Float64}, nominalmsr.steps.tf)
ts0 = ts .- t0

## Assemble bumps and impulse response
# Make a square Toeplitz matrix from the bumps
A = TriangularToeplitz(δs,:L)
pv = A[:,1:nsteps]*h
pv = pv[halfwindow:end] # remove first steps
plot(pv,label="predicted v")
plot!(vs0, label="demeaned v*") # so linearly predicted v resembles computed optimum
# So a convolution of bumps and impulse response
# does a good job of predicting the velocities

# so A*h = (v-mean(v)), and you can also go backward, yielding the impulse response
# we can invert the A (bump) matrix to get the impulse response
plot(h,label="optimal up",title="Optimal h and regressed h")
plot!(A[halfwindow:end,1:window] \ vs0[1:end-halfwindow+1],label="pinv response")


## Take a moving average of speeds
movingavekernel = ones(window)./window
movingavev = conv(movingavekernel, vs0)[halfwindow+0:end+1-halfwindow] # moving average speed error
plot(movingavev, title="moving average speed")
plot!(vs0)

## How well correlated is the moving average with the bumps?
plot(movingavev .* circshift(δs,20))
plot!(movingavev .* δs)


## we want to do gradient descent, where the
# error in speed, v-v0 is correlated with the bumps

matrixofcorr = zeros(nterrain,window)
newh = h.*0
mu = 0.01
for i in halfwindow:nterrain-halfwindow+1
    matrixofcorr[i,:] .= vs0[i]*δs[i-halfwindow+1:i+halfwindow-1]
    newh .= newh .- mu*(vs0[i-1]*δs[i-halfwindow+1:i+halfwindow-1])
end
plot(newh*1000) # it has roughly the right shape
plot!(h)

## Now do the same thing, except take a walk using the h0
vmprev = wstar4s.vm
for i in 1:nterrain
    # calculate the next v using the impulse response
    # transform the v into appropriate P, using onestepp
    # apply that step
    # evaluate the time, speed, energy cost
    # correlate time error with bump
    # correlate energy error with bump
    # make adjustment to the impulse response with
    # stochastic gradient descent
    osr = onestep(wstar4, vm=vmprev, P=Ps[i], δangle=δs[i])
    vees[i] = osr.vm
    taus[i] = osr.tf
    vmprev = osr.vm
    movingaveragev[i] = mean(vees[i-window:i])
    movingaveraget[i] = mean(taus[i-window:i])
end



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

at=Toeplitz([1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

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