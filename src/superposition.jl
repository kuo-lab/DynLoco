using DynLoco
using Plots, Statistics, ToeplitzMatrices
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

demeanandnormalize(response, normfactor, withmean=mean(response)) = (response .- withmean) ./ normfactor
wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
nsteps = 15
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



plot([hd[1]; hd[1:end-1]])
hlate = [hd[1]; hd[1:end-1]]
plot(hlate.-mean(hlate))
plot!(h+hlate.-mean(h+hlate))

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



nsteps = 100
δs = rand(nsteps)./10
coefs = hcat(0:nsteps-1,ones(nsteps))\δs # slope and offset
δs .= δs .- hcat(0:nsteps-1,ones(nsteps))*coefs # detrend the bumps
nominalmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δs)
plotvees(nominalmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)
vs = convert(Vector{Float64}, nominalmsr.steps.vm)
vs0 = vs .- nominalmsr.vm0
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

# Make a square Toeplitz matrix from the bumps
A = TriangularToeplitz(δs,:L)
plot(A[:,1:15]*h)
pv = A[:,1:15]*h
pv = pv[8:end] # remove first steps

plot(pv,label="predicted v")
plot!(vs0, label="demeaned v*") # so linearly predicted v resembles computed optimum



# so A*h = (v-mean(v)), and you can also go backward, yielding the impulse response
plot(A[8:end,1:15] \ (pv))

# we can invert the A (bump) matrix to get the impulse response
plot(h,label="optimal up")
plot!(A[8:end,1:15] \ vs0[1:end-7],label="pinv response")


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