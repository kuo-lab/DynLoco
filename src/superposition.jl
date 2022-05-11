using DynLoco
using Plots, Statistics
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

wstar4s = findgait(WalkRW2l(α=0.4,safety=true), target=:speed=>0.45, varying=:P)
nsteps = 15
δs = zeros(nsteps); δs[Int((nsteps+1)/2)] = 0.05
nominalmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δs)
plotvees(nominalmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)



nsteps = 100
δs = rand(nsteps)./10
coefs = hcat(0:nsteps-1,ones(nsteps))\δs # slope and offset
δs .= δs .- hcat(0:nsteps-1,ones(nsteps))*coefs # detrend the bumps
nominalmsr=optwalk(wstar4s, nsteps, boundarywork=false, δs=δs)
plotvees(nominalmsr,speedtype=:midstance,usespline=false,boundaryvels=(wstar4s.vm,wstar4s.vm),tchange=0)

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
A = Toeplitz(δs, [δs[1];zeros(nsteps-1)])
vs = convert(Vector{Float64}, nominalmsr.steps.vm)
h = A[:,1:8] \ vs

plot(h)