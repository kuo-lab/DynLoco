using Plots, Statistics, LsqFit
plotlyjs() # Use this backend to preserve fonts on export to SVG or PDF
default(grid=false, fontfamily="Helvetica") # no grid on plots

# Orendurff's data 
# Orendurff, M. S. How humans walk: Bout duration, steps per bout, and rest duration
# JRRD 45 (7): 1077-1090, 2008

# Here's the data from Fig. 3
x = 4:2:76
y = [17, 9.9, 7.3, 5.7, 4.9, 4.2, 3.7, 3.2, 2.9, 2.5, 2.5, 2.1, 1.8, 1.8, 1.5, 1.3, 1.2, 1.0, 1.0, 1.1, 0.9, 0.8, 0.7, 0.7, 0.7, 0.7, 0.7, 0.6, 0.6, 0.5, 0.4, 0.6, 0.5, 0.4,0.5, 0.4, 0.4]

plot(x,y)
plot!(x, 100 .* x.^-1.25 .+ 0.2) # pretty good but a bit low around 20

# Also tried Curvefit and found power fit is not bad
# Found LsqFit to work well:

@. model(x, p) = p[1]*x^p[2]+p[3]
p0 = [125., -1.3, 0.1]
fit = LsqFit.curve_fit(model, x, y, p0; autodiff=:forwarddiff)
#77.2390472692003
#-1.0978647695791948
#-0.24527437749659545

plot(x, model(x, fit.param), xlabel = "Steps", ylabel="% Total") # this is pretty darn good
plot!(x,y) # if you want to include the data
p = plot()

## For saving: Frequency of walking bouts vs # of steps, with curve fit
pyplot()
bar(x,y)
plot!(4:0.5:76, x->model(x, fit.param), xlabel = "Steps", ylabel="% Total")
savefig("orendurff.svg")

sse = sum(abs2, model(x, fit.param) .- y)
sst = sum(abs2, y .- mean(y))
r2 = 1 - sse / sst # 0.9964