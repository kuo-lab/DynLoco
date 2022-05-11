## Normalizationeffects to analyze data from short walking experiment

using Statistics, HypothesisTests, MAT, Plots

# From Elizabeth slack, the standard deviations of peak speeds and durations
# as a function of distance. This shows that normalization causes the variability
# to be reduced.

#Peak speed s.d. non-normalized: 0.1630    0.1351    0.1688    0.1399    0.1221    0.1414    0.1304    0.1252    0.1244    0.1377
#Peak speed s.d. normalized:         0.1120    0.0310    0.1024    0.1000    0.0314    0.0821    0.0434    0.0341    0.0651    0.0506
#Walk duration non-normalized: 0.4055    0.3026    0.3180    0.3345    0.3584    0.3733    0.5925    0.5563    0.8151    0.7486
#Walk duration normalized:         0.3779    0.2939    0.2931    0.2722    0.3154    0.3450    0

peakspeeds = [0.1630    0.1351    0.1688    0.1399    0.1221    0.1414    0.1304    0.1252    0.1244    0.1377;
              0.1120    0.0310    0.1024    0.1000    0.0314    0.0821    0.0434    0.0341    0.0651    0.0506]
durations = [0.4055    0.3026    0.3180    0.3345    0.3584    0.3733    0.5925    0.5563    0.8151    0.7486;
             0.3779    0.2939    0.2931    0.2722    0.3154    0.3450    0.4773    0.4441    0.6245    0.3084];

             # mean reduction in speed sd (s)
mean(peakspeeds[2,:].-peakspeeds[1,:]) # -0.07358999999999999
std(peakspeeds[2,:].-peakspeeds[1,:]) # 0.021076550951234882
mean((peakspeeds[2,:].-peakspeeds[1,:]) ./ peakspeeds[1,:]) # -0.5428240111838153, 54% reduction

# mean reduction in durations (s) 
mean(durations[2,:].-durations[1,:]) # -0.10529999999999999
std(durations[2,:].-durations[1,:]) # 0.13041473161503733
mean((durations[2,:].-durations[1,:])./durations[1,:]) # -0.17751404274419505, 18% reduction

pvalue(OneSampleTTest(peakspeeds[2,:], peakspeeds[1,:])) # p = 1.56e-6
pvalue(OneSampleTTest(durations[2,:], durations[1,:])) # p = 0.031


using Plots
plot(peakspeeds')
plot(durations')

# Read in Elizabeth's mat file

peakspeeds = [vec(std(vars["peak_speed"],dims=(1,3))), vec(std(vars["peak_speed_norm"],dims=(2)))]
# mean reduction in speed sd (s)
@show mean(peakspeeds[2] .- peakspeeds[1]) # in absolute amount
@show mean(peakspeeds[2] .- peakspeeds[1])/mean(peakspeeds[1]) # about half reduction

# mean reduction in durations
durations = [vec(std(vars["walk_duration"],dims=(1,3))), vec(std(vars["walk_duration_norm"],dims=(2)))]
@show mean(durations[2] .- durations[1]) # in absolute amount
@show mean(durations[2] .- durations[1])/mean(durations[1]) # about half reduction


# Here's another way to read the files, including getting variable names:
#==
file = matopen("short_walks.mat")
varnames = names(file)
peakspeed = read(file, "peak_speed")
walk_duration = read(file, "walk_duration")
close(file)
==#

using MAT
vars = matread("short_walks.mat")
# variables are peak_speed, walk_duration, and normalized versions
# where they are 10x10x4, subject x distance x trial

# peakspeed is 10 x 10 x 4, subject x cond x 4, which should be converted into a 10 x 40 cond x everything
peakspeed=permutedims(reshape(permutedims(vars["peak_speed"],(1,3,2)), (40,10)), (2,1)) # to yield 10x40
plot(peakspeed, linestyle=:dot)
plot!(mean(peakspeed, dims=2), linewidth=2, linecolor=:black, yerror=std(peakspeed, dims=2))

# peakspeednorm is stored differently, as 10x40 to begin with; needs to be de-normalized
peakspeednorm = vars["peak_speed_norm"]*1.5157 # to regular m/s
plot!(peakspeednorm,linestyle=:dot)
plot!(mean(peakspeednorm,dims=2), linewidth=2, linecolor=:black, yerror=std(peakspeednorm,dims=2))

# compare the peak speeds against normalized
@show mean(std(peakspeednorm,dims=2) .- std(peakspeed,dims=2)) # a reduction of 0.043 m/s (not much)
@show (mean(std(peakspeednorm,dims=2) .- std(peakspeed,dims=2)))/mean(std(peakspeed,dims=2)) # a reduction of 0.043 m/s (not much)
# or about 31.1% so not terrible
@show pvalue(OneSampleTTest(vec(std(peakspeed,dims=2) .- std(peakspeednorm,dims=2)))) # p = 0.0011468

# duration is 10 x 10 x 4, subject x cond x 4, which should be converted into a 10 x 40 cond x everything
duration=permutedims(reshape(permutedims(vars["walk_duration"],(1,3,2)), (40,10)), (2,1)) # to yield 10x40
plot(duration, linestyle=:dot)
plot!(mean(duration, dims=2), linewidth=2, linecolor=:black, yerror=std(duration, dims=2))

# peakspeednorm is stored differently, as 10x40 to begin with; needs to be de-normalized
durationnorm = vars["walk_duration_norm"]*1.4696/0.1389 # to regular s
plot!(durationnorm,linestyle=:dash)
plot!(mean(durationnorm,dims=2), linewidth=2, linecolor=:red, yerror=std(durationnorm,dims=2),legend=:none)

# compare the peak speeds against normalized
@show mean(std(duration,dims=2) .- std(durationnorm,dims=2)) # a reduction of 0.1055 s (not much)
@show (mean(std(duration,dims=2) .- std(durationnorm,dims=2)))/mean(std(duration,dims=2)) # a reduction of 0.043 m/s (not much)
# or about 24.3% so not terrible
@show pvalue(OneSampleTTest(vec(std(duration,dims=2) .- std(durationnorm,dims=2)))) # p = 0.0011468


# I guess I should be able to normalize myself, taking an average of the four trials of each person
# except taking each person's peak speed from their longest bout
allpeaks = vars["peak_speed"]
# to make things easier, let's average each person's four trapezoidresults
fewerpeaks = dropdims(mean(allpeaks, dims=3),dims=3) # 10x10 subject x cond
# personpeaks is the maximum for each subject across their four trials
personpeaks = maximum(allpeaks[:,10,:],dims=2)
normpeaks = [fewerpeaks[i,j]./personpeaks[i] for i in 1:10, j in 1:10]*mean(personpeaks)
plot(fewerpeaks')
plot!(mean(fewerpeaks,dims=1)', yerror=std(fewerpeaks,dims=1)')
plot!(normpeaks', legend=:none, linestyle=:dot)
plot!(mean(normpeaks,dims=1)', yerror=std(normpeaks,dims=1)', linestyle=:dot)

std(normpeaks,dims=1) .- std(fewerpeaks,dims=1) # this definitely also shows a reduction
@show mean(std(normpeaks,dims=1) .- std(fewerpeaks,dims=1)) # a reduction of 0.0587929 m/s
@show (mean(std(normpeaks,dims=1) .- std(fewerpeaks,dims=1)))/mean(std(fewerpeaks,dims=1)) # or 47.456789% 
@show pvalue(OneSampleTTest(vec(std(normpeaks,dims=1) .- std(fewerpeaks,dims=1)))) # p = 0.0005814

# same thing with durations
alldurations = vars["walk_duration"]
# to make things easier, let's average each person's four trials
personpeakdurations = maximum(alldurations[:,10,:],dims=2)
fewerdurations = dropdims(mean(alldurations, dims=3),dims=3) # 10x10 subject x cond
normdurations = [fewerdurations[i,j]./personpeakdurations[i] for i in 1:10, j in 1:10]*mean(personpeakdurations)
plot(fewerdurations')
plot!(mean(fewerdurations,dims=1)', yerror=std(fewerdurations,dims=1)')
plot!(normdurations', legend=:none, linestyle=:dot)
plot!(mean(normdurations,dims=1)', yerror=std(normdurations,dims=1)', linestyle=:dot)

std(normdurations,dims=1) .- std(fewerdurations,dims=1) # this definitely also shows a reduction
@show mean(std(normdurations,dims=1) .- std(fewerdurations,dims=1)) # a reduction of 0.1624 sec
@show (mean(std(normdurations,dims=1) .- std(fewerdurations,dims=1)))/mean(std(fewerdurations,dims=1)) # or 43.513% 
@show pvalue(OneSampleTTest(vec(std(normdurations,dims=1) .- std(fewerdurations,dims=1)))) # p = 0.025843

## My own normalization, taking each person's longest bout maximum and
# applying all subjects, all conditions
allpeaks = permutedims(vars["peak_speed"],(1,3,2)) # subject x trial x cond
# personpeaks is the maximum for each subject across their four trials
personpeaks = vec(maximum(allpeaks[:,:,10],dims=2)) # 10 subjects
normpeaks = [allpeaks[i,j,k]./personpeaks[i] for i in 1:10, j in 1:4, k in 1:10]*mean(personpeaks)
allpeakstogether = reshape(allpeaks,(40,10))
plot(allpeakstogether')
plot!(mean(allpeakstogether,dims=1)', yerror=std(allpeakstogether,dims=1)')
plot!(reshape(normpeaks,(40,10))', legend=:none, linestyle=:dot)
plot!(vec(mean(normpeaks,dims=(1,2))), yerror=vec(std(normpeaks,dims=(1,2))), linestyle=:dot)

vec(std(normpeaks,dims=(1,2)) .- std(allpeakstogether,dims=(1)))' # this definitely also shows a reduction
@show mean(std(normpeaks,dims=(1,2)) .- std(allpeakstogether,dims=(1))) # a reduction of 0.044634 m/s
@show mean(std(normpeaks,dims=(1,2)) .- std(allpeakstogether,dims=(1)))/mean(std(allpeakstogether,dims=(1))) # or 32.0876% 
@show pvalue(OneSampleTTest(vec(std(normpeaks,dims=(1,2)) .- std(allpeakstogether,dims=(1))))) # p = 1.202665e-18

# same thing with durations
alldurations = permutedims(vars["walk_duration"],(1,3,2)) # subject x trial x cond
personpeakdurations = vec(maximum(alldurations[:,:,10],dims=2)) # 10 subjects
normdurations = [alldurations[i,j,k]./personpeakdurations[i] for i in 1:10, j in 1:4, k in 1:10]*mean(personpeakdurations)
alldurationstogether = reshape(alldurations,(40,10))
plot(alldurationstogether')
plot!(mean(alldurationstogether,dims=1)', yerror=std(alldurationstogether,dims=1)')
plot!(reshape(normdurations,(40,10))', legend=:none, linestyle=:dot)
plot!(vec(mean(normdurations,dims=(1,2))), yerror=vec(std(normdurations,dims=(1,2))), linestyle=:dot)

@show vec(std(normdurations,dims=(1,2)) .- std(alldurationstogether,dims=(1)))' # this definitely also shows a reduction
@show mean(std(normdurations,dims=(1,2)) .- std(alldurationstogether,dims=(1))) # a reduction of 0.121772 sec
@show mean(std(normdurations,dims=(1,2)) .- std(alldurationstogether,dims=(1)))/mean(std(alldurationstogether,dims=(1))) # or 28.0087% 
@show pvalue(OneSampleTTest(vec(std(normdurations,dims=(1,2)) .- std(alldurationstogether,dims=(1))))) # p = 4.014572e-10
