
module DynLocoGraphics

using DynLoco, Plots

export drawmodel, drawmodel!

const StepData = Union{NamedTuple,StepResults}

drawmodel(w::DynLoco.AbstractWalkRW2, step::StepData, args...) = drawmodel!(plot(),
    w, step, args...)

"drawmodel!([p,] w::Walk, step::StepResults) adds to an existing plot."

drawmodel!(w::DynLoco.AbstractWalkRW2, step::StepData, args...) = 
    drawmodel!(Plots.CURRENT_PLOT.nullableplot, w, step, args...)

function drawmodel!(p::Union{Plots.Plot,Plots.Subplot}, w::DynLoco.AbstractWalkRW2, 
    step::StepData, (posnx, posny)=(0,0); scalev = 1)
    Ωminus = step.Ωminus
    Ωplus = step.Ωplus
    P = step.P
    C = step.C
    δ = step.δ
    L = w.L
    α = w.α
    δangle = step.δ
    q1 = -α + δangle 
    q2 = q1 + 2*α     # fixed step lengths for rimless wheel

    feetx = [0.,  L*(sin(q2)-sin(q1))] .+ posnx
    feety = [0.,  L*(cos(q1)-cos(q2))] .+ posny
    pelvisx = [-L*sin(q1)] .+ posnx
    pelvisy = [ L*cos(q1)] .+ posny
    legsx = [feetx[1]; pelvisx; feetx[end]]
    legsy = [feety[1]; pelvisy; feety[end]]

    plot!(p, legsx, legsy, color=:black, aspectratio = 1, labels=:none, showaxis=false)
    scatter!(p, feetx, feety, markershape=:circle, color=:black, markersize = 3, 
        markerstrokecolor=nothing, label=:none)
    # plot dots
    scatter!(p, pelvisx, pelvisy, markershape=:circle, color=:black, markersize=10, 
        markerstrokecolor=nothing, label=:none)

    # plot s2s transition
    vcomminusx = -L*Ωminus.*[0, cos(q1)] .* scalev .+ pelvisx
    vcomminusy = -L*Ωminus.*[0, sin(q1)] .* scalev .+ pelvisy
    pushoffx = P.*[0, -sin(q1)] .* scalev .+ vcomminusx[end]
    pushoffy = P.*[0, cos(q1)] .* scalev .+ vcomminusy[end]
    collisionx = C.*[0, -sin(q2)] .* scalev .+ pushoffx[end]
    collisiony = C.*[0, cos(q2)] .* scalev .+ pushoffy[end]
    vcomplusx = -L*Ωplus.*[0, cos(q2)] .* scalev .+ pelvisx
    vcomplusy = -L*Ωplus.*[0, sin(q2)] .* scalev .+ pelvisy
    plot!(p, vcomminusx, vcomminusy, color=:blue, arrow=:arrow, label=:none) # vcomminus
    plot!(p, pushoffx, pushoffy, color=:green, arrow=:arrow, label=:none) # pushoff
    plot!(p, collisionx, collisiony, color=:red, arrow=:arrow, label=:none) # collision
    plot!(p, vcomplusx, vcomplusy, color=:blue, arrow=:arrow, label=:none) # vcomplus

    # foot floor
    floortilex = feetx 
    floortiley = feety 
    plot!(p, floortilex, floortiley, linecolor=:gray, label=:none)

end

end