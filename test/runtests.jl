using DynLocoT
using Test

@testset "DynLocoT.jl" begin
    wrw = Walk(P=0.15, vm = 0.38)
    @test typeof(wrw) == Walk # basic constructor for Walk

    wstar = findgait(wrw, target=:speed=>0.4, varying=:P) # constraints olver
    @test typeof(wstar) == Walk

    @test typeof(onestep(wrw)) <: NamedTuple # onestep returns a named tuple

    @test onestep(wstar).vm ≈ wstar.vm # limit cycle
    @test islimitcycle(wstar) # a built-in limit cycle tester


    @test islimitcycle(findgait(wrw, P=0.2, target=:speed=>0.4, varying=:P)) # vary push-off
    @test islimitcycle(findgait(wrw; vm=0.35, P=0.15, α=0.35, :γ=>0.15)) # vary gravity to get limit cycle
    @test islimitcycle(findgait(wrw, target=:speed=>0.45, varying=:γ, P=0)) # gravity only, no push-off
    @test islimitcycle(findgait(wrw, target=:vm=>0.4, varying=:P)) # use vm as target
    @test islimitcycle(findgait(wrw, α=0.32, target=(:speed=>0.4,), varying=(:P,), vm=0.25, P=0.2)) # use tuple of targets

    # # test safe step with too little momentum
    @test typeof(onestep(wrw, P=0., vm=0.1, safety=true)) <: NamedTuple # should fail if safety=false
    @test_throws DomainError onestep(wrw, P=0., vm=0.1, safety=false)

    @test typeof(multistep(wstar, [0.15 0.17 0.18])) <: MultiStepResults
end
