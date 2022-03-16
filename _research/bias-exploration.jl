using DrWatson
@quickactivate "optimal-test-design"
using DEParamDistributions
using OrdinaryDiffEq
using Distributions
using Plots

## bias params are fixed: (extending to random will be good bit of refactoring)

function biussy(x::Vector, ntest, g=1/100, s=2/3, popsize=10_000)
    ret = []
    for t=eachindex(x)
        n1 = min(x[t] * popsize * s, ntest)
        push!(ret, n1*x[t]/(x[t]+g) + (ntest-n1)*x[t])
    end
    ret
end

biussy_chiu(x, ntest, b; popsize=10_000) = x * ntest * (popsize/ntest)^b

pd = SIRParamDistribution(S₀=0.6, β=1.25, α=0.2)
xtrue = solve(de_problem(pd; saveat=1, save_idxs=2), Tsit5()).u
y = rand(joint_poisson(biussy_chiu(xtrue, 1000, 0.2)))

plot(xtrue)
scatter!(y / 1000)

asim = simulate(SIRParamDistribution(S₀=Beta(3, 1), β=Uniform(0.4, 1.5), α=0.2), 500; keep=true, saveat=1, save_idxs=2)
plot(EnsembleSummary(asim))

ŷ = prior_predict(asim, x->joint_poisson(biussy_chiu(x, 1000, 0.2)))
μ = [mean(getindex.(ŷ, i)) for i=eachindex(xtrue)]
low = [quantile(getindex.(ŷ, i), 0.05) for i=eachindex(xtrue)]
hi = [quantile(getindex.(ŷ, i), 0.95) for i=eachindex(xtrue)]
plot(biussy_chiu(xtrue, 1000, 0.2), lc=:orange)
plot!(μ)
plot!(low, lc=:green, ls=:dash)
plot!(hi, lc=:green, ls=:dash)

## find some reasonable b for n=10, 100, 1000 total tests

# b1ussy(ntest, b=0.5; popsize=10_000) = ntest*((popsize-ntest)/popsize)^b
# b2ussy(x, ntest, b=0.5; g=1/100) = ntest^b * x/(x+g) + (ntest-ntest^b) * x
# b4ussy(x, ntest, g=1/100, s=2/3, popsize=10_000) = ntest^b * x/(x+g) + (ntest-ntest^b) * x

# nt = 10:1000
# plot(nt, nt; ls=:dash, lc=:gray)

# plot!(nt, b1ussy.(nt, 1.))
# plot!(nt, b1ussy.(nt, 10.))
# plot!(nt, b1ussy.(nt, 100.))

# plot!(nt, b2ussy.(nt, 0.1))
# plot!(nt, b2ussy.(nt, 0.5))
# plot!(nt, b2ussy.(nt, 0.9))

# let f=(n, b)->log(2-b, n)
#     plot(nt, nt; ls=:dash, lc=:gray)
#     plot!(nt, f.(nt, 0.1))
#     plot!(nt, f.(nt, 0.5))
#     plot!(nt, f.(nt, 0.9))
# end

let x=0.5, f=b3ussy
    plot(nt, nt; ls=:dash, lc=:gray)
    plot!(nt, f.(x, nt, 0.1))
    plot!(nt, f.(x, nt, 0.3))
    plot!(nt, f.(x, nt, 0.9))
end

