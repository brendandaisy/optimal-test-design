using Distributions, DEParamDistributions
using OrdinaryDiffEq
using Plots, StatsPlots
using Random

Random.seed!(1234)

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
dekwargs = (save_idxs=2, saveat=1)
ntest = 100
rate = 10

pdist = SIRParamDistribution(;θprior...)
xtrue = solve(de_problem(pdist, θtrue; dekwargs...), Tsit5()).u
sim = simulate(pdist, 500; keep=true, dekwargs...)
y = rand(joint_neg_binom(rate, xtrue .* ntest))
plot(EnsembleSummary(sim), title="S₀~Unif$(params(pdist.S₀)), β~Unif$(params(pdist.β)), α~Unif$(params(pdist.α))")
for _=1:100
    plot!(sim[rand(1:500)], linealpha=0.15, lc=:gray)
end
plot!(sim[1].t, xtrue; lc=:orange)
scatter!(sim[1].t, y / ntest, lab="data")

## posterior change from prior at several values

fit6, y6 = sample_mcmc(y[1:6], pdist, x->joint_neg_binom(rate, x[1:6] .* ntest); dekwargs..., arraylik=false, ŷ=true)
fit12, y12 = sample_mcmc(
    y[1:12], pdist, x->joint_neg_binom(rate, x[1:12] .* ntest); dekwargs..., 
    arraylik=false, iter=1000, ŷ=true
)

plot(fit12)

p1 = plot(sample(reshape(y6, :), 250), lc=:gray, linealpha=0.3, lab="")
scatter!(p1, sim[1].t[1:6], y[1:6] / ntest, lab="Observations up to t=6")
plot!(p1, sim[1].t, xtrue; lc=:orange, lab="True inf")

p2 = plot(sample(reshape(y12, :), 250), lc=:gray, linealpha=0.3, lab="")
scatter!(p2, sim[1].t[1:12], y[1:12] / ntest, lab="Observations up to t=12")
plot!(p2, sim[1].t, xtrue; lc=:orange, lab="True inf")

plot(p1, p2, layout=(1, 2))

histogram(fit12; normalize=true)
plot!(θprior.S₀, subplot=1, lc=:green, lab="S₀~Unif$(params(pdist.S₀))", lw=2.5)
plot!(θprior.β, subplot=2, lc=:green, lab="β~Unif$(params(pdist.β))", lw=2.5)
plot!(θprior.α, subplot=3, lc=:green, lab="α~Unif$(params(pdist.α))", lw=2.5)
vline!([θtrue.S₀], subplot=1, lc=:orange, lab="True S₀", ls=:dash, lw=2.5)
vline!([θtrue.β], subplot=2, lc=:orange, lab="True β", ls=:dash, lw=2.5)
vline!([θtrue.α], subplot=3, lc=:orange, lab="True α", ls=:dash, lw=2.5)
