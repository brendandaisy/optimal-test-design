using Distributions, DEParamDistributions
using OrdinaryDiffEq
using Plots, StatsPlots
using Random
using DataFrames
using RCall

Random.seed!(1234)

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
dekwargs = (save_idxs=2, saveat=1)
ntest = 100
rate = 10

pdist = SIRParamDistribution(;θprior...)
xtrue = solve(de_problem(pdist, θtrue; dekwargs...), Tsit5()).u
sim = simulate(pdist, 500; keep=true, save_idxs=2, dense=true)

plot(EnsembleSummary(sim); lc=:purple, fc=:purple, thickness_scaling=1.2, draw_arrow=true, grids=:y, yaxis="")
for _=1:19
    plot!(sim[rand(1:500)], linealpha=0.25, lc=:gray)
end
plot!(sim[rand(1:500)], linealpha=0.25, lc=:gray)
# y = rand(joint_neg_binom(rate, xtrue .* ntest))
# plot!(sim[1].t, xtrue; lc=:orange)
# scatter!(sim[1].t, y / ntest, lab="data")

function get_post(tmax, y=rand(joint_neg_binom(rate, xtrue[1:tmax] .* ntest)))
    fit = sample_mcmc(
        y, pdist, x->joint_neg_binom(rate, x[1:tmax] .* ntest); dekwargs..., 
        arraylik=false, iter=1000
    )
    return fitchain(fit, pdist)
end

function get_sig(tmax, y=rand(joint_neg_binom(rate, xtrue[1:tmax] .* ntest)))
    sims = simulate(pdist, 4000; dekwargs...)
    pri_ldists = map(x->joint_neg_binom(rate, x[1:tmax] .* ntest), sims)
    DEParamDistributions.sig(y, pri_ldists, joint_neg_binom(rate, xtrue[1:tmax] .* ntest))
end

# y = rand(joint_neg_binom(rate, xtrue[1:12] .* ntest))
# post = get_post(12, y)
# sig = get_sig(12, y)

# logpdf(post, vcat(θtrue...)) - logpdf(joint_prior(pdist), vcat(θtrue...))

logcomp = map(1:100) do _
    logpost = logpdf(get_post(12), vcat(θtrue...))
    logpri = logpdf(joint_prior(pdist), vcat(θtrue...))
    return (;logpost, logpri)
end |> DataFrame

@rput logcomp