using Distributions, DEParamDistributions
using OrdinaryDiffEq
using Plots, StatsPlots
using Random
using DataFrames
using RCall

Random.seed!(1234)

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
dekwargs = (save_idxs=2, saveat=0.1)
ntest = 100
rate = 10

pdist = SIRParamDistribution(;θprior...)
xtrue = solve(de_problem(pdist, θtrue; dekwargs...), Tsit5()).u
sim = simulate(pdist, 500; keep=true, dekwargs...)

qhi = EnsembleSummary(sim).qhigh
qlo = EnsembleSummary(sim).qlow
subsim = [sim[i].u for i ∈ sample(1:500, 50)]
curves = DataFrame(subsim, :auto)
curves.t = sim[1].t
@rput qhi qlo curves

function get_post(tmax, y=rand(joint_neg_binom(rate, xtrue[1:10:(tmax*10)] .* ntest)))
    sample_mcmc(
        y, pdist, x->joint_neg_binom(rate, x[1:10:((tmax-1)*10+1)] .* ntest); dekwargs..., 
        arraylik=false, iter=1000, ŷ=true
    )
end

function get_sig(tmax, y=rand(joint_neg_binom(rate, xtrue[1:tmax] .* ntest)))
    sims = simulate(pdist, 4000; dekwargs...)
    pri_ldists = map(x->joint_neg_binom(rate, x[1:tmax] .* ntest), sims)
    DEParamDistributions.sig(y, pri_ldists, joint_neg_binom(rate, xtrue[1:tmax] .* ntest))
end

y=rand(joint_neg_binom(rate, xtrue[1:10:end] .* ntest))
_, _c6 = get_post(6, y[1:6])
_, _c12 = get_post(12, y[1:12])

c6 = DataFrame(sample(reshape(_c6, :), 50), :auto)
c6.t = sim[1].t
c12 = DataFrame(sample(reshape(_c12, :), 50), :auto)
c12.t = sim[1].t

@rput c6 c12 y xtrue

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