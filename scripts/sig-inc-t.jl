using DEParamDistributions
using Distributions
using CSV, DataFrames
using Distributed, ClusterManagers
import Dates: today

function inc_tspan_exper(σ, θtrue, θprior, param_combs; N=100, ts=1:30, obs_args)
    ret = []
    parms = keys(θprior)
    for comb ∈ param_combs
        vals = (θ ∈ comb ? getindex(θprior, θ) : getindex(θtrue, θ) for θ ∈ parms)
        pdist = SIRParamDistribution(;zip(parms, vals)...)
        precomps = all_designs_precomps(θtrue, pdist; umap=(SIG=20_000,), obs_args...)
        for t=ts
            util = local_utility(t, θtrue, pdist, x->single_obs(t, x, σ); N, precomps, obs_args...)
            push!(ret, (t=t, free=comb, sd=σ, util...))
        end
    end
    return ret
end

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
obs_args = (saveat=1, save_idxs=2) # observations may occur at Δt=1 intervals at comparment 2 (infectious)
# param_combs = [(:S₀,), (:β,), (:α,), (:S₀, :β), (:β, :α), keys(θtrue)]
param_combs = [(:S₀, :β), (:β, :α), keys(θtrue)]

remote = length(ARGS) > 0
if remote
    @info "*************Beginning remote process*************"
    addprocs(SlurmManager(parse(Int, ARGS[1])), partition="short")
    @everywhere using Distributions
    @everywhere using DEParamDistributions
else # make some diagnostic plots
    using OrdinaryDiffEq, Plots, StatsPlots
    pd = SIRParamDistribution(;θprior...)
    xtrue = solve(de_problem(pd, θtrue; obs_args...), Tsit5()).u
    sim = simulate(pd, 500; keep=true, obs_args...)
    display(plot(EnsembleSummary(sim), title="S₀~Unif$(params(pd.S₀)), β~Unif$(params(pd.β)), α~Unif$(params(pd.α))"))
    display(plot!(sim[rand(1:500, 100)], linealpha=0.15, lc=:gray))
    display(plot!(xtrue; lc=:orange))
end

@everywhere single_obs(t, x, σ) = Normal(100 * x[t], σ * x[t]) # now make noise ∝ inf size
# @everywhere single_obs(t, numtest, x) = Poisson(numtest * x[t])
# @everywhere single_obs(d, x) = single_obs(d[1], d[2], x)

ex1 = uoft_exper(10., θtrue, θprior, param_combs; N=20_000, obs_args)
ex2 = uoft_exper(40., θtrue, θprior, param_combs; N=20_000, obs_args)
CSV.write("noise-of-inf-$(today()).csv", DataFrame(append!(ex1, ex2)))


# ex1 = uoft_exper(3, 1, pdist, precomps; N=20_000)
# ex2 = uoft_exper(3, 10, pdist, precomps; N=20_000)
# ex3 = uoft_exper(3, 100, pdist, precomps; N=20_000)
# CSV.write("single-obs-uoft.csv", DataFrame(append!(ex1, ex2, ex3)))