using DrWatson
using DEParamDistributions
using Distributions
using CSV, DataFrames
using Distributed, ClusterManagers
import Dates: today, format

todaystr() = format(today(), "mm-dd")

function uoft_exper(d; N=100, M=100, ts=1:30)
    @unpack θtrue, θprior, obs_args, param_comb, obs_type, obs_param = d
    ret = []
    parms = keys(θprior)
    vals = (θ ∈ param_comb ? getindex(θprior, θ) : getindex(θtrue, θ) for θ ∈ parms)
    pdist = SIRParamDistribution(;zip(parms, vals)...)
    precomps = all_designs_precomps(θtrue, pdist; umap=(SIG=M,), obs_args...)
    for t=ts
        util = local_utility(t, θtrue, pdist, x->single_obs(t, x; obs_param...); N, precomps, obs_args...)
        push!(ret, (t=t, free=param_comb, obs_param..., util...))
    end
    return ret
end

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
obs_args = (saveat=1, save_idxs=2) # observations may occur at Δt=1 intervals at comparment 2 (infectious)
param_comb = [(:S₀, :β), (:β, :α), keys(θtrue)]
obs_type = "normal"
obs_param = [(n=1000, σ=10.,), (n=1000, σ=40.,)]

factors = Dict{Symbol, Any}()
@pack! factors = θtrue, θprior, obs_args, param_comb, obs_type, obs_param

vacc = length(ARGS) > 0
if vacc
    addprocs(SlurmManager(parse(Int, ARGS[1])), partition="short")
    @everywhere using Distributions
    @everywhere using DEParamDistributions
    fname = datadir("sim", "single-observation")
else # make some diagnostic plots
    using OrdinaryDiffEq, Plots, StatsPlots
    pd = SIRParamDistribution(;θprior...)
    xtrue = solve(de_problem(pd, θtrue; obs_args...), Tsit5()).u
    sim = simulate(pd, 500; keep=true, obs_args...)
    display(plot(EnsembleSummary(sim), title="S₀~Unif$(params(pd.S₀)), β~Unif$(params(pd.β)), α~Unif$(params(pd.α))"))
    display(plot!(sim[rand(1:500, 100)], linealpha=0.15, lc=:gray))
    display(plot!(xtrue; lc=:orange))
    fname = "_research/tmp/"
end

@everywhere single_obs(t, x; σ, n) = Normal(n * x[t], σ) # now make noise ∝ inf size

res = []
for d ∈ dict_list(factors)
    append!(res, uoft_exper(d; N=50_000, M=20_000))
end
CSV.write(joinpath(fname, "$obs_type-$(todaystr()).csv"), DataFrame(res))

for i in workers()                                                             
    rmprocs(i)                                                                 
end