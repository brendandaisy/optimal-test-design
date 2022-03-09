using DrWatson
@quickactivate "optimal-test-design"
using DEParamDistributions
using Distributions
using CSV, DataFrames
using Distributed, ClusterManagers
import Dates: today, format

include(srcdir("observation-dicts.jl"))
ENV["JULIA_WORKER_TIMEOUT"] = 120.

todaystr() = format(today(), "mm-dd")
mysavename(d) = replace(
    savename(d,  "jld2"; connector=" || ", allowedtypes=(String, Tuple, NamedTuple), ignores = (:dekwargs, :θprior)), 
    " = " => "="
)

function uoft_exper(d; N=100, M=100)
    @unpack θtrue, θprior, dekwargs, param_comb, obs_model, obs_params = d
    ret = []
    parms = keys(θprior) # param names
    vals = (θ ∈ param_comb ? getindex(θprior, θ) : getindex(θtrue, θ) for θ ∈ parms) # merge priors and fixed vals
    pdist = SIRParamDistribution(;zip(parms, vals)...) # make a SIR prob with priors and fixed vals
    precomps = all_designs_precomps(θtrue, pdist; umap=(SIG=M,), dekwargs...) # simulate SIR curves from prior
    ts = Int.(1:dekwargs.saveat:pdist.stop)
    for t ∈ ts
        util = local_utility(t, θtrue, pdist, x->obs_func(t, x; obs_params...); N, precomps, dekwargs...)
        push!(ret, util.SIG)
    end
    return ret
end

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
dekwargs = (saveat=1, save_idxs=2) # observations may occur at Δt=1 intervals at comparment 2 (infectious)
param_comb = [(:S₀, :β), (:β, :α), keys(θtrue)]
obs_model = "neg_binom"
obs_params = [(r=rate, n=ntest) for rate ∈ [1, 20] for ntest ∈ [10, 100]]

factors = @strdict θtrue θprior dekwargs param_comb obs_model obs_params

vacc = length(ARGS) > 0
if vacc
    addprocs(SlurmManager(parse(Int, ARGS[1])), partition="short")
    @everywhere using DrWatson
    @everywhere begin
        @quickactivate "optimal-test-design"
        using Distributions, DEParamDistributions
        include(srcdir("observation-dicts.jl"))
    end
    fname = datadir("sims", "single-observation")
else # make some diagnostic plots
    using OrdinaryDiffEq, Plots, StatsPlots
    pd = SIRParamDistribution(;θprior...)
    xtrue = solve(de_problem(pd, θtrue; dekwargs...), Tsit5()).u
    sim = simulate(pd, 500; keep=true, dekwargs...)
    display(plot(EnsembleSummary(sim), title="S₀~Unif$(params(pd.S₀)), β~Unif$(params(pd.β)), α~Unif$(params(pd.α))"))
    display(plot!(sim[rand(1:500, 100)], linealpha=0.15, lc=:gray))
    display(plot!(xtrue; lc=:orange))
    fname = "_research/tmp"
end

@everywhere obs_func(t, x; kw...) = single_obs_dict[$obs_model](t, x; kw...)

for d ∈ dict_list(factors)
    res = uoft_exper(d; N=40_000, M=20_000)
    d["utils"] = res
    @tagsave("$fname/$(mysavename(d))", d)
end

for i in workers()
    rmprocs(i)
end
