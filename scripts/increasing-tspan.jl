using DrWatson
@quickactivate "optimal-test-design"
using DEParamDistributions
using Distributions
using CSV, DataFrames
using Distributed, ClusterManagers
import Dates: today, format

include(srcdir("watson-tools.jl"))
ENV["JULIA_WORKER_TIMEOUT"] = 120.

function inct_exper(d; N=100, M=100)
    @unpack θtrue, θprior, dekwargs, param_comb, obs_model, obs_params = d
    ret = []
    parms = keys(θprior) # param names
    vals = (θ ∈ param_comb ? getindex(θprior, θ) : getindex(θtrue, θ) for θ ∈ parms) # merge priors and fixed vals
    pdist = SIRParamDistribution(;zip(parms, vals)...) # make a SIR prob with priors and fixed vals
    precomps = all_designs_precomps(θtrue, pdist; umap=(SIG=M,), dekwargs...) # simulate SIR curves from prior
    imaxs = 1:floor(Int, pdist.stop / dekwargs.saveat) # these are ~indexes~ to x(t) where t=i*saveat
    for imax ∈ imaxs
        # note first arg may never actually be needed? Consider rethinking
        util = local_utility(-1, θtrue, pdist, x->obs_func(imax, x, obs_model; obs_params...); N, precomps, dekwargs...)
        push!(ret, util.SIG)
    end
    return ret
end

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
dekwargs = (saveat=2, save_idxs=2) # observations may occur at Δt=2 intervals at comparment 2 (infectious)
param_comb = [(:S₀, :β), (:β, :α), keys(θtrue)]
obs_model = "neg_binom"
# obs_params = (n=1000,)
obs_params = [(r=rate, n=1000) for rate ∈ [1, 10]]

factors = @strdict θtrue θprior dekwargs param_comb obs_model obs_params

vacc = length(ARGS) > 0
if vacc
    addprocs(SlurmManager(parse(Int, ARGS[1])), partition="short")
    @everywhere using DrWatson
    @everywhere begin
        @quickactivate "optimal-test-design"
        using Distributions, DEParamDistributions
    end
    fname = datadir("sims", "increasing-tspan")
else
    fname = "_research/tmp"
end

@everywhere include(srcdir("observation-dicts.jl"))
@everywhere obs_func(maxt, x, mod; kw...) = inct_dict[mod](maxt, x; kw...)

for d ∈ dict_list(factors)
    res = inct_exper(d; N=8000, M=2000)
    # res = inct_exper(d; N=200, M=200)
    d["utils"] = res
    @tagsave("$fname/$(mysavename(d))", d)
end

for i in workers()
    rmprocs(i)
end
