using DrWatson
@quickactivate "optimal-test-design"
using DEParamDistributions
using Distributions
using CSV, DataFrames
using NamedTupleTools
using Distributed
import Dates: today, format

include(srcdir("watson-tools.jl"))
include(srcdir("cond-simulations.jl"))

function disp_vs_test_exper!(d, cond_sims; N, M)
    @unpack θtrue, obs_model, obs_param_grid = d
    true_sim = cond_sims[Set(keys(θtrue))].u
    pri_sims = cond_sims[Set{Symbol}()] # all vars are free
    obs_params = merge(obs_param_grid, (maxt=length(true_sim),))
    likelihood = inct_dict[obs_model]
    d["SIG"] = local_utility(true_sim, pri_sims, likelihood; N, M, obs_params)
end

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
dekwargs = (saveat=5, save_idxs=2)
obs_model = "neg_binom"
obs_param_grid = [(r=100, n=ntest) for ntest ∈ 100:100:5000]

cond_sims = get_cond_sims(θtrue, θprior, 40_000; dekwargs...)

factors = @strdict θtrue obs_model obs_param_grid

vacc = Threads.nthreads() > 4
if vacc
    fname = datadir("sims", "dispersion-vs-tests")
    savef = d->tagsave("$fname/$(mysavename(d))", d; safe=true)
else
    fname = "_research/tmp"
    savef = d->wsave("$fname/$(mysavename(d))", d)
end

include(srcdir("observation-dicts.jl"))

for d ∈ dict_list(factors)
    disp_vs_test_exper!(d, cond_sims; N=3000, M=3000)
    savef(d)
end