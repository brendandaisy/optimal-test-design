using DrWatson
@quickactivate "optimal-test-design"
using DEParamDistributions
using Distributions
using CSV, DataFrames
using NamedTupleTools
using Distributed, ClusterManagers
import Dates: today, format

include(srcdir("watson-tools.jl"))
include(srcdir("cond-simulations.jl"))
ENV["JULIA_WORKER_TIMEOUT"] = 120.

function inct_exper!(d, cond_sims; N=100)
    @unpack θtrue, known, obs_model, obs_params = d
    pset = Set(keys(θtrue))
    true_sim = cond_sims[pset].u
    pri_sims = cond_sims[known]

    for θᵢ ∈ setdiff(pset, known) # for each unknown var
        cond_simᵢ = cond_sims[union(Set([θᵢ]), known)] # get sims x∣known, θᵢ
        lab = "sig-"*string(θᵢ)
        d[lab] = []
        for imax ∈ eachindex(true_sim) # for each timespan
            obsp = merge(obs_params, (maxt=imax,))
            likelihood = inct_dict[obs_model]
            push!(d[lab], local_marginal_utility(true_sim, cond_simᵢ, pri_sims, obsp, likelihood; N))
        end
    end
end

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
dekwargs = (saveat=2, save_idxs=2) # observations may occur at Δt=2 intervals at comparment 2 (infectious)
known = [Set([:α]), Set([:β]), Set([:S₀]), Set{Symbol}()]
obs_model = "poisson"
# obs_params = (n=1000,)
obs_params = [(n=ntest,) for ntest ∈ [10, 100, 1000]]

cond_sims = get_cond_sims(θtrue, θprior, 2500; dekwargs...)

factors = @strdict θtrue known obs_model obs_params

vacc = length(ARGS) > 0
if vacc
    addprocs(SlurmManager(parse(Int, ARGS[1])), partition="short")
    @everywhere using DrWatson
    @everywhere begin
        @quickactivate "optimal-test-design"
        using Distributions, DEParamDistributions
    end
    fname = datadir("sims", "increasing-tspan")
    safe = true
else
    fname = "_research/tmp"
    safe = false
end

@everywhere include(srcdir("observation-dicts.jl"))

for d ∈ dict_list(factors)
    inct_exper!(d, cond_sims; N=10_000)
    tagsave("$fname/$(mysavename(d))", d; safe)
end

for i in workers()
    rmprocs(i)
end
