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
# include(srcdir("observation-dicts.jl"))

function inct_marg_exper!(d, cond_sims; N=100, M=100)
    @unpack θtrue, known, obs_model = d
    pset = Set(keys(θtrue))
    true_sim = cond_sims[pset].u
    pri_sims = cond_sims[known]

    for θᵢ ∈ setdiff(pset, known) # for each unknown var
        cond_simᵢ = cond_sims[union(Set([θᵢ]), known)] # get sims x∣known, θᵢ
        lab = "sig-"*string(θᵢ)
        d[lab] = []
        for imax ∈ 2:length(true_sim) # for each timespan
            likelihood = (sim, m)->obs_tspan(sim, m, imax)
            push!(
                d[lab], 
                local_marginal_utility(true_sim, cond_simᵢ, pri_sims, likelihood; N, M, obs_mod=obs_model)
            )
        end
    end
end

subfolder = ""
θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3.), α=Uniform(0.05, 0.3))
dekwargs = (saveat=1, save_idxs=2) # observations may occur at Δt=1 intervals at comparment 2 (infectious)
known = Set{Symbol}()
# obs_model = "neg_binom"
# obs_params = [(r=rate, n=ntest) for rate ∈ [1, 10] for ntest ∈ [10, 100, 1000]]
obs_model = PoissonTests(1000.)

cond_sims = get_cond_sims(θtrue, θprior, 40_000; dekwargs...)

factors = @strdict θtrue known obs_model

vacc = Threads.nthreads() > 4
if vacc
    fname = datadir("sims", "increasing-tspan-marg", subfolder)
    safe = true
else
    fname = "_research/tmp/"*subfolder
    safe = false
end

for d ∈ dict_list(factors)
    inct_marg_exper!(d, cond_sims; N=2500, M=2500)
    tagsave("$fname/$(mysavename(d))", d; safe)
end
