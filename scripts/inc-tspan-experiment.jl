using DrWatson
@quickactivate "optimal-test-design"
using Revise
using DiffEqInformationTheory
using Distributions, MonteCarloMeasurements
using CSV, DataFrames
using NamedTupleTools
using Distributed
import Dates: today, format
using Plots

include(srcdir("watson-tools.jl"))
include(srcdir("cond-simulations.jl"))
# include(srcdir("observation-dicts.jl"))

function inct_marg_exper!(d; N=100, dekwargs...)
    @unpack θtrue, known, lat_mod, obs_mod = d
    ϕs = setdiff(Set(keys(θtrue)), Set(known)) # each variable to compute δ for, i.e. each unknown var
    
    # particle vectors that are constant wrt ϕ:
    inf_pri = solve(lat_mod; dekwargs..., saveat=1).u # TODO produce or load all of these solves
    inf_true = solve(lat_mod, θtrue; dekwargs..., saveat=1).u
    y = Particles(40_000, observe_dist(obs_mod; observe_params(obs_mod, inf_true)...))

    for θᵢ in ϕs # for each δ to compute
        lab = "md-"*string(θᵢ)
        d[lab] = []
        θcond = NamedTupleTools.select(θtrue, (known..., θᵢ))
        inf_cond = solve(lat_mod, θcond; dekwargs..., saveat=1).u # simulate x∣known, θᵢ
        
        for imax in 2:length(inf_true) # for each timespan (0, t)
            ts_idx = 1:imax
            ynew = bootstrap(y[ts_idx], N)
            ic_new = inf_cond[ts_idx]
            ip_new = inf_pri[ts_idx]
            push!(d[lab], marginal_divergence(ynew, ic_new, ip_new, obs_mod; M))
        end
    end
end

subfolder = ""
θtrue = (S₀=0.6f0, β=1.25f0, α=0.2f0)
θprior = (S₀=Uniform(0.1f0, 0.9f0), β=Uniform(0.3f0, 3f0), α=Uniform(0.05f0, 0.3f0))
lat_mod = SIRModel{Float32}(
    S₀=Particles(10_000, θprior.S₀), 
    β=Particles(10_000, θprior.β), 
    α=Particles(10_000, θprior.α)
)
dekwargs = (save_idxs=2,) # observations may occur at Δt=1 intervals at comparment 2 (infectious)
# known = [(), (:α,)]
known = ()
obs_mod = PoissonTests(1000)

# factors = @strdict θtrue known obs_model
factors = @strdict θtrue known lat_mod obs_mod

vacc = Threads.nthreads() > 4
if vacc
    fname = datadir("sims", "increasing-tspan-marg", subfolder)
    safe = true
else
    fname = "_research/tmp/"*subfolder
    safe = false
end

for d ∈ dict_list(factors)
    inct_marg_exper!(d; N=1500, M=2500, dekwargs...)
    print(d)
    # tagsave("$fname/$(mysavename(d))", d; safe)
end
