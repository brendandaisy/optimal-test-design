using DrWatson
@quickactivate "optimal-test-design"
using Revise
using DiffEqInformationTheory
using Distributions, MonteCarloMeasurements
using NamedTupleTools
using Distributed

include(srcdir("watson-tools.jl"))

function inct_single_param_exper!(d; N=100, dekwargs...)
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
            push!(d[lab], marginal_divergence(ynew, ic_new, ip_new, obs_mod))
        end
    end
end

addprocs(8)

subfolder = ""

lat_mod = SIRModel{Float32}(
    S₀=Particles(60_000, θprior.S₀), 
    β=Particles(60_000, θprior.β), 
    α=Particles(60_000, θprior.α)
)
dekwargs = (save_idxs=2,) # observations may occur at Δt=1 intervals at comparment 2 (infectious)
# known = [(), (:α,)]
known = ()
obs_mod = PoissonTests(1000)

# factors = @strdict θtrue known obs_model
factors = @strdict θtrue known lat_mod obs_mod

if nprocs() > 1
    @everywhere using DrWatson
    @everywhere begin
        @quickactivate "optimal-test-design"
        using DiffEqInformationTheory
        using Distributions, MonteCarloMeasurements
    end
    fdir = datadir("sims", "increasing-tspan-marg", subfolder)
    safe = true
else
    fdir = "_research/tmp/"*subfolder
    safe = false
end

for d ∈ dict_list(factors)
    fname = mysavename(d)
    inct_single_param_exper!(d; N=3000, dekwargs...)
    tagsave("$fdir/$fname", d; safe)
end
