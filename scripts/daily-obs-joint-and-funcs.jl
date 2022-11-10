using DrWatson
@quickactivate "optimal-test-design"
using Revise
using DiffEqInformationTheory
using Distributions, MonteCarloMeasurements
using NamedTupleTools
using Distributed

using NLsolve

include(srcdir("watson-tools.jl"))

function prob_αS₀_cond_outbreak_size(α, S₀, size; prior, I₀=0.01)
    if size > S₀+I₀ # impossible since more than available individuals
        return 0
    end
    J = abs(d_inv_outbreak_size(α, S₀, size))
    pdf(prior.α, α) * pdf(prior.S₀, S₀) * J * pdf(prior.β, inv_outbreak_size(α, S₀, size))
end

function prob_αβ_cond_R(α, β; R, prior)
    pdf(prior.α, α) * pdf(prior.β, β) * α/β * pdf(prior.S₀, α*R/β)
end

function sir_cond_R(prior, f, R, nsamples=60_000; acc_rate=100)
    target_pdf = x->f(x[1], x[2]; R, prior)
    αcond, βcond = accept_reject(product_distribution([prior.α, prior.β]), target_pdf; acc_rate, nsamples)
    Scond = (αcond * R) ./ βcond
    SIRModel{Float32}(
        S₀=Particles(Float32.(Scond)), 
        β=Particles(Float32.(βcond)), 
        α=Particles(Float32.(αcond))
    )
end

function sir_cond_outbreak_size(prior, size, nsamples=60_000; acc_rate=100)
    target_pdf = x->prob_αS₀_cond_outbreak_size(x[1], x[2], size; prior)
    αcond, Scond = accept_reject(product_distribution([prior.α, prior.S₀]), target_pdf; acc_rate, nsamples)
    βcond = inv_outbreak_size.(αcond, Scond, size)
    SIRModel{Float32}(
        S₀=Particles(Float32.(Scond)), 
        β=Particles(Float32.(βcond)), 
        α=Particles(Float32.(αcond))
    )
end

function inct_joint_func_exper!(config; N=100, M=60_000, acc_rate=100)
    @unpack θtrue, θprior, obs_mod = config
    
    # particle vectors that are constant wrt ϕ:
    lat_mod = SIRModel{Float32}(
        S₀=Particles(M, θprior.S₀), 
        β=Particles(M, θprior.β), 
        α=Particles(M, θprior.α)
    )
    inf_pri = solve(lat_mod; save_idxs=2, saveat=1).u # TODO produce or load all of these solves
    inf_true = solve(lat_mod, θtrue; save_idxs=2, saveat=1).u

    # particle vector simulations for the conditional function models
    samp_func = x->prob_αβ_cond_R(x[1], x[2]; R=rnot(θtrue), prior=θprior)
    Rtrue = rnot(θtrue)
    inf_cond_R = solve(sir_cond_func(θprior, prob_αβ_cond_R, Rtrue, M; acc_rate); save_idxs=2, saveat=1).u

    y = Particles(40_000, observe_dist(obs_mod; observe_params(obs_mod, inf_true)...))
    labs = ["md-joint", "md-R"]
    for l in labs
        config[l] = []
    end
        
    for imax in 2:length(inf_true) # for each timespan (0, t)
        ts_idx = 1:imax
        ynew = bootstrap(y[ts_idx], N)
        push!(config["md-joint"], marginal_divergence(ynew, inf_true[ts_idx], inf_pri[ts_idx], obs_mod))
        push!(config["md-R"], marginal_divergence(ynew, inf_cond_R[ts_idx], inf_pri[ts_idx], obs_mod))
    end
end

addprocs(7)

subfolder = ""
θtrue = (S₀=0.6f0, β=1.25f0, α=0.2f0)
# θtrue = (S₀=1f0, β=0.2f0, α=0.2f0)
θprior = (S₀=Uniform(0.1f0, 1f0), β=Uniform(0.3f0, 3f0), α=Uniform(0.05f0, 0.5f0))
obs_mod = PoissonTests(1000)

factors = @strdict θtrue θprior obs_mod

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

# TODO I'm tired of collect_results and produce_or_load...rewrite alternatives
for d ∈ dict_list(factors)
    fname = mysavename(d)
    inct_joint_func_exper!(d; N=3000)
    tagsave("$fdir/$fname", d; safe)
end
