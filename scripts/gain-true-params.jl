using DrWatson
@quickactivate "optimal-test-design"
using Revise
using DiffEqInformationTheory
using ConditionalTransform
using Distributions, MonteCarloMeasurements
using NamedTupleTools
using Optim

include(srcdir("watson-tools.jl"))
include(srcdir("transform-funcs.jl"))

function sir_from_samples(α, β, S₀)
    SIRModel{Float32}(
        S₀=Particles(Float32.(S₀)), 
        β=Particles(Float32.(β)), 
        α=Particles(Float32.(α))
    )
end

function solve_inf_cond(θtrue, θprior, lat_mod)
    nsamples=nparticles(lat_mod.β)
    inf_conds = Dict()
    #= Peak intensity =#
    imax = max_inf(θtrue...)
    if imax > lat_mod.I₀
        αcond, βcond, Scond = sample_cond_f(
            θprior, imax, inv_max_inf, d_inv_max_inf; 
            pivot=:S₀, nsamples, silent=true
        )
    else
        αcond, βcond, Scond = sample_trunc_f(θprior, (a, b, S)->reff(a, b, S) > 1; nsamples)
    end
    sir_cond = sir_from_samples(αcond, βcond, Scond)
    inf_conds["peak-intensity"] = solve(sir_cond; save_idxs=2, saveat=1).u

    #= Simulations for ind. variables =#
    for θᵢ in keys(θtrue) # setup for ind params
        lab = string(θᵢ)
        θcond = NamedTupleTools.select(θtrue, (θᵢ,))
        inf_conds[lab] = solve(lat_mod, θcond; save_idxs=2, saveat=1).u # simulate x∣θᵢ
    end

    #= Basic Reproductive Number =#
    Rtrue = rnot(θtrue.α, θtrue.β)
    αcond, βcond = sample_cond_f(NamedTupleTools.select(θprior, (:α, :β)), Rtrue, inv_rnot, d_inv_rnot; pivot=:β, nsamples, silent=true)
    sir_cond = sir_from_samples(αcond, βcond, rand(θprior.S₀, nsamples))
    inf_conds["rep-number"] = solve(sir_cond; save_idxs=2, saveat=1).u

    #= Outbreak Size =#
    osize = outbreak_size(θtrue.α, θtrue.β, θtrue.S₀)
    αcond, βcond, Scond = sample_cond_f(
        θprior, osize, inv_outbreak_size, d_inv_outbreak_size; 
        pivot=:β, cond=(α, S₀, size) -> S₀+0.01 < size, nsamples, silent=true
    )
    sir_cond = sir_from_samples(αcond, βcond, Scond)
    inf_conds["outbreak-size"] = solve(sir_cond; save_idxs=2, saveat=1).u

    #= Initial growth rate =#
    grate = growth_rate(θtrue...)
    αcond, βcond, Scond = sample_cond_f(
        θprior, grate, inv_growth_rate, d_inv_growth_rate; 
        pivot=:α, nsamples, silent=true
    )
    sir_cond = sir_from_samples(αcond, βcond, Scond)
    inf_conds["growth-rate"] = solve(sir_cond; save_idxs=2, saveat=1).u
    
    return inf_conds
end

function ident_over_θtrue(θtrue_set, θprior, obs_mod; N=100, M=60_000)
    lat_mod = SIRModel{Float32}(
        S₀=Particles(M, θprior.S₀), 
        β=Particles(M, θprior.β), 
        α=Particles(M, θprior.α)
    )
    inf_pri = solve(lat_mod; save_idxs=2, saveat=1).u

    ret = []
    for θtrue in θtrue_set
        @show θtrue
        sol_true = solve(lat_mod, θtrue; save_idxs=2, dense=true)
        f = t -> -sol_true(t)
        peak_t = ceil(Int, optimize(f, 0, 30).minimizer)

        inf_true = sol_true(0:peak_t).u
        y = Particles(10_000, observe_dist(obs_mod; observe_params(obs_mod, inf_true)...))
        inf_conds = solve_inf_cond(θtrue, θprior, lat_mod)
        
        ts_idx = 1:(peak_t+1)
        for lab in keys(inf_conds)
            ynew = bootstrap(y, N)
            δ = marginal_divergence(ynew, inf_conds[lab][ts_idx], inf_pri[ts_idx], obs_mod)
            push!(ret, merge(θtrue, (lab=lab, t=peak_t, md=δ)))
        end
    end
    return ret
end

θprior = (α=Uniform(0.05f0, 0.85f0), β=Uniform(0.3f0, 1.5f0), S₀=Uniform(0.1f0, 0.99f0))
# θtrue_set = [(α=0.2, β=b, S₀=S) for b in Float32.(0.3:0.2:1.5) for S in Float32.(0.1:0.2:0.9)]
θtrue_set = [(α=0.2f0, β=0.7f0, S₀=0.3f0), (α=0.2f0, β=0.3f0, S₀=0.9f0), (α=0.2f0, β=0.5f0, S₀=0.5f0)]
# tmp = solve(SIRModel{Float32}(α=0.2f0, β=0.3, S₀=0.7), save_idxs=2, dense=true)
# plot(tmp)
# θtrue_set = filter(x->reff(x...) > 0.9, θtrue_set)
obs_mod = PoissonTests(1000)

res = ident_over_θtrue(θtrue_set, θprior, obs_mod; N=3000, M=60_000)
CSV.write(datadir("sims", "gain-by-peak-check-out.csv"), DataFrame(res))
