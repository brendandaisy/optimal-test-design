using DrWatson
@quickactivate "optimal-test-design"
using Revise
using DiffEqInformationTheory
using Distributions, MonteCarloMeasurements
using NamedTupleTools
using Distributed

using NLsolve

include(srcdir("watson-tools.jl"))

function find_sol(f, x_init=1.)
    sol = nlsolve(x->[x[1] - f(x[1])], [x_init]; iterations=5_000)
    if (!sol.f_converged)
        @warn "Solver for implicit function did not converge; do not trust results"
    end
    sol.zero[1]
end

rnot(α, β, S₀) = S₀*β/α
rnot(θ) = rnot(θ.α, θ.β, θ.S₀)

# rec_inf(α, β, S₀; r_inf, I₀=0.01) = 1 - S₀*exp(-β/α * (r_inf+S₀+I₀-1))
function outbreak_size(α, β, S₀; I₀=0.01)
    sol = solve(SIRModel{Float64}(;stop=10_000, α, β, S₀, I₀); save_idxs=1)
    S₀ + I₀ - last(sol.u)
end
# inv_rec_inf(α, β, r_inf; S₀, I₀=0.01) = (1 - r_inf) / exp(-β/α*(r_inf + S₀ + I₀ - 1))

function inv_outbreak_size(α, S₀, size; I₀=0.01)
    R₀ = 1 - S₀ - I₀
    -α / size * log((1-R₀-size) / S₀)
end

function d_inv_outbreak_size(α, S₀, size; I₀=0.01)
    R₀ = 1 - S₀ - I₀
    e1 = 1 - R₀ - size
    α / size * (1 / size * log(e1/S₀) + 1/e1)
end

inv_outbreak_size(0.2, 0.6, 0.609999999999999)

d_inv_outbreak_size(0.2, 0.6, 0.5928956712760867)

S₀sol = find_sol(x->inv_rec_inf(0.2, 1.2, 0.7; S₀=x))

rec_inf(0.2, 1.2, 0.6; r_inf=0.9)

rec_inf_ode(0.2, 1.2, 0.6)

find_sol(x->inv_rec_inf(0.2, 1.2, 0.4; S₀=x))
find_sol(x->rec_inf(0.2, 1.2, 0.6; r_inf=x))

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

function accept_reject(sampler, target; acc_rate=100, nsamples=100)
    ret = []
    attempts = 0
    for _ in 1:nsamples
        s = rand(sampler)
        u = rand()
        attempts += 1
        while u >= target(s) / (acc_rate*pdf(sampler, s))
            s = rand(sampler)
            u = rand()
            attempts += 1
        end
        push!(ret, s)
    end
    @info "Acceptance rate was $(nsamples/attempts)"
    return (first.(ret), last.(ret))
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
θtrue = (S₀=1f0, β=0.2f0, α=0.2f0)
θprior = (S₀=Uniform(0.1f0, 0.9f0), β=Uniform(0.3f0, 3f0), α=Uniform(0.05f0, 0.3f0))
# known = [(), (:α,)]
known = ()
obs_mod = PoissonTests(1000)

solve(SIRModel{Float32}(;θtrue...))

outbreak_size(sir_cond.α.particles[12], sir_cond.β.particles[12], sir_cond.S₀.particles[12])

# using Plots
osize = outbreak_size(θtrue.α, θtrue.β, θtrue.S₀)
sir_cond = sir_cond_outbreak_size(θprior, osize)
inf_cond = solve(sir_cond; save_idxs=2, saveat=1)
ribbonplot(inf_cond.t, inf_cond.u, 0.05)
mcplot!(inf_cond.t, inf_cond.u, 200, alpha=0.3)

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
