using DrWatson
@quickactivate "optimal-test-design"
using Revise
using DiffEqInformationTheory
using ConditionalTransform
using Distributions, MonteCarloMeasurements
using NamedTupleTools

include(srcdir("watson-tools.jl"))
include(srcdir("transform-funcs.jl"))

function mean_ml(θtrue, θprior, obs_mod; N=100, M=60_000)
    lat_mod = SIRModel{Float32}(
        S₀=Particles(M, θprior.S₀), 
        β=Particles(M, θprior.β), 
        α=Particles(M, θprior.α)
    )
    inf_pri = solve(lat_mod; save_idxs=2, saveat=1).u
    inf_true = solve(lat_mod, θtrue; save_idxs=2, saveat=1).u
    y = Particles(N, observe_dist(obs_mod; observe_params(obs_mod, inf_true)...))
    inf_conds = solve_inf_cond(θtrue, θprior, lat_mod)
    
    ret = []
    for tmax in [3, 8]
        @show tmax
        ts_idx = 1:(tmax+1)
        Eℓpy = mean_marginal_likelihood(y, inf_pri[ts_idx], obs_mod)
        push!(ret, (lab="prior", t=tmax, eml=Eℓpy))
        for lab in keys(inf_conds)
            Eℓpy_cond = mean_marginal_likelihood(y, inf_conds[lab][ts_idx], obs_mod)
            push!(ret, (lab=lab, t=tmax, eml=Eℓpy_cond))
        end
    end
    return ret
end

θtrue = (α=0.2f0, β=1.25f0, S₀=0.6f0)
θprior = (α=Uniform(0.05f0, 0.85f0), β=Uniform(0.3f0, 1.5f0), S₀=Uniform(0.1f0, 0.99f0))
obs_mod = PoissonTests(1000)

res = mean_ml(θtrue, θprior, obs_mod; N=10, M=60_000)
CSV.write(datadir("sims", "ex-marg-likelihoods.csv"), DataFrame(res))
