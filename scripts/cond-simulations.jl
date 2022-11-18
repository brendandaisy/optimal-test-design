using DrWatson
@quickactivate "optimal-test-design"
using Revise
using DiffEqInformationTheory
using ConditionalTransform
using NamedTupleTools
using Distributions, MonteCarloMeasurements
using DataFrames, CSV

include(srcdir("transform-funcs.jl"))

function save_samples(varname, α, β, S₀)
    CSV.write(
        datadir("sims", "cond-samples", "samp-$(varname).csv"),
        DataFrame([α, β, S₀], ["alpha", "beta", "S0"])
    )
end

function save_sims(varname, inf)
    CSV.write(
        datadir("sims", "cond-simulations", "sim-$(varname).csv"),
        DataFrame(Matrix(inf.u), string.(inf.t))
    )
end

# DANGER: these values must match the ones in daily-obs script!!
θtrue = (α=0.2f0, β=1.25f0, S₀=0.6f0)
# θtrue = (S₀=1f0, β=0.2f0, α=0.2f0)
θprior = (α=Uniform(0.05f0, 0.85f0), β=Uniform(0.3f0, 1.5f0), S₀=Uniform(0.1f0, 0.99f0))

asamp, bsamp, Ssamp = rand.(values(θprior), 5000)
Rpri = reff.(asamp, bsamp, Ssamp)
sum(Rpri .> 1) / 5000

#= Simulations for ind. variables =#
lat_mod = SIRModel{Float32}(
    S₀=Particles(500, θprior.S₀), 
    β=Particles(500, θprior.β), 
    α=Particles(500, θprior.α)
)
for θᵢ in keys(θtrue) # setup for ind params
    θcond = NamedTupleTools.select(θtrue, (θᵢ,))
    inf_cond = solve(lat_mod, θcond; save_idxs=2, saveat=1) # simulate x∣θᵢ
    save_sims(string(θᵢ), inf_cond)
end

#= Sample a single y|θtrue for post fitting =#
inf_true = solve(lat_mod, θtrue; save_idxs=2, saveat=1).u
ydist = observe_dist(obs_mod; observe_params(obs_mod, inf_true)...)
ysamp = rand(ydist)
println(ysamp)

#= Basic Reproductive Number =#
Rtrue = rnot(θtrue.α, θtrue.β)
αcond, βcond = sample_cond_f(NamedTupleTools.select(θprior, (:α, :β)), Rtrue, inv_rnot, d_inv_rnot; pivot=:β, nsamples=60_000)
sir_cond = SIRModel{Float32}(
    S₀=Particles(60_000, θprior.S₀), 
    β=Particles(Float32.(βcond)), 
    α=Particles(Float32.(αcond))
)
inf_cond = solve(resample(sir_cond, 500); save_idxs=2, saveat=0.2)

save_samples("rep-number", αcond, βcond, sir_cond.S₀.particles)
save_sims("rep-number", inf_cond)

#= Outbreak Size =#
osize = outbreak_size(θtrue.α, θtrue.β, θtrue.S₀)

αcond, βcond, Scond = sample_cond_f(
    θprior, osize, inv_outbreak_size, d_inv_outbreak_size; 
    pivot=:β, cond=(α, S₀, size) -> S₀+0.01 < size, nsamples=60_000
)
sir_cond = SIRModel{Float32}(
    S₀=Particles(Float32.(Scond)), 
    β=Particles(Float32.(βcond)), 
    α=Particles(Float32.(αcond))
)
inf_cond = solve(resample(sir_cond, 500); save_idxs=2, saveat=1)

save_samples("outbreak-size", αcond, βcond, Scond)
save_sims("outbreak-size", inf_cond)

#= Peak intensity =#
imax = max_inf(θtrue...)

αcond, βcond, Scond = sample_cond_f(
    θprior, imax, inv_max_inf, d_inv_max_inf; 
    pivot=:S₀, nsamples=60_000
)
αcond, βcond, Scond = sample_trunc_f(θprior, (a, b, S)->reff(a, b, S) > 1; nsamples=60_000)
sir_cond = SIRModel{Float32}(
    S₀=Particles(Float32.(Scond)), 
    β=Particles(Float32.(βcond)), 
    α=Particles(Float32.(αcond))
)
inf_cond = solve(resample(sir_cond, 500); save_idxs=2, saveat=0.2)

save_samples("peak-intensity2", αcond, βcond, Scond)
save_sims("peak-intensity2", inf_cond)

#= Initial growth rate =#
grate = growth_rate(θtrue...)

αcond, βcond, Scond = sample_cond_f(
    θprior, grate, inv_growth_rate, d_inv_growth_rate; 
    pivot=:α, nsamples=60_000
)
sir_cond = SIRModel{Float32}(
    S₀=Particles(Float32.(Scond)), 
    β=Particles(Float32.(βcond)), 
    α=Particles(Float32.(αcond))
)
inf_cond = solve(resample(sir_cond, 500); save_idxs=2, saveat=0.2)

save_samples("growth-rate", αcond, βcond, Scond)
save_sims("growth-rate", inf_cond)