using DrWatson
@quickactivate "optimal-test-design"
using Revise
using DiffEqInformationTheory
using ConditionalTransform
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

θtrue = (α=0.2f0, β=1.25f0, S₀=0.6f0)
# θtrue = (S₀=1f0, β=0.2f0, α=0.2f0)
θprior = (α=Uniform(0.05f0, 0.5f0), β=Uniform(0.3f0, 3f0), S₀=Uniform(0.1f0, 1f0))

#=Basic Reproductive Number=#
Rtrue = rnot(values(θtrue)...)
αcond, βcond, Scond = sample_cond_f(θprior, Rtrue, inv_rnot, d_inv_rnot; pivot=:S₀, nsamples=40_000)
sir_cond = SIRModel{Float32}(
    S₀=Particles(Float32.(Scond)), 
    β=Particles(Float32.(βcond)), 
    α=Particles(Float32.(αcond))
)
inf_cond = solve(resample(sir_cond, 500); save_idxs=2, saveat=0.2)

save_samples("rep-number", αcond, βcond, Scond)
save_sims("rep-number", inf_cond)

#= Outbreak Size =#
osize = outbreak_size(θtrue.α, θtrue.β, θtrue.S₀)

αcond, βcond, Scond = sample_cond_f(
    θprior, osize, inv_outbreak_size, d_inv_outbreak_size; 
    pivot=:β, cond=(α, S₀, size) -> S₀+0.01 < size, nsamples=40_000
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
    pivot=:S₀, nsamples=40_000
)
sir_cond = SIRModel{Float32}(
    S₀=Particles(Float32.(Scond)), 
    β=Particles(Float32.(βcond)), 
    α=Particles(Float32.(αcond))
)
inf_cond = solve(resample(sir_cond, 500); save_idxs=2, saveat=0.2)

save_samples("peak-intensity", αcond, βcond, Scond)
save_sims("peak-intensity", inf_cond)

#= Initial growth rate =#
grate = growth_rate(θtrue...)

αcond, βcond, Scond = sample_cond_f(
    θprior, grate, inv_growth_rate, d_inv_growth_rate; 
    pivot=:α, nsamples=40_000
)
sir_cond = SIRModel{Float32}(
    S₀=Particles(Float32.(Scond)), 
    β=Particles(Float32.(βcond)), 
    α=Particles(Float32.(αcond))
)
inf_cond = solve(resample(sir_cond, 500); save_idxs=2, saveat=0.2)

save_samples("growth-rate", αcond, βcond, Scond)
save_sims("growth-rate", inf_cond)