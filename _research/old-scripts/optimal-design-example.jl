using Revise
using DEParamDistributions
using Distributions
using OrdinaryDiffEq
using Plots
using StatsPlots
using Turing
import IterTools: fieldvalues

function expected_util(η, pfixed, prior::AbstractODEParamDistribution; M=100, cache=true)
    prob = ode_problem(pfixed)
    θ = [getproperty(pfixed, x) for x ∈ fieldnames(random_vars(prior))] # params that go into utility
    inf_curve = solve(prob, Tsit5(), save_idxs=1, saveat=1.).u
    ret = 0
    fits = []
    for i=1:M
        y = rand.(Binomial.(η, inf_curve))

        ret += U1(θ, fit)
    end
    if cache
        return ret / M, fits
    end
    return ret / M
end

function U1(θ, fit)
    mean_post = summarystats(fit).nt.mean
    err = (θ .- mean_post).^2 |> sum
    -err
end

## Run some prior predictive viz

## idea for better check with PI PT: think of as regression, and instread of heatmap do regression line with 95 CI, plus random subset of pts

## Test some utility funs
pfixed = SIRParamDistribution(60., 0.7, 0.8, 0.1)
prior = SIRParamDistribution(60., 0.7, TruncatedNormal(0.4, 0.5, 0.1, 2), 0.1)
plot(prior.inf_rate)
inf_curve = solve(ode_problem(pfixed), Tsit5(), save_idxs=1, saveat=1.).u
plot(inf_curve)
y = rand.(Binomial.(100, inf_curve))
mod = odemodel(y, prior, x->DEParamDistributions.joint_binom(100, x), )
fit = sample(mod, NUTS(500, 0.65), MCMCThreads(), 200, 2)

summarystats(fit).nt.mean
U1([0.8], fit)

## Do things for good priors

U, fits = expected_util(fill(100, 61), pfixed, prior, M=10, cache=true)
means = [summarystats(f).nt.mean for f ∈ fits]
base_prob = ode_problem(pfixed)
post_curves = map(means) do m
    base_prob.p[1] = m[1]
    solve(base_prob, Tsit5(), save_idxs=1, saveat=1.).u
end

plot(inf_curve)
plot!(post_curves, linealpha=0.15, lab="")

bad_eta = sample([1, 2, 4], 61)
U2, fits2 = expected_util(bad_eta, pfixed, prior, M=10, cache=true)
means = [summarystats(f).nt.mean for f ∈ fits2]
base_prob = ode_problem(pfixed)
post_curves = map(means) do m
    base_prob.p[1] = m[1]
    solve(base_prob, Tsit5(), save_idxs=1, saveat=1.).u
end

plot(inf_curve)
plot!(post_curves, linealpha=0.15, lab="")

## GAMEPLAN FOR TOMORROW
## Layout for prior predict plots
## clean up prior-predict to satisfy these plots
## come up with a better set of priors
## demonstrate spreading out tests is better than placing them all before peak

# using Statistics
# using Plots
# using Lazy

# include("DiffEqDistributions.jl")

# T = 100
# B = 100
# design = round.(Int, 10 .* randn(T) .+ 100)

# 𝑈(𝜂, 𝑦, 𝑥) = sum((𝑦 ./ 𝜂 .- 𝑥).^2)

# pop_size = 1.
# rec_rate = 0.1
# sir_params = (
#     start=1., stop=T, inf_init=0.01, rec_rate=rec_rate, pop_size=pop_size,
#     rec_init=truncated(Normal(0.5*pop_size, 0.25*pop_size), 0, pop_size), 
#     inf_rate=truncated(Normal(rec_rate, 1.5), 0, Inf)
# )
# epi_curves = de_distribution(10_000, sir_params, blankSIR)

# function 𝑈(𝜂, epi_curves)
#     ret = 0
#     et = ensemble_timeseries(epi_curves, 1., (p, n) -> Binomial(n, p), 𝜂)
#     for e ∈ et
#         ret += 𝑈(𝜂, e.y, e.series)
#     end
#     ret / length(et)
# end

# 𝑈(η, epi_curves)

# odes = [sim[2] for sim ∈ epi_curves];
# expect_inf = map(1:T) do t
#     inf_set_t = componentwise_vectors_timepoint(odes, t)[1]
#     mean(inf_set_t .* (1 .- inf_set_t))
# end

# plot(expect_inf)

# sum(expect_inf ./ η)