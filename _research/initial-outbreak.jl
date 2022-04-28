using Distributions, DEParamDistributions
using OrdinaryDiffEq
using Random
using DataFrames
using RCall

Random.seed!(1234)

θtrue = (S₀=0.6, β=1.25, α=0.2)
pdist = SIRParamDistribution(;S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=0.2)
lik(x, t=3) = obs_tspan(x, PoissonTests(1000.), t)

xtrue = solve_de_problem(pdist, θtrue; saveat=1, save_idxs=2).u
y = rand(lik(xtrue))

pri_draws = map(1:500) do _
    s, b = rand(joint_prior(pdist))
    (S₀=s, β=b)
end
sim = simulate(pri_draws, pdist; saveat=1, save_idxs=2)
liks = [pdf(lik(x), y) for x ∈ sim]
y2 = rand(lik(xtrue, 8))
liks2 = [pdf(lik(x, 8), y2) for x ∈ sim]

likdf = DataFrame(pri_draws)
likdf.lik3 = liks
likdf.lik8 = liks2
@rput likdf