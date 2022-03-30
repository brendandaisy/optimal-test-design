using Distributions, DEParamDistributions
using OrdinaryDiffEq
using Plots, StatsPlots
using Random
using DataFrames
using RCall
using DiffEqSensitivity, ForwardDiff

Random.seed!(1234)

θtrue = (S₀=0.6, β=1.25, α=0.2)
θprior = (S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
dekwargs = (save_idxs=2, saveat=0.1)
ntest = 100
rate = 10

pdist = SIRParamDistribution(;θprior...)
xtrue = solve(de_problem(pdist, θtrue; dekwargs...), Tsit5()).u
sim = simulate(pdist, 500; keep=true, dekwargs...)

qhi = EnsembleSummary(sim).qhigh
qlo = EnsembleSummary(sim).qlow
subsim = [sim[i].u for i ∈ sample(1:500, 50)]
pri_curves = DataFrame(subsim, :auto)
curves.t = sim[1].t

y=rand(joint_neg_binom(rate, xtrue[1:10:end] .* ntest))
pri_dist = [joint_neg_binom(rate, s[1:10:(12*10)] .* ntest) for s ∈ subsim]
pri_lpdf = [logpdf(dist, y[1:12]) for dist ∈ pri_dist]
@rput qhi qlo pri_curves pri_lpdf


βpd = SIRParamDistribution(S₀=Uniform(0.1, 0.9), β=1.25, α=Uniform(0.05, 0.3))
βcurves = simulate(βpd, 50; saveat=0.1, save_idxs=2)
βdist = [joint_neg_binom(rate, s[1:10:(12*10)] .* ntest) for s ∈ βcurves]
βlpdf = [logpdf(dist, y[1:12]) for dist ∈ βdist]
bdf = DataFrame(βcurves.u, :auto)
@rput y bdf βlpdf

function get_post(tmax, y=rand(joint_neg_binom(rate, xtrue[1:10:(tmax*10)] .* ntest)))
    sample_mcmc(
        y, pdist, x->joint_neg_binom(rate, x[1:10:((tmax-1)*10+1)] .* ntest); dekwargs..., 
        arraylik=false, iter=1000, ŷ=true
    )
end

function get_sig(tmax, y=rand(joint_neg_binom(rate, xtrue[1:tmax] .* ntest)))
    sims = simulate(pdist, 4000; dekwargs...)
    pri_ldists = map(x->joint_neg_binom(rate, x[1:tmax] .* ntest), sims)
    DEParamDistributions.sig(y, pri_ldists, joint_neg_binom(rate, xtrue[1:tmax] .* ntest))
end

y=rand(joint_neg_binom(rate, xtrue[1:10:end] .* ntest))
fit6, _c6 = get_post(6, y[1:6])
fit12, _c12 = get_post(12, y[1:12])

c6 = DataFrame(sample(reshape(_c6, :), 50), :auto)
c6.t = sim[1].t
c12 = DataFrame(sample(reshape(_c12, :), 50), :auto)
c12.t = sim[1].t

@rput c6 c12 y xtrue

S0_samp = vec(fit12.value[:,1,:])
beta_samp = vec(fit12.value[:,2,:])
alpha_samp = vec(fit12.value[:,3,:])
@rput S0_samp beta_samp alpha_samp

## Sensitivity plot

pdist = SIRParamDistribution(;S₀=Uniform(0.1, 0.9), β=Uniform(0.3, 3), α=Uniform(0.05, 0.3))
u₀ = match_initial_values(pdist, θtrue)
p = match_parameters(pdist, θtrue)
prob = ODEProblem(fiip,u0,(0.0,10.0),p)

true_prob = de_problem(pdist, θtrue)
y = rand(lik(solve_de_problem(pdist, θtrue; saveat=1, save_idxs=2).u, 31))

function ℓoft(x, y, t)
    _prob = remake(true_prob, u0=[x[1], pdist.I₀], p=x[2:end])
    sol = solve(_prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=1, save_idxs=2)
    logpdf(lik(sol.u, t), y[1:t])
end

lik(x, t) = joint_poisson(1000 * x[1:t])

dels = []
for i=1:250_000
    y = rand(lik(solve_de_problem(pdist, θtrue; saveat=1, save_idxs=2).u, 31))
    ∇ℓoft = map(1:31) do t
        ForwardDiff.gradient(x->ℓoft(x, y, t), vcat(θtrue.S₀, p))
    end
    push!(dels, ∇ℓoft)
end

μdel = reduce((a, b)->a .+ b, dels) ./ length(dels)
plot(vecvec_to_mat(μdel))

function inf_of_t(x, t)
    _prob = remake(true_prob, u0=[x[1], pdist.I₀], p=x[2:end])
    sol = solve(_prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=0.1, save_idxs=2)
    sol.u[t]
end

∇inf = map(1:301) do t
    ForwardDiff.gradient(x->inf_of_t(x, t), vcat(θtrue.S₀, p))
end

plot(vecvec_to_mat(∇inf))

sens_prob = ODEForwardSensitivityProblem(de_func(pdist), u₀, timespan(pdist), p)
sol = solve(sens_prob, Tsit5())
x, dp = extract_local_sensitivities(sol)
