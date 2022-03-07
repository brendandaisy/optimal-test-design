using Revise
using DEParamDistributions
using DifferentialEquations.EnsembleAnalysis
using OrdinaryDiffEq
using Distributions
using Plots
using Optim
using BenchmarkTools
using Profile

function optimal_design(
        xtrue, demodel::AbstractODEParamDistribution; 
        B=1000, U, jointlik=DEParamDistributions.joint_poisson, M=100
    )

    
end

sse(θ, θest) = -sum((θ .- θest).^2)

function util(
        design, xtrue, demodel::AbstractODEParamDistribution; 
        U, jointlik=DEParamDistributions.joint_poisson, Q=100, M=100 # Q is num prior draws (importance samples) since not part of standard "NxM" loop (i.e. inner posterior loop)
    )
    # precompute sims from prior for importance sampling 
    rvs = random_vars(demodel)
    gdraws = map(_->NamedTuple{keys(rvs)}(rand.(values(rvs))), 1:Q)
    sols = prior_predict(gdraws, demodel)
    ret = mapreduce(+, 1:M) do _ # compute sum of utilities, over y1,...,yM
        y = rand(jointlik(design .* xtrue)) # take a draw from y | xtrue, design
        W = importance_weights(y, sols, jointlik) # jointlik will be called on y here
        μpost = importance_mean(W, gdraws)
        U(μpost) # up to user to specify other needed quants (i.e. true θ)
    end
    ret / M
end

function util(
        design, xtrue; 
        B=1000, U, jointlik=DEParamDistributions.joint_poisson, M=100,
        gdraws, simdists = map(jointlik, sols) # precomputed stuff
    )
    # if sum(design) > B # put the bounds check in this one since this'll be the one to use for optim (for now)
    #     return -10e6
    # end
    if any(design .< 0)
        return -10e6
    end
    ntest = design * B
    ydist = jointlik(ntest .* xtrue)
    ret = zeros(M)
    Threads.@threads for i=1:M
        y = rand(ydist) # take a draw from y | xtrue, design
        W = importance_weights(y, simdists) # jointlik will be called on y here
        μpost = importance_mean(W, gdraws)
        ret[i] = U(μpost) # up to user to specify other needed quants (i.e. true θ)
    end
    mean(ret)
    # ret = mapreduce(+, 1:M) do _ # compute sum of utilities, over y1,...,yM
    #     y = rand(ydist) # take a draw from y | xtrue, design
    #     W = importance_weights(y, simdists) # jointlik will be called on y here
    #     μpost = importance_mean(W, gdraws)
    #     U(μpost) # up to user to specify other needed quants (i.e. true θ)
    # end
    # ret / M
end

## setup and viz inputs for util
pfixed = SIRParamDistribution(60., 0.1, 0.5, 0.1)
prior = SIRParamDistribution(60., 0.1, TruncatedNormal(0.4, 0.5, 0.1, 2), 0.1)
inf_curve = solve(ode_problem(pfixed), Tsit5(), save_idxs=1, saveat=7.).u
pri_curves = prior_predict(prior, 1000; sparse=false, saveat=7.)
plot(EnsembleSummary(pri_curves))
plot!(pri_curves[1].t, inf_curve, linecol="orange", lab="True outbreak")

## double check importance sampling (NOT RUN)
y = rand.(Poisson.(10 * inf_curve))
# in this case, we know only need 1 set of draws for g (prior, for imp sampling) and inf curves (for generating y and imp sampling)
rvs = random_vars(prior)
gdraws = map(_->NamedTuple{keys(rvs)}(rand.(values(rvs))), 1:1_000_000)
sols = prior_predict(gdraws, prior; saveat=7.)

turingfit = sample_mcmc(y, prior, x->DEParamDistributions.array_poisson(10*x); saveat=7.)
mcmc_mean(turingfit)

W = importance_weights(y, sols, x->DEParamDistributions.joint_poisson(10 * x))
sum(W)
importance_mean(W, gdraws)
importance_ess(W)
@btime importance_weights!(W, y, sols, x->DEParamDistributions.joint_poisson(2 * x))

# all seems in order! now I can use this flow in util
@btime util(d0 * 100, inf_curve, prior; U=θ->sse([0.5], θ), Q=1000, M=1000)
simdists = map(DEParamDistributions.joint_poisson, sols);
U = θ->sse([0.5], θ)
Profile.clear()
@btime util(d0, inf_curve; B=100, U, gdraws, simdists, M=1000)
# @btime [DEParamDistributions.logweight(y, DEParamDistributions.joint_poisson(100 * inf_curve)) for i=1:100_000]
Profile.print()

d0 = fill(floor(1 / 10; digits=2), 10) # dumb initial design (make sure < 1)
f = d -> log(-util(d, inf_curve; B=100, U, gdraws, simdists, M=2000))
manif = Optim.Sphere() # issue here is allows negative numbers
f0 = f(d0)
f(fill(10, 10))
f(fill(100, 10))

opt = optimize(f, d0, NelderMead(manifold=manif))
opt.minimum < f(d0)
opt.minimizer