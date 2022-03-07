using Revise
using DEParamDistributions
using OrdinaryDiffEq
using Distributions
using CSV, DataFrames
## For plotting:
# using DifferentialEquations.EnsembleAnalysis
# using Plots
# using StatsPlots

function new_lik(t1, t2, prop, x; budget, array=false)
    n = [budget * prop, budget * (1 - prop)]
    if array
        return array_poisson(n, [x[t1], x[t2]])
    end
    joint_poisson(n, [x[t1], x[t2]])
end

new_lik(d::Tuple, x; budget=100, array=false) = new_lik(d..., x; budget, array)

function uofp_exper(dvec, xtrue, pdist::AbstractODEParamDistribution; 
    N=10_000, M=100, reps=1, θtrue, dekwargs...
)
    ## precomputations for importance sampling (nse)
    postfit = initial_fit(
        pdist.start:pdist.stop, xtrue, pdist; 
        alik=DEParamDistributions.array_poisson, dekwargs...
    )
    rvs = random_vars(pdist)
    postdraw = map(eachcol(rand(postfit, N))) do draw
        NamedTuple{keys(rvs)}(draw)
    end
    gsims = prior_predict(postdraw, pdist; dekwargs...)
    pd = joint_prior(pdist)
    gdraw = map(x->vcat(values(x)...), postdraw)

    ## precomputations for importance sampling (nse)
    
    pdraws = map(_->NamedTuple{keys(rvs)}(rand.(values(rvs))), 1:N)
    psims = prior_predict(pdraws, pdist; dekwargs...)
    
    # compute SIG/SSE|xtrue for each design
    ret = Vector{Tuple}(undef, length(dvec) * reps)
    Threads.@threads for i=eachindex(dvec)
        design = dvec[i]
        ld_is = [DEParamDistributions.new_lik(design, gx) for gx ∈ gsims]
        ld_me = [DEParamDistributions.new_lik(design, px) for px ∈ psims]
        for r=1:reps
            ret[(i-1)*reps+r] = (design..., util(design, xtrue; gdraw, ld=ld_is, pd, gd=postfit, M, θtrue), 
                util2(design, xtrue; ld=ld_me, M), r)
        end
    end
    ret
end

function util_grid_exper(
    dgrid, xtrue, pdist::AbstractODEParamDistribution;
    N=10_000, Ng=1500, M=100, θtrue=[0.3, 0.7], dekwargs...
)
    ## precomputations for importance sampling (nse)
    postfit = initial_fit(
        pdist.start:pdist.stop, xtrue, pdist; 
        alik=DEParamDistributions.array_poisson, dekwargs...
    )
    rvs = random_vars(pdist)
    ## Get Ng samples from the target distribution g
    postdraw = map(eachcol(rand(postfit, Ng))) do draw
        NamedTuple{keys(rvs)}(draw)
    end
    gsims = prior_predict(postdraw, pdist; dekwargs...)
    pd = joint_prior(pdist)
    gdraw = map(x->vcat(values(x)...), postdraw)

    # setup simulation curves for all designs
    rvs = random_vars(pdist)
    pdraws = map(_->NamedTuple{keys(rvs)}(rand.(values(rvs))), 1:N)
    psims = prior_predict(pdraws, pdist; dekwargs...)

    # compute SIG/SSE|xtrue for each design
    ret = Vector{Tuple}(undef, length(dgrid))
    sig = zeros(M)
    nsse = zeros(M)
    for i=eachindex(dgrid) # each design
        design = dgrid[i]
        ld_is = [DEParamDistributions.new_lik(design, gx) for gx ∈ gsims] # Ng distributions for IS
        ld_me = [DEParamDistributions.new_lik(design, px) for px ∈ psims] # N dist. for MI
        Threads.@threads for j=1:M # expectation over ys
            y = rand(DEParamDistributions.new_lik(design, xtrue))
            sig[j] = get_sig(y, design, xtrue; ld=ld_me)
            nsse[j] = get_nsse(y; gdraw, ld=ld_is, pd, gd=postfit, θtrue)
        end
        ret[i] = (design..., mean(nsse), mean(sig))
    end
    ret
end

## Setup
pfixed = SIRParamDistribution(60., 0.3, 0.7, 0.1)
pdist = SIRParamDistribution(60., Beta(1.4, 6), TruncatedNormal(0.4, 1, 0.15, 3), 0.1)

## GRID EXPERIMENT

range = 1:3:pdist.stop
design_grid = [(Int(i), Int(j), 0.5) for i=range for j=range if j > i]
𝔾 = util_grid_exper(design_grid, inf_curve, pdist; N=1000, Ng=1000, M=50, save_idxs=1, saveat=1)
# CSV.write("output/sig-pairs-2-4.csv", DataFrame(𝔾), append=true)

## UofP EXPERIMENT + MC ERROR
# The comparison to sse, and testing error, can be made in here

# props = 0:0.1:1
# dofp = [(10, 20, p) for p ∈ props]
# 𝔾 = uofp_exper(dofp, inf_curve, pdist; N=100_000, M=1000, reps=1, save_idxs=1, saveat=1, θtrue=[pfixed.rec_init, pfixed.inf_rate])
# CSV.write("uofp-2-3.csv", DataFrame(𝔾), header=false)

## CONFIRMATION FIGURES
## what do posteriors look like for some good designs?
# des = (10, 20, 0.4)
# lf(x) = DEParamDistributions.new_lik(des, x; array=true)
# y = rand.(lf(inf_curve))
# fit = sample_mcmc(y, pdist, lf; save_idxs=1, saveat=1, iter=1000)
# using Turing
# postsim = generated_quantities(turingode(y, pdist, lf; save_idxs=1, saveat=1), fit)
# postsim = reshape(postsim, (4000,))

# uofp() = [util2((5, 20, p), inf_curve, sols; M=800) for p ∈ 0:0.1:1]
# utils = [uofp() for _ in 1:50]

# CSV.write("SIGofP-2-22.txt", DataFrame(utils, :auto), header=false)
# um = vecvec_to_mat(utils)'
# μ_um = mapslices(mean, um; dims=2)
# plot(0:0.1:1, μ_um)
# plot!(0:0.1:1, mapslices(minimum, um; dims=2))
# plot!(0:0.1:1, mapslices(maximum, um; dims=2))
# savefig("mc-error-sigofp.pdf")