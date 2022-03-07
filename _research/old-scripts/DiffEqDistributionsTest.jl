using Test

include("DiffEqDistributions.jl")

pop_size = 10.
rec_rate = 0.1
sir_params = (
    start=0.0, stop=30.0, inf_init=1.0, rec_rate=rec_rate, pop_size=pop_size,
    rec_init=truncated(Normal(0.5*pop_size, 0.25*pop_size), 0, pop_size), 
    inf_rate=truncated(Normal(rec_rate, 1.5), 0, Inf)
)
ensemble_prob = EnsembleProblem(blankSIR(), prob_func=prob_func, output_func=output_func)
sim = solve(ensemble_prob, trajectories=10)

@test blankSIR() isa ODEProblem
@test keys(param_sample(sir_params)) == keys(sir_params)
@test [getindex(sir_params, i) for i=1:(lastindex(sir_params)-2)] == 
    [getindex(param_sample(sir_params), i) for i=1:(lastindex(sir_params)-2)]
@test let sir = blankSIR()
    sir !== prob_func(sir, 1, 1)
end
@test let sir = blankSIR()
    ensemble_prob = EnsembleProblem(sir, prob_func=prob_func, output_func=output_func)
    sim = solve(ensemble_prob, trajectories=2)
    sim[1][2].prob.tspan == (sir_params.start, sir_params.stop)
end
@test let vv = fill([1, 2, 3], 5)
    pluck_vecvec(vv) == fill(1, 5)
    pluck_vecvec(vv, 1:2) == fill([1, 2], 5)
end
ensemble_timeseries(sim, 0:0.1:10)

ts_pois = ensemble_timeseries(sim, 1., Poisson)

## TODO don't nec need to have lost res on the series part
plot([x.series for x ∈ ts_pois])
scatter!([x.y for x ∈ ts_pois])
