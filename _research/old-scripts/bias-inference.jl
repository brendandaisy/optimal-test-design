using CSV, DataFrames
using Dates
using Plots
using RCall
using Distributions
using Revise
using DEParamDistributions
using DifferentialEquations.EnsembleAnalysis
using Lazy

n1(n, b) = n^b

plot(n1.(1:1000, .9))
plot!(1:1000, linecol="gray")

# bias parameters known, single simulated x=(x1,...,xT)
function U1(x::Vector, y, n; b=1/2, z=3/4)
    n₁ = n1.(n, b)
    x̂ = @. (y - n₁ * z) / (n - n₁)
    Us = @. -(x - x̂)^2
    sum(Us)
end

# bias parameters known, analytic result
function U1(μx, n; b=1/2, z=3/4)
    n₁ = n1.(n, b) # number of biased tests per day
    Us = @. -(n₁ * z + (n - n₁) * μx) / (n - n₁)^2
    sum(Us)
end

function optimal_design(pdist::AbstractODEParamDistribution)
    xs = prior_predict(pdist, N = 100)
    ys = predict_yhat
    μx = EnsembleAnalysis.componentwise_vectors_timestep(xs, )
end

sir_pdist = SIRParamDistribution(60., 0.7, TruncatedNormal(0.4, 0.5, 0.1, 2), 0.1)

(n₁ * z + (n - n₁) * x) / (n - n₁) ^ 2
xs = prior_predict(sir_pdist, 100)
plot(xs.u, lab="")
𝔏 = (x, z, n1, n2) -> DEParamDistributions.array_poisson(n1 .* z .+ n2 .* x)
d = fill(100, 61)
rand.(𝔏(xs[1], 3/4, n1.(d, 1/2), d .- n1.(d, 1/2)))
@run ys = predict_yhat(xs, x->𝔏(x, 3/4, n1.(d, 1/2), d .- n1.(d, 1/2)))
mean([U1(x, y, d) for (x, y) in zip(xs, ys)])

(ys[1] - n1(100, 1/2) * )

μx = timestep_mean(xs, 1:61)
U1(μx, d)

# bias(n, N, b) = (N / n)^b

# posrate(x, n; N=1000, b=1/2) = x*n*bias(n, N, b)

jhu = CSV.read("us-states.csv", DataFrame)
@rput jhu
R"""
library(tidyverse)
# jhu_ts = jhu |>
#     pivot_longer(`4/1/20`:`11/10/21`) |>
#     mutate(name = as.Date(name, "%m/%d/%y")) |>
#     group_by(name, Province_State) |>
#     summarize(state_cases = sum(value))

ggplot(jhu, aes(date, cases_avg, col=state, group=state)) +
    geom_line(alpha = 0.3)
"""
@rget jhu

begin
    jhu = stack(jhu, ["4/1/20", "11/10/21"], ["Admin2", "Province_State"]; variable_name=:Date, value_name=:Cases)
    fmt = dateformat"m/d/y"
    transform!(jhu, :Date => ByRow(x->Date(x, fmt)) => :Date)
    combine(groupby(jhu, [:Province_State, :Date]), :Cases => sum)
end

plot(jhu_ts.name, jhu_ts.state_cases, group=jhu_ts.Province_State, lab="")

mean_bias_obs(n, x; g=2/100, s=2/3, b=exp(1)) = log(b, n) / (g + s*x + 1) + (n - log(b, n)) * x

let n=5_000
    μs = [[a, b, mean_bias_obs(a, b)] for a in 1:n for b in [0.05, 0.2, 0.5]]
    pdat = reduce(hcat, μs)'
    plot(pdat[:,1], pdat[:,3], group=pdat[:,2])
    plot!(pdat[:,1], pdat[:,1], linealpha=0.5, linecolor=:gray)
    plot!(pdat[:,1], pdat[:,1] * 0.05, linealpha=0.5, linecolor=:gray, ylims=(0, 500))
end