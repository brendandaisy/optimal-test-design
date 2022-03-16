### To be included in scripts

const single_obs_dict = Dict(
    "normal" => (t, x; σ, n) -> Normal(n * x[t], σ),
    "normal_inf_noise" => (t, x; σ, n) -> Normal(n * x[t], σ * x[t]),
    "poisson" => (t, x; n) -> Poisson(n * x[t]),
    "neg_binom" => (t, x; r, n) -> NegativeBinomial(r, r / (r + n*x[t]))
    # "poisson_bias" => (t, x; b, n, g) -> Poisson(b*n*(x[t]/(g+x[t])(1- b)*n * x[t]),
)

const inct_dict = Dict(
    "normal" => (maxt, x; σ, n) -> Normal(n * x[1:maxt], σ),
    "normal_inf_noise" => (maxt, x; σ, n) -> Normal(n * x[1:maxt], σ * x[1:maxt]),
    "poisson" => (maxt, x; n) -> joint_poisson(n, x[1:maxt]),
    "neg_binom" => (maxt, x; r, n) -> joint_neg_binom(r, n .* x[1:maxt])
)