### To be included in scripts

single_obs_dict = Dict(
    "normal" => (t, x; σ, n) -> Normal(n * x[t], σ),
    "normal_inf_noise" => (t, x; σ, n) -> Normal(n * x[t], σ * x[t]),
    "poisson" => (t, x; n) -> Poisson(n * x[t]),
    "neg_binom" => (t, x; r, n) -> NegativeBinomial(r, r / (r + n*x[t]))
)