using LambertW

rnot(α, β, S₀) = S₀*β/α
inv_rnot(α, β, R) = α*R/β
d_inv_rnot(α, β, R) = α/β

function outbreak_size(α, β, S₀; I₀=0.01)
    sol = solve(SIRModel{Float64}(;stop=10_000, α, β, S₀, I₀); save_idxs=1)
    S₀ + I₀ - last(sol.u)
end
# inv_rec_inf(α, β, r_inf; S₀, I₀=0.01) = (1 - r_inf) / exp(-β/α*(r_inf + S₀ + I₀ - 1))

function inv_outbreak_size(α, S₀, size; I₀=0.01)
    R₀ = 1 - S₀ - I₀
    -α / size * log((1-R₀-size) / S₀)
end

function d_inv_outbreak_size(α, S₀, size; I₀=0.01)
    R₀ = 1 - S₀ - I₀
    e1 = 1 - R₀ - size
    α / size * (1 / size * log(e1/S₀) + 1/e1)
end

function max_inf(α, β, S₀; I₀=0.01)
    iR = α/β
    I₀ + S₀ - iR*log(S₀) - iR*(1-log(iR))
end

function inv_max_inf(α, β, imax; I₀=0.01)
    A = β/α
    B = exp(-1 + log(α/β) - β/α*(imax-I₀))
    -1/A*lambertw(-A*B, -1)
end

function d_inv_max_inf(α, β, imax; I₀=0.01)
    A = β/α
    B = exp(-1 + log(α/β) - β/α*(imax-I₀))
    1 / (1 - A*B*exp(lambertw(-A*B, -1)))
end

growth_rate(α, β, S₀) = β*S₀ - α
inv_growth_rate(β, S₀, grate) = β*S₀ - grate
d_inv_growth_rate(β, S₀, grate) = 1
