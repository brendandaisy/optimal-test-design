using LambertW

function get_var_transforms()
    Dict(
        "rep-number" => [rnot, inv_rnot, d_inv_rnot],
        "outbreak-size" => [outbreak_size, inv_outbreak_size, d_inv_outbreak_size],
        "peak-intensity" => [max_inf, inv_max_inf, d_inv_max_inf],
        "growth-rate" => [growth_rate, inv_growth_rate, d_inv_growth_rate]
    )
end

rnot(α, β) = β/α
inv_rnot(α, R) = α*R
d_inv_rnot(α, R) = α

reff(α, β, S₀) = β*S₀/α

function outbreak_size(α, β, S₀; I₀=0.01)
    sol = solve(SIRModel{Float64}(;stop=10_000, α, β, S₀, I₀); save_idxs=1)
    S₀ + I₀ - last(sol.u)
end

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
    if reff(α, β, S₀) <= 1
        return I₀
    end
    iR = α/β
    I₀ + S₀ - iR*log(S₀) - iR*(1-log(iR))
end

function inv_max_inf(α, β, imax; I₀=0.01)
    A = β/α
    B = exp(-1 - β/α*(imax-I₀))
    -1/A*lambertw(-B, -1)
end

function d_inv_max_inf(α, β, imax; I₀=0.01)
    # A = β/α
    B = exp(-1 - β/α*(imax-I₀))
    1 / (1 - B*exp(lambertw(-B, -1)))
end

growth_rate(α, β, S₀) = β*S₀ - α
inv_growth_rate(β, S₀, grate) = β*S₀ - grate
d_inv_growth_rate(β, S₀, grate) = 1
