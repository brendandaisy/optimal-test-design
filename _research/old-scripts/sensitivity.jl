using OrdinaryDiffEq
using ForwardDiff
using DiffEqSensitivity
using Distributions
using Plots

function newsir!(dx, x, p, t) # fact that pop size doesn't matter as long put Beta dist on S_0
    S, I = x
    β, α = p
    dx[1] = -β * I * S
    dx[2] = β * I * S - α * I
end

function nℓlik(y, n, t1, t2, x)
    d = Poisson.([x[t1]*n[1], x[t2]*n[2]])
    -logpdf(product_distribution(d), y)
end

θ = [0.7, 0.1]
x₀ = [0.7, 0.01]
tspan = (0., 60.)
prob = ODEProblem(newsir!, x₀, tspan, θ)
xtrue = solve(prob; saveat=1., save_idxs=2)
plot(xtrue)

t1 = 7
t2 = 10
n = [50, 50]
y = rand.(Poisson.([xtrue[t1]*n[1], xtrue[t2]*n[2]]))

function _nℓlik(θ)
    x = solve(prob, Vern9(); p=[θ[2], 0.1], u0=[θ[1], 0.01], saveat=1., save_idxs=2, abstol=1e-12, reltol=1e-12)
    nℓlik(y, n, t1, t2, x)
end
H = ForwardDiff.hessian(_nℓlik, [0.7, 0.7])
Σ = inv(H)

P̂ = MvNormal([0.7, 0.7], Σ)

# H  = second_order_sensitivities(x->ℓlik(y, n, t1, t2, x), prob, Vern9(), saveat=1., abstol=1e-12, reltol=1e-12)