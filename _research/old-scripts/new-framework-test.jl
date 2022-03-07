using MonteCarloMeasurements
using OrdinaryDiffEq
using ModelParameters

# mutable struct ODEModel <: AbstractModel
#     u0
#     p
# end

# Base.parent(m::ODEModel) = [m.u0, m.p]
# ode_func(::Type{T}) where T <: ODEModel = error("$T requires an associated system of equations")
# ode_func(m::T) where T <: ODEModel = ode_func(T)

function sir!(dx, x, p, t)
    S, I = x
    β, α = p
    dx[1] = -β * I * S
    dx[2] = β * I * S - α * I
end

Base.@kwdef struct InitialValues
    S₀ = Param(0.99 ± 0, fixed=true)
    I₀ = Param(0.01 ± 0, fixed=true)
end

Base.@kwdef struct Parameters
    β = Param(0.3 ± 0, fixed=true)
    α = Param(0.1 ± 0, fixed=true)
end

function odeprob(m)
    u0, p = groupparams(m, :component)
    ODEProblem(sir!, stripparams.(u0), (0., 30.), stripparams.(p))
end

u0 = InitialValues(S₀=Param(1. ⊠ Beta(1, 1), fixed=false))
p = Parameters(β=Param(0.1..1.0, fixed=false))
sir_model = Model((u0, p))
prob = odeprob(sir_model)
solve(prob, Tsit5())