using DrWatson
@quickactivate "optimal-dest-design"
using DEParamDistributions
using NamedTupleTools

function prob_αβ_cond_R(α, β; R, prior)
    pdf(prior.α, α) * pdf(prior.β, β) * α/β * pdf(prior.S₀, α*R/β)
end

function accept_reject(sampler, target; M=100, N=100)
    ret = []
    for _ in 1:N
        s = rand(sampler)
        u = rand()
        while u >= target(s) / (M*pdf(sampler, s))
            s = rand(sampler)
            u = rand()
        end
        push!(ret, s)
    end
    return (first.(ret), last.(ret))
end

α = Uniform(0.1, 0.4)
β = Uniform(0.3, 2.)

αcond, βcond = accept_reject(product_distribution([α, β]), x->joint_cond(x[1], x[2], 3.))
Scond = (αcond * 3) ./ βcond

(βcond .* Scond) ./ αcond

density(βcond)

tup2str(nt::Tuple) = replace(string(nt), r"\(|\)|\:|,\)" => "")

separate(nt::NamedTuple) = [(k, ) for k ∈ keys(nt)]

function cross(nt::NamedTuple)
    k = keys(nt)
    [(k[i], k[j]) for i=eachindex(k) for j=eachindex(k) if i < j]
end

function cond_simulations(θfix, pdist::AbstractDEParamDistribution, M=100; dekwargs...)
    ret = Dict(
        tup2str(keys(θfix)) => solve_de_problem(pdist, θfix; dekwargs...),
        "" => simulate(pdist, M; dekwargs...)
    )
    for s ∈ separate(θfix)
        ret[tup2str(s)] = simulate(NamedTupleTools.select(θfix, s), pdist, M; dekwargs...)
    end
    for c ∈ cross(θfix)
        ret[tup2str(c)] = simulate(NamedTupleTools.select(θfix, c), pdist, M; dekwargs...)
    end
    return ret
end

str2set(s) = length(s) == 0 ? Set{Symbol}() : Symbol.(split(s,", ")) |> Set
set2str(s) = join(map(string, vcat(s...)), ", ")

function get_cond_sims(θfix, θprior, M=100; subdir="cond-simulations", dekwargs...)
    config = (;θfix, θprior, M, dekwargs...)
    cond_sims, _ = produce_or_load(
        datadir("sims", subdir), config; 
        filename=mysavename(config; ignores=()), tag=false
    ) do c
        pdist = SIRParamDistribution(;c.θprior...)
        cond_simulations(c.θfix, pdist, c.M; dekwargs...)
    end
    return Dict(str2set(k) => v for (k, v) ∈ cond_sims)
end