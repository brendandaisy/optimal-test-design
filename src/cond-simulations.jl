using DrWatson
@quickactivate "optimal-dest-design"
using DEParamDistributions
using NamedTupleTools

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
    config = (;θfix, θprior, M)
    cond_sims, _ = produce_or_load(
        datadir("sims", subdir), config; 
        filename=mysavename(config; ignores=()), tag=false
    ) do c
        pdist = SIRParamDistribution(;c.θprior...)
        cond_simulations(c.θfix, pdist, c.M; dekwargs...)
    end
    return Dict(str2set(k) => v for (k, v) ∈ cond_sims)
end