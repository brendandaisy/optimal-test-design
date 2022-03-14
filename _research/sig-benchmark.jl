## what is going on with slow scaling?

using DrWatson
@quickactivate "optimal-test-design"
using Distributions
using BenchmarkTools
using Statistics
using StatsFuns

@btime [Poisson(5) for _ in 1:10000];
@btime [Poisson.(fill(5, 10)) for _ in 1:10000];
@btime [product_distribution(Poisson.(fill(5, 10))) for _ in 1:10000];
@btime @inbounds [logpdf(product_distribution(Poisson.(fill(5, 10))), fill(5,10)) for _ in 1:10000];
@btime @inbounds [pdf(Poisson(5), 5) for _ in 1:10000];

dlist_naive(N=10_000, p=16) = [product_distribution(Poisson.(fill(5, p))) for _ in 1:N]

logy(λ, k) = k*log(λ) - log(factorial(k)) - λ
logyvec(λs, ks) = mapreduce((λ, k)->logy(λ, k), +, λs, ks)

me_naive(dlist, y) = map(d->pdf(d, y), dlist) |> mean |> log
me_statsfuns(dlist, y, M) = -log(M) + logsumexp(map(d->logpdf(d, y), dlist))
me_hardcode(λs, y, M) = -log(M) + logsumexp(map(x->logyvec(x, y), λs))

function me_mem_setup(M=10_000)

end

@benchmark dlist = dlist_naive(20_000)
λs = [fill(5, 16) for _ ∈ 1:20_000]
y = rand(4:7, length(dlist[1]))

@benchmark me_naive($dlist, $y)
@benchmark me_statsfuns($dlist, $y, $(length(dlist)))
@benchmark me_hardcode($λs, $y, $(length(λs)))

### Try to fix closure type inference (couldn't get it)
include(srcdir("observation-dicts.jl"))

function fake_util(xfunc)
    x = 1:10
    xfunc(x)
end

obs_func(maxt, x, mod; kw...) = DDD[mod](maxt, x; kw...)

function outer(mod; params...)
    i = 2
    fake_util(x->obs_func(i, x, mod; params...))
end

const DDD = inct_dict
obs_model = "neg_binom"

@code_warntype obs_func(2, [1, 2, 3], mod; r=10, n=10)
@code_warntype outer(obs_model; r=10, n=10)
@code_warntype fake_util(x->obs_func(2, x, mod; r=10, n=10))

