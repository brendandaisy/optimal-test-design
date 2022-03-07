using Plots

f(n, x, y) = (n * x)^y * exp(-n * x)
f(n, y; xs=0:0.02:1) = sum(x->f(n, x, y), xs)

fnoy(n, x) = f(n, x, n*x)
fnoy(n; xs=0:0.02:1) = sum(x->fnoy(n, x), xs)

xs = 0:0.02:1
# plot(xs, [fnoy(10, a) / fnoy(10) for a ∈ xs])
# plot!(xs, [fnoy(100, a) / fnoy(100) for a ∈ xs])
# plot!(xs, [fnoy(500, a) / fnoy(500) for a ∈ xs])

plot(xs, [f(10, a, 4) / f(10, 4) for a ∈ xs], lab="n=10, y=4")
plot!(xs, [f(10, a, 7) / f(10, 7) for a ∈ xs], lab="n=10, y=7", ls=:dash)
plot!(xs, [f(100, a, 35) / f(100, 35) for a ∈ xs], lab="n=100, y=35")
plot!(xs, [f(100, a, 25) / f(100, 25) for a ∈ xs], lab="n=100, y=25", ls=:dash)
plot!(xs, [f(500, a, 100) / f(500, 100) for a ∈ xs], lab="n=500, y=100")
plot!(xs, [f(500, a, 25) / f(500, 25) for a ∈ xs], lab="n=500, y=25", ls=:dash)

savefig("poisson-obs.png")