
θprior = (α=Uniform(0.05f0, 0.85f0), β=Uniform(0.3f0, 1.5f0), S₀=Uniform(0.1f0, 0.99f0))

tpeak_pri = map(_->peak_timing(rand(θprior.α), rand(θprior.β), rand(θprior.S₀)), 1:10_000)
CSV.write("_research/tpeak-pri.csv", DataFrame([tpeak_pri], [:V]))

post3 = CSV.read("_research/stan-samp-post3.csv", DataFrame)
tpeak_post3 = peak_timing.(post3.alpha, post3.beta, post3.S0)
CSV.write("_research/tpeak-post3.csv", DataFrame([tpeak_post3], [:V]))

post8 = CSV.read("_research/stan-samp-post8.csv", DataFrame)
tpeak_post8 = peak_timing.(post8.alpha, post8.beta, post8.S0)
CSV.write("_research/tpeak-post8.csv", DataFrame([tpeak_post8], [:V]))