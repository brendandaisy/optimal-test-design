## To be included in scripts

todaystr() = format(today(), "mm-dd")

mysavename(d) = replace(
    savename(d,  "jld2"; connector=" || ", allowedtypes=(String, Tuple, NamedTuple), ignores = (:dekwargs, :θprior)), 
    " = " => "="
)

function collect_csv(simdir=""; kwargs...)
    res = collect_results(datadir("sims", simdir); black_list=[:gitpatch, :script], kwargs...)
    CSV.write(datadir("sims", simdir, "results-$(todaystr()).csv"), string.(res))
end