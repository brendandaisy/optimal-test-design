## To be included in scripts
using Dates
using CSV

todaystr() = Dates.format(today(), "mm-dd")

function mysavename(d; ignores=(:dekwargs,)) 
    replace(
        savename(
            d,  "jld2"; 
            connector=" || ", 
            allowedtypes=(Any,),
            ignores
        ), 
    " = " => "=",
    " ± " => "±",
    r"\{\w+\d+\}" => ""
    )
end

function collect_csv(simdir=""; kwargs...)
    res = collect_results(datadir("sims", simdir); black_list=[:gitpatch, :script], kwargs...)
    CSV.write(datadir("sims", simdir, "results-$(todaystr()).csv"), string.(res))
end