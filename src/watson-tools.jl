## To be included in scripts

todaystr() = format(today(), "mm-dd")

mysavename(d) = replace(
    savename(d,  "jld2"; connector=" || ", allowedtypes=(String, Tuple, NamedTuple), ignores = (:dekwargs, :θprior)), 
    " = " => "="
)