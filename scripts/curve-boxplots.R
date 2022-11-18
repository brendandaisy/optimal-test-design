library(tidyverse)

inf_all <- bind_rows(
    read_cond_sim("α"),
    read_cond_sim("β"),
    read_cond_sim("S₀"),
    read_cond_sim("rep-number"),
    read_cond_sim("outbreak-size"),
    read_cond_sim("peak-intensity"),
    read_cond_sim("growth-rate")
)

envelope <- function(curves) {
    ts <- unique(curves$t)
    
    curveg <- curves |> 
        group_by(t)
    
    cmin <- slice_min(curveg, inf, n=1, with_ties=FALSE)
    cmax <- slice_max(curveg, inf, n=1, with_ties=FALSE)
    
    tibble(t=ts, ymin=cmin$inf, ymax=cmax$inf)
}

curve_in_env <- function(curve, cmin, cmax) all(cmin <= curve & curve <= cmax)

curve_scores <- function(curves, n_env, n_samples) {
    grcurves <- curves |> 
        mutate(score=0) |> 
        group_by(id)
    
    len <- n_groups(grcurves)
    
    for (i in 1:n_samples) {
        env <- curves |> 
            filter(id %in% sample(len, n_env)) |> 
            envelope()
        
        grcurves <- mutate(grcurves, score=score + curve_in_env(inf, env$ymin, env$ymax))
    }
    ungroup(grcurves)
}

#TODO: OK, let's give these settings a shot. Get envelope from each one with a reduce over varnames
tmp <- read_cond_sim("β")

cscores <- curve_scores(tmp, 50, 100)

cscores |> 
    slice_max(score, prop=0.75) |> 
    envelope() |> 
    ggplot(aes(t)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax))
    # geom_line(aes(t, inf, group=id), data=tmp2, alpha=0.25, inherit.aes=FALSE)
