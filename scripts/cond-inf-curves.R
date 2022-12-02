library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)

read_cond_sim <- function(varname) {
    read_csv(paste0("data/sims/cond-simulations/sim-", varname, ".csv")) |> 
        mutate(id=1:n()) |> 
        pivot_longer(-id, names_to="t", values_to="inf") |> 
        mutate(t=as.double(t), var=varname)
}
    
theme_inf <- theme_half_open() +
    theme(
        # strip.background=element_blank(),
        # strip.text=element_blank(),
        axis.title=element_blank(),
        text=element_text(face="plain"),
        # # axis.title.x=element_text(size=rel(0.5)),
        # axis.text=element_blank(),
        # panel.grid=element_blank(),
        legend.position="none",
        plot.margin=unit(c(0.01, 0.01, 0.01, 0.01), "in")
    )

inf_plot_var <- function(var, idx) {
    sdf <- read_cond_sim(var)
    
    sdf <- filter(sdf, id %in% sample(length(unique(sdf$id)), 50)) |> 
        mutate(inf_lab="$I(t)$")
    
    sdf_summ <- sdf |> group_by(t) |> summarise(ymin=quantile(inf, 0.05), ymax=quantile(inf, 0.95))
    
    p <- ggplot(sdf, aes(t, inf)) +
        geom_ribbon(aes(t, ymin=ymin, ymax=ymax), data=sdf_summ, col="gray70", fill="gray70", inherit.aes=FALSE) +
        geom_line(aes(group=id), alpha=0.2, col="darkblue") +
        geom_line(data=true_inf, col="orange", size=1.2) +
        ylim(0, 0.65) +
        labs(x="Time $t$", y="$I(t)$", col=NULL) +
        theme_inf
        # geom_text(
        #     x=24, y=0.58, 
        #     label=str_extract(names(texlab)[which(texlab == var)], "\\$.+\\$"), 
        #     fontface="plain",
        #     size=5
        # )
    
    # if (idx %% 2 == 0)
    #     p <- p + theme(plot.background=element_rect(fill="lightblue", color="lightblue"))
    return(p)
}

inf_plot_var("peak-timing", 1)
