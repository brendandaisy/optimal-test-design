library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)

read_cond_sim <- function(varname) {
    read_csv(paste0("data/sims/cond-simulations/sim-", varname, ".csv")) |> 
        mutate(id=1:n()) |> 
        pivot_longer(-id, names_to="t", values_to="inf") |> 
        mutate(t=as.double(t))
}

inf_rnot <- read_cond_sim("rep-number")
inf_osize <- read_cond_sim("outbreak-size")
inf_imax <- read_cond_sim("peak-intensity")
inf_grate <- read_cond_sim("growth-rate")

inf_all <- bind_rows(
    mutate(inf_rnot, func="Reproductive number"),
    mutate(inf_osize, func="Outbreak size"),
    mutate(inf_imax, func="Peak intensity"),
    mutate(inf_grate, func="Growth rate")
)
    
gg <- ggplot(inf_all, aes(t, inf, group=id, col=func)) +
    geom_line(alpha=0.15) +
    facet_wrap(~fct_inorder(func), ncol=1) +
    scale_x_continuous(guide = guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
    labs(x="$t$", y="$I(t)$") +
    theme_bw() +
    theme(legend.position="none")

ggsave("plots/cond-inf-curves.pdf", width=4, height=9)
