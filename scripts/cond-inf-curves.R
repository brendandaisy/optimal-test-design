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

inf_all <- bind_rows(
    read_cond_sim("α"),
    read_cond_sim("β"),
    read_cond_sim("S₀"),
    read_cond_sim("rep-number"),
    read_cond_sim("outbreak-size"),
    read_cond_sim("peak-intensity"),
    read_cond_sim("growth-rate")
)
    # mutate(var=fct_relevel(fct_recode(var, !!!texlab), !!!names(texlab)))
    
theme_inf <- theme_bw() +
    theme(
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.title=element_blank(),
        # axis.title.x=element_text(size=rel(0.5)),
        # axis.text=element_text(size=rel(0.5)),
        panel.grid=element_blank(),
        legend.position="none",
        plot.margin=unit(c(0, 0, 0, 0.15), "in")
    )

inf_plot_var <- function(var, idx) {
    sdf <- read_cond_sim(var)
    
    sdf <- filter(sdf, id %in% sample(length(unique(sdf$id)), 50)) |> 
        mutate(inf_lab="$I(t)$")
    
    sdf_summ <- sdf |> group_by(t) |> summarise(ymin=quantile(inf, 0.05), ymax=quantile(inf, 0.95))
    
    p <- ggplot(sdf, aes(t, inf)) +
        geom_ribbon(aes(t, ymin=ymin, ymax=ymax), data=sdf_summ, col="gray70", inherit.aes=FALSE) +
        geom_line(aes(group=id, col=inf_lab), alpha=0.2) +
        # scale_x_continuous(guide = guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
        labs(x="$t$", y=NULL, col=NULL) +
        theme_inf
    
    if (idx %% 2 == 0)
        p <- p + theme(plot.background=element_rect(fill="lightblue", color="lightblue"))
    return(p)
}
    
vars <- c("α", "β", "S₀", "rep-number", "outbreak-size", "peak-intensity", "growth-rate")
inf_plot_var("rep-number")
plot_inf <- imap(vars, inf_plot_var)

ggsave("plots/cond-inf-curves.pdf", width=4, height=9)
