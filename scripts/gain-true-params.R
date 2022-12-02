library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)
library(cowplot)

library(haven)
library(tidyverse)

yrs <- as.character(2010:2022)

read_sas(paste0("~/Downloads/prds_hosp10_yr2010.sas7bdat"))
walk(yrs, 
     ~write_csv(read_sas(paste0("~/Downloads/prds_hosp10_yr", .x, ".sas7bdat")), paste0("~/Downloads/prds_hosp10_yr", .x, ".csv")))

theme_fig2 <- function() {
    theme_bw() +
        theme(
            strip.background=element_blank(),
            strip.text=element_blank(),
            panel.spacing.x=unit(0.1, "cm")
        )
}

res <- read_csv("data/sims/gain-by-peak.csv") |> 
    rename(alpha=α, beta=β, S0=`S₀`) |> 
    mutate(
        Reff=beta*S0/alpha,
        grate=beta*S0 - alpha,
        tpeak=t,
        lab=fct_relevel(fct_recode(lab, !!!texlab), !!!names(texlab))
    ) |> 
    filter(tpeak < 30)
# bind_rows(
#     mutate(res, out=FALSE),
#     mutate(res_out, out=TRUE)
#     )|> 

cor_labs1 <- res |> 
    group_by(lab) |> 
    summarise(cor=round(cor(tpeak, md), 2), y=max(md*0.95)) |> 
    mutate(x=19, label=str_c("$\\rho=", cor, "$"))

res |> 
    group_by(lab) |> 
    summarise(cor=round(cor(grate, md), 2))

cor_labs2 <- res |> 
    group_by(lab) |> 
    summarise(cor=round(cor(Reff, md), 2), y=min(md)+0.2) |> 
    mutate(x=4.5, label=str_c("$\\rho=", cor, "$"))

p1 <- ggplot(res, aes(grate, md, col=beta)) +
    # geom_line(aes(group=as.factor(S0), col=S0), alpha=0.4) +
    geom_point() +
    geom_text(aes(x, y, label=label), data=cor_labs1, inherit.aes=FALSE, size=2.5) +
    facet_wrap(~lab, scales="free_y", nrow=1) +
    labs(x="$t_{\\textnormal{peak}}^*$", y="$\\delta_u$", col="$\\beta^*$") +
    theme_fig2()

p2 <- res |> 
    ggplot(aes(Reff, md, col=beta)) +
    # geom_line(aes(group=as.factor(S0), col=S0), alpha=0.4) +
    geom_point() +
    geom_text(aes(x, y, label=label), data=cor_labs2, inherit.aes=FALSE, size=2.5) +
    facet_wrap(~lab, scales="free_y", nrow=1) +
    labs(x="$\\mathcal{R}_e^*$", y="$\\delta_u$", col="$\\beta^*$") +
    theme_fig2()

gg <- plot_grid(p1, p2, align="h", nrow=2)

tikz_plot(gg, "tmp", 9.8, 3.8)

# TODO: rerun this one
filter(res, md < 2.5, lab %in% c("growth-rate", "peak-intensity"))
filter(res, tpeak > 24)

res |> 
    filter(Reff > 0.9) |> 
    ggplot(aes(beta, S0, fill=md)) +
    geom_tile() +
    stat_function(fun=~0.2/.x, col="tomato", xlim=c(0.2, 1.6)) +
    facet_wrap(~lab, nrow=2) +
    scale_fill_viridis_c()
    
ggsave("plots/gain-var-beta1.pdf", width=6.6, height=4.5)