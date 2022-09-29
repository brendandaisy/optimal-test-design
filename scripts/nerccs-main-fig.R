library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)
library(cowplot)

tikz_plot <- function(ggp, fname = 'plots/tikz', w = 8.5, h = 4) {
  tikz(paste0(fname, '.tex'), standAlone = TRUE, width = w, height = h)
  print(ggp)
  dev.off()
  system(paste0('pdflatex ', fname, '.tex'))
}

theme_sig <- theme(
    panel.grid = element_blank(),
    axis.line = element_line(arrow=arrow(type="closed", angle=17, length=unit(0.13, "in")), size=2),
    panel.background = element_rect(fill="white"),
    axis.ticks = 
    axis.text = element_blank(),
)

true_inf <- c(
    0.01, 0.017234986028274994, 0.029226916688214286, 0.0482619354074302, 0.07644676001950913, 0.11397801673090434, 0.15694203765403592,
    0.197068585361762, 0.22556384374277325, 0.23791293991026669, 0.2351898235396991, 0.22162345300496464, 0.2018696959400936, 
    0.1795447345733238, 0.15704515578805922, 0.13574849445435141, 0.11635024064541434, 0.09910377645846799, 0.08402199574081326, 0.0709903078468306, 0.059817465836561626, 0.05029848904490269, 0.04222657248401512, 0.03540215277236312, 
    0.029650387465228304, 0.024814410446680653, 0.02075294687536625, 0.017345265065778624, 0.014491436294653753, 0.012103619476697134, 0.010105741535796938
)
peak <- which.max(true_inf)
tsteps <- seq(0, length(true_inf)-1, by=1)

known_lab <- list(
    "Set([:α])" = "$\\alpha$ known",
    "Set([:β])" = "$\\beta$ known",
    "Set([:S₀])" = "$S_0$ known",
    "Set{Symbol}()" = "All unknown"
)

recover_sig <- function(df, var) {
    df |>
        mutate("{{var}}":=str_replace_all({{var}}, c("Any\\[" = "0,", "\\]" = ""))) |>
        separate({{var}}, str_c("u", as_label(enquo(var)), tsteps), sep=",", convert=TRUE, fill="right")
}

(marg_org <- read_csv("_research/tmp/res.csv", na="missing"))
(marg_org <- read_csv("data/sims/increasing-tspan-marg/results-04-28.csv", na="missing"))

marg <- marg_org |>
    mutate(known=map_chr(known, ~known_lab[[.x]])) |>
    select(-path, true=θtrue, sigbeta=`sig-β`, sigalpha=`sig-α`, sigS0=`sig-S₀`)

wmarg <- marg |>
    recover_sig(sigbeta) |>
    recover_sig(sigalpha) |>
    recover_sig(sigS0) |>
    pivot_longer(contains("usig"), values_to="SIG") |>
    mutate(var=str_extract(name, "sig[a-zS]+"), t=as.double(str_extract(name, "\\d+"))) |>
    drop_na() |>
    group_by(t, var, obs_model, known) |>
    summarise(SIG=mean(SIG)) |>
    ungroup() |>
    mutate(ntest=str_extract(obs_model, "\\d+"))

pdat <- wmarg |>
    filter(str_detect(ntest, "1000"), !str_detect(known, "alpha"))

texlab = c("$\\beta$" = "sigbeta", "$\\alpha$" = "sigalpha", "$S_0$" = "sigS")

gg <- pdat |>
    ggplot(aes(t, SIG, col=fct_recode(var, !!!texlab)), group=var) +
    geom_vline(xintercept=peak, col="#ffa24b", linetype="dashed", size=1.3) +
    geom_line(size=1.5, alpha=0.7) +
    facet_wrap(~known, nrow=1) +
    labs(
        x="{\\fontfamily{cmss}\\selectfont\\large Days of observation}", 
        y="{\\fontfamily{cmss}\\selectfont\\large Marginal Information Gain}",
        col="{\\fontfamily{cmss}\\selectfont Parameter}"
    ) +
    theme_bw() +
    theme(
        # legend.background = element_blank(),
        # legend.key = element_blank(),
        # panel.grid.major.x = element_blank(),
        strip.text = element_text(size=9, margin=margin(2, 0, 2, 0)),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom"
    ) +
    scale_color_manual(values=c("#d175ab", "#78cf8c", "#7f82db")) +
    scale_x_continuous(guide = guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
    scale_y_continuous(breaks=seq(0, 4, 0.5))

tikz_plot(gg, "main-fig", w=7.2, h=3.6)

### Learning stages chart

###

pdat <- wmarg |>
    filter(str_detect(known, "alp.*bet.*S0"), str_detect(ntest, "1000")) |>
    bind_rows(tibble(t=rep(0, 3), var=c("sigalpha", "sigbeta", "sigS"), SIG=c(0, 0, 0)))

inset <- wmarg |>
    filter(known %in% c("Unknown: alpha, S0", "Unknown: alpha, beta"), str_detect(ntest, "1000")) |>
    bind_rows(tibble(
        t=rep(0, 4), 
        known=c(rep("Unknown: alpha, S0", 2), rep("Unknown: alpha, beta", 2)),
        var=c("sigalpha", "sigS", "sigalpha", "sigbeta"), 
        SIG=rep(0, 4)
    )) |>
    mutate(known=fct_recode(known, "$\\beta^*$ known"="Unknown: alpha, S0", "$S_0^*$ known"="Unknown: alpha, beta"))

ggmain <- pdat |>
    ggplot(aes(t, SIG, col=fct_recode(var, !!!texlab)), group=var) +
    geom_vline(xintercept=peak, col="#ffa24b", linetype="dashed", size=1.5) +
    geom_line(size=1.5, alpha=0.7) +
    labs(
        x="{\\fontfamily{cmss}\\selectfont\\large Days of observation}", 
        y="{\\fontfamily{cmss}\\selectfont\\large Marginal Information Gain}",
        col="{\\fontfamily{cmss}\\selectfont\\footnotesize Parameter}"
    ) +
    # xlim(0, 33) +
    theme_bw() +
    theme(
        legend.background = element_blank(),
        legend.key = element_blank(),
        # panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(1.2, 0.18)
    ) +
    scale_color_manual(values=c("#d175ab", "#78cf8c", "#7f82db")) +
    scale_x_continuous(guide = guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
    scale_y_continuous(breaks=seq(0, 3, 0.5))

ggin <- ggplot(inset, aes(t, SIG, col=fct_recode(var, !!!texlab)), group=var) +
    geom_vline(xintercept=peak, col="#ffa24b", linetype="dashed", size=1.5) +
    geom_line(size=1.1, alpha=0.7) +
    labs(
        x="{\\fontfamily{cmss}\\selectfont\\large Days}", 
        y="{\\fontfamily{cmss}\\selectfont\\large MIG}",
        col="{\\fontfamily{cmss}\\selectfont\\footnotesize Parameter}"
    ) +
    # xlim(0, 33) +
    facet_wrap(~known, nrow=1) +
    theme_bw() +
    theme(
        text = element_text(size=3),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        # panel.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=5),
        # strip.background = element_rect(size=5)
    ) +
    scale_color_manual(values=c("#d175ab", "#78cf8c", "#7f82db"))
    # scale_x_continuous(breaks=seq(0, 30, 10)) +
    # scale_y_continuous(breaks=seq(0, 3, 0.5))

gg <- ggdraw() + draw_plot(ggmain) + draw_plot(ggin, x=0.12, y=0.64, width=0.45, height=.3)

tikz_plot(gg, "main-fig", w=4.5, h=3.5)
