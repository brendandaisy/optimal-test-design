library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)

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
tsteps <- seq(1, length(true_inf)-1, by=1)

known_lab <- list(
    "Set([:α])" = "Unknown: beta, S0",
    "Set([:β])" = "Unknown: alpha, S0",
    "Set([:S₀])" = "Unknown: alpha, beta",
    "Set{Symbol}()" = "Unknown: alpha, beta, S0"
)

recover_sig <- function(df, var) {
    df |> mutate("{{var}}":=str_replace_all({{var}}, "Any|\\[|\\]", "")) |>
        separate({{var}}, str_c("u", as_label(enquo(var)), tsteps), sep=",", convert=TRUE, fill="right")
}

(marg_org <- read_csv("data/sims/increasing-tspan-marg/results-03-25.csv", na="missing"))

marg <- marg_org |>
    mutate(known=map_chr(known, ~known_lab[[.x]])) |>
    select(-path, true=θtrue, sigbeta=`sig-β`, sigalpha=`sig-α`, sigS0=`sig-S₀`)

wmarg <- marg |>
    recover_sig(sigbeta) |>
    recover_sig(sigalpha) |>
    recover_sig(sigS0) |>
    pivot_longer(contains("usig"), values_to="SIG") |>
    mutate(var=str_extract(name, "sig[a-zS]+"), t=as.double(str_extract(name, "\\d+")))

pdat <- wmarg |>
    drop_na() |>
    group_by(t, var, obs_params, known) |>
    summarise(SIG=mean(SIG)) |>
    ungroup() |>
    mutate(ntest=str_extract(obs_params, "n = \\d+")) |>
    filter(str_detect(known, "alp.*bet.*S0"), str_detect(ntest, "1000")) |>
    bind_rows(tibble(t=rep(0, 3), var=c("sigalpha", "sigbeta", "sigS"), SIG=c(0, 0, 0)))

texlab = c("$\\beta$" = "sigbeta", "$\\alpha$" = "sigalpha", "$S_0$" = "sigS")

gg <- pdat |>
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
        panel.grid.major.x = element_blank(),
        # panel.grid.minor.x = element_blank(),
        legend.position = c(.88, 0.15)
    ) +
    scale_color_manual(values=c("#d175ab", "#d175ab", "#d175ab"))
    scale_x_continuous(guide = guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
    scale_y_continuous(breaks=seq(0, 3, 0.5))

tikz_plot(gg, "main-fig", w=4.1, h=3.5)
