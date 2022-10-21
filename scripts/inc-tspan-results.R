library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)
library(cowplot)

tikz_plot <- function(ggp, fname='tikz', w=8.5, h=4, dir="plots") {
  cur_wd <- getwd()
  setwd(paste0(cur_wd, "/", dir))
  tikz(paste0(fname, '.tex'), standAlone = TRUE, width = w, height = h)
  print(ggp)
  dev.off()
  system(paste0('pdflatex ', fname, '.tex'))
  setwd(cur_wd)
}

theme_md <- theme(
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
    "(:α)" = "$\\alpha$ known",
    "(:β)" = "$\\beta$ known",
    "(:S₀)" = "$S_0$ known",
    "()" = "All unknown"
)

recover_md <- function(df, var) {
  df |>
    mutate("{{var}}":=str_replace_all({{var}}, c("Any\\["="0,", "f"="e", "\\]"=""))) |>
    separate({{var}}, str_c("u", as_label(enquo(var)), tsteps), sep=",", convert=TRUE, fill="right")
}

# (res1 <- read_csv("_research/tmp/res.csv", na="missing"))
(res1 <- read_csv("data/sims/increasing-tspan-marg/results-10-02.csv", na="missing"))
(res2 <- read_csv("data/sims/increasing-tspan-marg/results-10-11.csv", na="missing"))

res <- res1 |> 
  rename(mdbeta=`md-β`, mdalpha=`md-α`, mdS0=`md-S₀`) |> 
  bind_cols(select(res2, mdR=`md-R`, mdjoint=`md-joint`))
    # mutate(known=map_chr(known, ~known_lab[[.x]])) |>
    # select(-path, true=θtrue, mdbeta=`md-β`, mdalpha=`md-α`, mdS0=`md-S₀`)

reduce(list(mdalpha, mdbeta, mdS0, mdR, mdjoint), recover_md, .init=res) |> 
  pivot_longer(contains("umd"), values_to="mdiv")

res_long <- res |>
  recover_md(mdalpha) |>
  recover_md(mdbeta) |>
  recover_md(mdS0) |>
  recover_md(mdR) |>
  recover_md(mdjoint) |>
  pivot_longer(contains("umd"), values_to="mdiv") |>
  mutate(var=str_extract(name, "md[a-zRS]+"), t=as.double(str_extract(name, "\\d+"))) |>
  mutate(ntest=str_extract(obs_mod, "\\d+")) |> 
  select(-obs_mod, -name)

texlab = c(
  "$\\beta$"="mdbeta", "$\\alpha$"="mdalpha", "$S_0$"="mdS", 
  "$\\mathcal{R}$"="mdR", "$(\\alpha, \\beta, S_0)$"="mdjoint"
)

gg <- res_long |>
    ggplot(aes(t, mdiv, col=fct_recode(var, !!!texlab)), group=var) +
  # ggplot(aes(t, mdiv, col=var, group=var)) +
    geom_vline(xintercept=peak, col="#ffa24b", linetype="dashed", size=1.3) +
    geom_line(size=1.5, alpha=0.7) +
    # facet_wrap(~known, nrow=1) +
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
    # scale_color_manual(values=c("#d175ab", "#78cf8c", "#7f82db")) +
    scale_x_continuous(guide = guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29)
    # scale_y_continuous(breaks=seq(0, 4, 0.5))

tikz_plot(gg, "increasing-tspan-new", w=4.2, h=3.6)
