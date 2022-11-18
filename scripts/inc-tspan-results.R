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
  system(paste0('lualatex ', fname, '.tex'))
  setwd(cur_wd)
}

true_inf <- c(
    0.01, 0.017234986028274994, 0.029226916688214286, 0.0482619354074302, 0.07644676001950913, 0.11397801673090434, 0.15694203765403592,
    0.197068585361762, 0.22556384374277325, 0.23791293991026669, 0.2351898235396991, 0.22162345300496464, 0.2018696959400936, 
    0.1795447345733238, 0.15704515578805922, 0.13574849445435141, 0.11635024064541434, 0.09910377645846799, 0.08402199574081326, 0.0709903078468306, 0.059817465836561626, 0.05029848904490269, 0.04222657248401512, 0.03540215277236312, 
    0.029650387465228304, 0.024814410446680653, 0.02075294687536625, 0.017345265065778624, 0.014491436294653753, 0.012103619476697134, 0.010105741535796938
)
peak <- which.max(true_inf)
tsteps <- seq(0, length(true_inf)-1, by=1)

recover_md <- function(df, var) {
  df |>
    mutate("{{var}}":=str_replace_all({{var}}, c("Any\\["="0,", "f"="e", "\\]"=""))) |> 
    separate({{var}}, str_c("u", as_label(enquo(var)), tsteps), sep=",", convert=TRUE, fill="right")
}

# (res1 <- read_csv("_research/tmp/res.csv", na="missing"))
(res <- read_csv("data/sims/increasing-tspan-marg/results-11-10.csv", na="missing"))
(res2 <- read_csv("data/sims/results-11-10.csv", na="missing"))

# res <- res |> 
#   rename(mdbeta=`md-β`, mdalpha=`md-α`, mdS0=`md-S₀`)
  # bind_cols(
  #   select(
  #     res2, mdR=`rep-number`, mdosize=`outbreak-size`, mdimax=`peak-intensity`, mdgrate=`growth-rate`
  #     )
  #   )
    # mutate(known=map_chr(known, ~known_lab[[.x]])) |>
    # select(-path, true=θtrue, mdbeta=`md-β`, mdalpha=`md-α`, mdS0=`md-S₀`)

res_long <- res |>
  recover_md(`md-α`) |>
  recover_md(`md-β`) |>
  recover_md(`md-S₀`) |>
  recover_md(`md-rep-number`) |>
  recover_md(`md-outbreak-size`) |>
  recover_md(`md-peak-intensity`) |>
  recover_md(`md-growth-rate`) |>
  pivot_longer(contains("umd"), values_to="mdiv") |>
  mutate(var=str_extract(name, "(?<=md-)[^\\d]+"), t=as.double(str_extract(name, "\\d+"))) |>
  select(mdiv, var, t) |> 
  mutate(var=fct_relevel(fct_recode(var, !!!texlab), !!!names(texlab)))

texlab = c(
  "Recovery rate $\\alpha$"="α", 
  "Infection rate $\\beta$"="β", 
  "Initial susceptible $S_0$"="S₀", 
  "Basic rep. number $\\mathcal{R}$"="rep-number", 
  "Outbreak size $\\mathcal{O}$"="outbreak-size",
  "Peak intensity $\\mathcal{P}$"="peak-intensity",
  "Growth rate $\\mathcal{G}$"="growth-rate"
)

md_plot_var <- function(var, idx) {
  p <- res_long |>
    filter(var == !!var) |> 
    ggplot(aes(t, mdiv)) +
    # geom_vline(xintercept=peak, col="#ffa24b", linetype="dashed", size=1.3) +
    geom_segment(x=3, xend=3, y=-0.1, yend=0.5, col="#f53db5", size=1.3, alpha=0.7, linetype="dashed") +
    geom_segment(x=8, xend=8, y=-0.1, yend=0.5, col="purple", size=1.3, alpha=0.7, linetype="dashed") +
    geom_line(size=1.5, col="gray30") +
    labs(title=var, x=NULL) +
    ylim(0, 4.5) +
    scale_x_continuous(guide=guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
    theme_bw() +
    theme(
      strip.background=element_blank(),
      axis.title=element_blank(),
      axis.text=element_text(size=rel(1.1)),
      panel.grid.minor.x=element_blank(),
      legend.position="none",
      plot.margin = unit(c(0.02, 0.15, 0, 0.02), "in")
    )
  
  if (idx %% 2 == 0)
    p <- p + theme(plot.background=element_rect(fill="lightblue", color="lightblue"))
  return(p)
}

plot_md <- imap(names(texlab), md_plot_var)
plot_cols <- map(1:7, ~{
  bot <- plot_grid(plot_dens[[.x]], plot_inf[[.x]], align="h", rel_widths=c(1, 1.2))
  ret <- plot_grid(plot_md[[.x]], bot, align="hv", axis="r", nrow=2, rel_heights=c(1, 0.6))
  if (.x %% 2 == 0) {
    ret <- ret + theme(
      plot.background=element_rect(fill="lightblue", color="lightblue")
      # panel.background=element_rect(fill="lightblue", color="lightblue")
      )
  }
  ret
})

md_rows <- map(1:7, ~{
  gdraw <- ggdraw(plot_md[[.x]])
  if (.x < 6)
    gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.09, y=0.88, width=0.5, height=0.3, vjust=1)
  else
    gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.88, y=0.08, width=0.5, height=0.3, hjust=1)
  gdraw
})

gg <- plot_grid(plotlist=md_rows, align="h", axis="l", nrow=1)

gg <- plot_grid(plotlist=plot_cols, align="h", axis="t", nrow=1)
gg <- plot_grid(plot_md, bot_row, nrow=2, align="v", axis="r", rel_heights=c(1, 0.75))
tikz_plot(gg, "md-rows", w=15.2, h=3)
