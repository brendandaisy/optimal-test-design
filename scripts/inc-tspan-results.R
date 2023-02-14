library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)
library(cowplot)
library(grid)
library(gridExtra)

tikz_plot <- function(ggp, fname='tikz', w=8.5, h=4, dir="plots") {
  cur_wd <- getwd()
  setwd(paste0(cur_wd, "/", dir))
  tikz(paste0(fname, '.tex'), standAlone = TRUE, width = w, height = h)
  print(ggp)
  dev.off()
  system(paste0('lualatex ', fname, '.tex'))
  setwd(cur_wd)
}

true_inf <- read_csv("data/sims/inf-true.csv") |> 
  mutate(t=seq(0, 30, 0.2))
peak <- which.max(true_inf$inf)
tsteps <- seq(0, 30, 1)

recover_md <- function(df, var) {
  df |>
    mutate("{{var}}":=str_replace_all({{var}}, c("Any\\["="0,", "f"="e", "\\]"=""))) |> 
    separate({{var}}, str_c("u", as_label(enquo(var)), tsteps), sep=",", convert=TRUE, fill="right")
}

(res <- read_csv("data/sims/increasing-tspan-marg/results-11-30.csv", na="missing"))
# (res <- read_csv("data/sims/inc-tspan-slow.csv", na="missing"))

res <- mutate(res[2, ], `md-peak-timing`=res[1, "md-peak-timing"])

texlab = c(
  "Recovery rate $\\alpha$"="α", 
  "Infection rate $\\beta$"="β", 
  "Initial susceptible $S_0$"="S₀", 
  "Basic rep. number $\\mathcal{R}$"="rep-number", 
  "Outbreak size $\\mathcal{O}$"="outbreak-size",
  "Peak intensity $\\mathcal{P}$"="peak-intensity",
  "Peak timing $\\mathcal{T}$"="peak-timing",
  "Growth rate $\\mathcal{G}$"="growth-rate"
)

res_long <- res |>
  recover_md(`md-α`) |>
  recover_md(`md-β`) |>
  recover_md(`md-S₀`) |>
  recover_md(`md-rep-number`) |>
  recover_md(`md-outbreak-size`) |>
  recover_md(`md-peak-intensity`) |>
  recover_md(`md-peak-timing`) |>
  recover_md(`md-growth-rate`) |>
  pivot_longer(contains("umd"), values_to="mdiv") |>
  mutate(var=str_extract(name, "(?<=md-)[^\\d]+"), t=as.double(str_extract(name, "\\d+"))) |>
  select(mdiv, var, t) |> 
  mutate(var=fct_relevel(fct_recode(var, !!!texlab), !!!names(texlab)))

md_plot_var <- function(var, idx) {
  p <- res_long |>
    filter(var == !!var) |> 
    ggplot(aes(t, mdiv)) +
    # geom_vline(xintercept=peak, col="#ffa24b", linetype="dashed", size=1.3) +
    geom_segment(x=3, xend=3, y=-0.1, yend=1.1, col="lightblue", size=0.9, alpha=0.5, linetype="dashed") +
    geom_segment(x=8, xend=8, y=-0.1, yend=1.1, col="#f53db5", size=0.9, alpha=0.5, linetype="dashed") +
    geom_line(size=1.6, col="gray30") +
    labs(title=var, x=NULL) +
    ylim(0, 4.5) +
    # scale_x_continuous(guide=guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
    theme_bw() +
    theme(
      strip.background=element_blank(),
      axis.title=element_blank(),
      axis.text=element_text(size=rel(1.1)),
      panel.grid.minor.x=element_blank(),
      legend.position="none",
      plot.margin = unit(c(0.03, 0, 0, 0.1), "in")
    )
  
  # if (idx %% 2 == 0)
  #   p <- p + theme(plot.background=element_rect(fill="lightblue", color="lightblue"))
  return(p)
}

plot_md <- imap(names(texlab), md_plot_var)

md_rows <- map(1:8, ~{
  gdraw <- ggdraw(plot_md[[.x]])
  if (.x < 6)
    # gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.09, y=0.88, width=0.49, height=0.39, vjust=1)
    gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.16, y=0.88, width=0.47, height=0.37, vjust=1)
  else
    gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.96, y=0.1, width=0.47, height=0.37, hjust=1)
    # gdraw <- gdraw + draw_plot(plot_dens[[.x]], x=0.9, y=0.1, width=0.49, height=0.39, hjust=1)
  gdraw
})

panel1 <- plot_grid(plotlist=md_rows, align="h", axis="l", nrow=2)
y_lab <- textGrob("Identifiability $\\delta_u$", rot=90, gp=gpar(fontsize=16))
x_lab <- textGrob("Days of observation", gp=gpar(fontsize=16))
panel1 <- grid.arrange(arrangeGrob(panel1, left=y_lab, bottom=x_lab))

##

inf_pri <- read_csv("data/sims/sim-full-prior.csv") |> 
  mutate(id=1:n()) |> 
  pivot_longer(-id, names_to="t", values_to="inf") |> 
  mutate(t=as.double(t))

inf_pri_summ <- inf_pri |> 
  group_by(t) |> 
  summarise(ymin=quantile(inf, 0.025), ymax=quantile(inf, 0.975))

panel2 <- 
  # inf_pri |> 
  # filter(id %in% sample(max(id), 50)) |> 
  # ggplot(aes(t, inf)) +
  # geom_ribbon(
  #   aes(t, ymin=ymin, ymax=ymax), data=inf_pri_summ, 
  #   alpha=0.8, col="gray80", fill="gray80", inherit.aes=FALSE
  # ) +
  ggplot(true_inf, aes(t, inf)) +
  # stat_function(fun=~0.01 * exp(.x * (beta_true*Strue - alpha_true)), col="gray30", alpha=0.8, xlim=c(0, 6)) +
  # geom_line(aes(group=id), alpha=0.2, col="gray30") +
  annotate("rect", xmin=-1, ymin=0, xmax=6, ymax=0.5, fill="lightblue", alpha=0.35) +
  annotate("rect", xmin=6, ymin=0, xmax=14, ymax=0.5, fill="#f53db5", alpha=0.2) +
  annotate("rect", xmin=14, ymin=0, xmax=31, ymax=0.5, fill="blue4", alpha=0.2) +
  geom_line(col="orange", size=1.4) +
  geom_point(data=tibble(t=0:30, inf=stan_dat$y/1000), col="#f53db5", size=0.9) +
  scale_x_continuous(guide=guide_prism_minor(), breaks=seq(0, 30, 5), minor_breaks=1:29) +
  labs(x="Days $t$", y="$I(t)$", col=NULL) +
  theme_half_open() +
  theme(plot.margin=unit(c(0, 0, 0, 0.25), "in"))

gg <- plot_grid(
  panel1, 
  plot_grid(NULL, panel2, NULL, ncol=1, rel_heights=c(0.25, 0.5, 0.25), labels=c("", "B", "")),
  rel_widths=c(0.67, 0.33),
  labels=c("A", "")
)

tikz_plot(plot_grid(panel1), "increasing-tspan-p1", w=7.95, h=5.42)
 
##
vars <- c("α", "β", "S₀", "rep-number", "outbreak-size", "peak-intensity", "growth-rate")
pi <- imap(vars, inf_plot_var)
pi_top <- plot_grid(plotlist=pi[1:4], nrow=1)
pi_bottom <- plot_grid(NULL, pi[[5]], NULL, pi[[6]], NULL, pi[[7]], NULL, nrow=1, rel_widths=c(1/4, 1, 1/4, 1, 1/4, 1, 1/4))
panel3 <- plot_grid(pi_top, pi_bottom, ncol=1)

panel3 <- plot_grid(
  NULL, plot_grid(plotlist=pi, nrow=1), NULL, 
  nrow=3, rel_heights=c(0.15, 0.7, 0.15)
)

##
plot_grid(panel1, NULL, panel2, panel3, nrow=2, rel_widths=c(0.99, 0.01, 0.4, 0.6))

gg <- plot_grid(
  panel1, NULL,
  plot_grid(panel2, panel3, NULL, align="v", axis="b", nrow=1, rel_widths=c(0.5, 1, 0.2), labels=c("B", "C", "")),
  nrow=3,
  rel_heights=c(1, 0.08, 0.92),
  labels=c("A", ""), vjust=0
)

gg <- plot_grid(
  NULL,
  panel1, 
  NULL,
  plot_grid(panel2, panel3, nrow=1, labels=c("B", "C", ""), rel_widths=c(0.6, 1)),
  nrow=4,
  rel_heights=c(0.1, 1, 0.08, 0.92),
  labels=c("A", "")
)

tikz_plot(gg, "increasing-tspan", w=15.2, h=6)
