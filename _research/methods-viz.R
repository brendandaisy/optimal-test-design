library(tidyverse)
library(tikzDevice)

# options(tikzDocumentDeclaration = c(
#     "\\documentclass[10pt]{article}",
#     "\\fontfamily{cmss}\\selectfont"
# ))

tikz_plot <- function(ggp, fname = 'plots/tikz', w = 8.5, h = 4) {
  tikz(paste0(fname, '.tex'), standAlone = TRUE, width = w, height = h)
  print(ggp)
  dev.off()
  system(paste0('pdflatex ', fname, '.tex'))
}

theme_svg <- theme(
    panel.grid = element_blank(),
    axis.line = element_line(arrow=arrow(type="closed", angle=17, length=unit(0.13, "in")), size=2),
    panel.background = element_rect(fill="white"),
    axis.ticks = element_blank(),
    axis.text = element_blank()
)

p1dat <- tibble(
    qhi = qhi,
    qlo = qlo,
) |>
    bind_cols(curves)

gg <- p1dat |>
    pivot_longer(contains("x")) |>
    ggplot(aes(t, value, groups=name)) +
    geom_line(alpha=0.15, size=1.17) +
    theme_svg +
    labs(y="{\\fontfamily{cmss}\\selectfont\\Large prevalence}", x="{\\fontfamily{cmss}\\selectfont\\Large time}")

tikz_plot(gg, "pri-curve", w=4.5, h=4)

p2disc <- tibble(
    y = y[1:12],
    t = c6$t[seq(1, 111, 10)]
)

p2curves <- pivot_longer(c12, contains("x"))

gg <- ggplot(p2curves, aes(t, value, groups=name)) +
    geom_line(alpha=0.15, size=1.17) +
    geom_line(
        aes(t, xtrue), tibble(t=seq(0, 11, 0.1), xtrue=xtrue[seq(1, 111, 1)]), 
        col="#ffa24b", linetype="dashed", inherit.aes=FALSE, size=2.2
    ) +
    geom_point(aes(t, y / 100), p2disc, col="white", inherit.aes=FALSE, shape=21, fill="#d175ab", size=2.4) +
    theme_svg +
    labs(y="{\\fontfamily{cmss}\\selectfont\\Large prevalence}", x="{\\fontfamily{cmss}\\selectfont\\Large time}")

tikz_plot(gg, "plots/post-curve", w=4.5, h=4)