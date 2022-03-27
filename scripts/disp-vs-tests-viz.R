library(tidyverse)
library(ggthemes)
library(tikzDevice)

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
    axis.text = element_blank()
)

(resgrid <- read_csv("_research/tmp/res.csv"))

pdat <- resgrid |>
    separate(obs_param_grid, c("r", "n"), ",") |>
    mutate(r=as.double(str_extract(r, "\\d+")), n=as.double(str_extract(n, "\\d+")))

ggplot(pdat, aes(r, n, fill=SIG, z=SIG)) +
    geom_raster(interpolate=TRUE) +
    geom_contour(alpha=0.6, col="gray80") +
    scale_fill_viridis_c(option = "magma") +
    scale_x_continuous(expand = expansion()) +
    scale_y_continuous(expand = expansion())
