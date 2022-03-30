library(tidyverse)
library(ggthemes)
library(tikzDevice)

tikz_plot <- function(ggp, fname = 'tikz-plot', w = 8.5, h = 4, cd="plots") {
    cwd <- getwd()
    system(paste0("mkdir ", cd, "/tikz-work"))
    setwd(paste0(cd, "/tikz-work"))
    tryCatch({
        tikz(paste0(fname, '.tex'), standAlone = TRUE, width = w, height = h)
        print(ggp)
        }, finally=dev.off()
    )
    system(paste0('pdflatex ', fname, '.tex'))
    system("mv *.pdf ..")
    setwd(cwd)
}

theme_sig <- theme(
    panel.grid = element_blank(),
    axis.line = element_line(arrow=arrow(type="closed", angle=17, length=unit(0.13, "in")), size=2),
    panel.background = element_rect(fill="white"),
    axis.ticks = 
    axis.text = element_blank()
)

(resgrid <- read_csv("data/sims/dispersion-vs-tests/results-03-28.csv"))

pdat <- resgrid |>
    separate(obs_param_grid, c("r", "n"), ",") |>
    mutate(r=as.double(str_extract(r, "\\d+")), n=as.double(str_extract(n, "\\d+")))

gg <- ggplot(pdat, aes(r, n, z=SIG)) +
    # geom_raster(interpolate=TRUE) +
    # geom_contour(alpha=0.6, col="lightblue") +
    geom_contour_filled() +
    scale_fill_viridis_d(option = "magma") +
    scale_x_continuous(expand = expansion(), breaks=c(1,  seq(20, 100, 20))) +
    scale_y_continuous(expand = expansion()) +
    theme_bw() +
    # theme(
    #     legend.position = "bottom",
    #     axis.text = element_text(size=6)
    # ) +
    labs(
        x="{\\fontfamily{cmss}\\selectfont Dispersion $r$}", 
        y="{\\fontfamily{cmss}\\selectfont Testing rate $\\eta$}", 
        fill="{\\fontfamily{cmss}\\selectfont\\footnotesize Joint Information Gain}"
    )

tikz_plot(gg, "disp-vs-tests", w=5.4, h=3.85)

pdat |>
    filter(r == 100) |>
    ggplot(aes(n, SIG)) +
    geom_line() +
    stat_function(fun=log)
