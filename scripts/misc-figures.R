library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(grid)
library(gridExtra)


samp_rnot <- read_csv("data/sims/cond-samples/samp-rep-number.csv")

p1 <- ggplot(samp_rnot, aes(alpha, beta)) +
    geom_density_2d_filled(adjust=1.2) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\alpha$", y="$\\beta$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

p2 <- ggplot(samp_rnot, aes(beta, S0)) +
    geom_density_2d_filled(adjust=1.5) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\beta$", y="$S_0$", fill="Density") +
    theme(legend.position="none")

# Outbreak size-------------------------------------------------------------------
samp_osize <- read_csv("data/sims/cond-samples/samp-outbreak-size.csv")

p3 <- ggplot(samp_osize, aes(alpha, beta)) +
    geom_density_2d_filled(adjust=2) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\alpha$", y="$\\beta$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

p4 <- ggplot(samp_osize, aes(beta, S0)) +
    geom_density_2d_filled(adjust=1.5) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\beta$", y="$S_0$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

# Peak intensity------------------------------------------------------------------
samp_imax <- read_csv("data/sims/cond-samples/samp-peak-intensity.csv")

p5 <- ggplot(samp_imax, aes(alpha, beta)) +
    geom_density_2d_filled(adjust=1.6) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\alpha$", y="$\\beta$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

p6 <- ggplot(samp_imax, aes(beta, S0)) +
    geom_density_2d_filled(adjust=1.6) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\beta$", y="$S_0$", fill="Density") +
    theme(legend.position="none")

# Assemble plot-------------------------------------------------------------------

row_titles <- str_c("\\textbf{\\large", c("Reproductive number", "Outbreak size", "Peak intensity"), "}")
plots <- list(p1, p2, p3, p4, p5, p6)
# Add row titles
plots[c(1, 3, 5)] <- map2(plots[c(1, 3, 5)], row_titles, ~arrangeGrob(.x, left=.y))

tikz_plot(grid.arrange(grobs=plots, nrow=3), "cond-samples.pdf", 5.5, 8.5)

# TODO: for method paper, show "number of fit ODEs" for both RMD and variance of an estimator using MC+Maximum likelihood. Between fitting for each gradient (and possibly being more biased) the MC method may take more time
