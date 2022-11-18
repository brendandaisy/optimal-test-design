library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(grid)
library(gridExtra)

samp_rnot <- read_csv("data/sims/cond-samples/samp-rep-number.csv")

p1 <- ggplot(samp_rnot, aes(alpha, beta)) +
    geom_density_2d_filled(adjust=1) +
    geom_point(data=slice_sample(samp_rnot, n=100), col="gray40", shape=1) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\alpha$", y="$\\beta$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

p2 <- ggplot(samp_rnot, aes(beta, S0)) +
    geom_density_2d_filled(adjust=2) +
    geom_point(data=slice_sample(samp_rnot, n=100), col="gray40", shape=1) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\beta$", y="$S_0$", fill="Density") +
    theme(legend.position="none")

# Outbreak size-------------------------------------------------------------------
samp_osize <- read_csv("data/sims/cond-samples/samp-outbreak-size.csv")

p3 <- ggplot(samp_osize, aes(alpha, beta)) +
    geom_density_2d_filled(adjust=2) +
    geom_point(data=slice_sample(samp_osize, n=100), col="gray40", shape=1) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\alpha$", y="$\\beta$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

p4 <- ggplot(samp_osize, aes(beta, S0)) +
    geom_density_2d_filled(adjust=1.5) +
    geom_point(data=slice_sample(samp_osize, n=100), col="gray40", shape=1) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\beta$", y="$S_0$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

# Peak intensity------------------------------------------------------------------
samp_imax <- read_csv("data/sims/cond-samples/samp-peak-intensity.csv")

p5 <- ggplot(samp_imax, aes(alpha, beta)) +
    geom_density_2d_filled(adjust=1.4) +
    geom_point(data=slice_sample(samp_imax, n=100), col="gray40", shape=1) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\alpha$", y="$\\beta$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

p6 <- ggplot(samp_imax, aes(beta, S0)) +
    geom_density_2d_filled(adjust=2) +
    geom_point(data=slice_sample(samp_imax, n=100), col="gray40", shape=1) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\beta$", y="$S_0$", fill="Density") +
    theme(legend.position="none")

# Peak intensity------------------------------------------------------------------
samp_grate <- read_csv("data/sims/cond-samples/samp-growth-rate.csv")

p7 <- ggplot(samp_grate, aes(alpha, beta)) +
    geom_density_2d_filled(adjust=1.6) +
    geom_point(data=slice_sample(samp_grate, n=100), col="gray40", shape=1) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\alpha$", y="$\\beta$", fill="Density") +
    theme_bw() +
    theme(legend.position="none")

p8 <- ggplot(samp_grate, aes(beta, S0)) +
    geom_density_2d_filled(adjust=2.2) +
    geom_point(data=slice_sample(samp_grate, n=100), col="gray40", shape=1) +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\beta$", y="$S_0$", fill="Density") +
    theme(legend.position="none")

# Assemble plot-------------------------------------------------------------------

rt <- str_c(
    "\\textbf{\\large ", 
    c("Reproductive number", "Outbreak size", "Peak intensity", "Growth rate"), 
    "}"
)
row_titles <- rep("", 8)
row_titles[seq(1, 7, 2)] <- rt
plots <- list(p1, p2, p3, p4, p5, p6, p7, p8)
# Add row titles
plots <- map2(plots, row_titles, ~arrangeGrob(.x, top=.y))

tikz_plot(grid.arrange(grobs=plots, nrow=4), "cond-samples.pdf", 5.5, 9)

# TODO: for method paper, show "number of fit ODEs" for both RMD and variance of an estimator using MC+Maximum likelihood. Between fitting for each gradient (and possibly being more biased) the MC method may take more time
