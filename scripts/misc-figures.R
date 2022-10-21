library(tidyverse)
library(ggthemes)
library(tikzDevice)

samp <- read_csv("data/samples-cond-R=3.csv")

samp |> 
    mutate(across(everything(), ~ifelse(.x == 0, NA, .x)))

gg <- ggplot(samp, aes(alpha, beta)) +
    geom_density_2d_filled() +
    scale_fill_viridis_d() +
    scale_x_continuous(expand=expansion(0.0001)) +
    scale_y_continuous(expand=expansion(0.0001)) +
    labs(x="$\\alpha$", y="$\\beta$", fill="Density") +
    theme_bw()

tikz_plot(gg, "plots/samples-cond-R=3", 4, 3.2)
