library(tidyverse)
library(plotly)

colnames(likdf) <- c("S", "b", "lik3", "lik8")
btrue <- 1.25
Strue <- 0.6
atrue <- 0.2

gg <- ggplot(likdf, aes(S, b, col=lik3)) +
    geom_point(size=2.5, alpha=0.7) +
    stat_function(fun=~btrue*Strue / .x) +
    geom_point(aes(Strue, btrue), tibble(S=Strue, b=btrue), col="#ffa24b", size=2.7, shape=8) +
    scale_color_viridis_c() +
    scale_x_continuous(expand = expansion(), limits=c(0.09, 0.91)) +
    scale_y_continuous(expand = expansion(), limits=c(0.29, 3)) +
    theme_bw() +
    labs(col="$P(y\\mid S_0, \\beta, \\alpha^*)$", x="$S_0$", y="$\\beta$")

tikz_plot(gg, "early-likelihood", w=4.3, h=3.2)

gg <- ggplot(likdf, aes(S, b, col=lik8)) +
    geom_point(size=2.5, alpha=0.7) +
    stat_function(fun=~btrue*Strue / .x) +
    geom_point(aes(Strue, btrue), tibble(S=Strue, b=btrue), col="#ffa24b", size=2.7, shape=8) +
    scale_color_viridis_c() +
    scale_x_continuous(expand = expansion(), limits=c(0.09, 0.91)) +
    scale_y_continuous(expand = expansion(), limits=c(0.29, 3)) +
    theme_bw() +
    labs(col="$P(y\\mid S_0, \\beta, \\alpha^*)$", x="$S_0$", y="$\\beta$")

tikz_plot(gg, "mid-likelihood", w=4.3, h=3.2)

# beta <- seq(0.1, 3, 0.1)
# S0 <- seq(0, 1, 0.1)

# btrue <- 1.25
# Strue <- 0.6
# atrue <- 0.2

# dist(rbind(c(btrue, Strue, atrue), c(0.6, 1, 0.0500)))

# surf <- crossing(beta, S0) |>
#     mutate(alpha = beta*S0 - (btrue*Strue - atrue)) |>
#     filter(alpha > 0) |>
#     mutate(lik = beta - S0 + alpha)

# plot_ly(x=~surf$beta, y=~surf$S0, z=~surf$alpha, type="surface", mode="markers", color=surf$lik)
