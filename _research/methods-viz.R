library(tidyverse)
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

theme_svg <- theme(
    panel.grid = element_blank(),
    axis.line = element_line(arrow=arrow(type="closed", angle=17, length=unit(0.13, "in")), size=2),
    panel.background = element_rect(fill="white"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
)

pri_cdat <- pri_curves |>
    pivot_longer(contains("x")) |>
    mutate(lpdf=rep(pri_lpdf, 301))

gg <- ggplot(pri_cdat, aes(t, value, groups=name, col=lpdf)) +
    geom_line(alpha=0.7, size=1.17) +
    geom_line(
        aes(t, xtrue), tibble(t=seq(0, 11, 0.1), xtrue=xtrue[seq(1, 111, 1)]), 
        col="#ffa24b", linetype="dashed", inherit.aes=FALSE, size=2.5
    ) +
    geom_point(aes(t, y / 100), ydat, col="white", inherit.aes=FALSE, shape=21, fill="#d175ab", size=2.7) +
    # scale_color_gradientn(colors=c("#132B43", "#56B1F7", "#86f2fb"), values=c(0, 0.8, 1)) +
    scale_color_viridis_c() +
    theme_svg +
    theme(legend.position="none") +
    labs(y="{\\fontfamily{cmss}\\selectfont\\Large prevalence}", x="{\\fontfamily{cmss}\\selectfont\\Large time}")

tikz_plot(gg, "pri-curve", w=4.5, h=4)

############

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

#############

ydat <- tibble(
    y = y[1:12],
    t = 0:11    
)

cond_cdat <- mutate(bdf, t=seq(0, 30, 0.1)) |>
    pivot_longer(contains("x")) |>
    mutate(lpdf=rep(βlpdf, 301))

gg <- ggplot(cond_cdat, aes(t, value, col=lpdf, group=name)) +
    geom_line(alpha=0.7, size=1.17) +
    geom_line(
        aes(t, xtrue), tibble(t=seq(0, 11, 0.1), xtrue=xtrue[seq(1, 111, 1)]), 
        col="#ffa24b", linetype="dashed", inherit.aes=FALSE, size=2.2
    ) +
    geom_point(aes(t, y / 100), ydat, col="white", inherit.aes=FALSE, shape=21, fill="#d175ab", size=2.7) +
    # scale_color_gradientn(colors=c("#132B43", "#56B1F7", "#86f2fb"), values=c(0, 0.8, 1)) +
    scale_color_viridis_c() +
    theme_svg +
    theme(legend.position="none") +
    labs(y="{\\fontfamily{cmss}\\selectfont\\Large prevalence}", x="{\\fontfamily{cmss}\\selectfont\\Large time}")

tikz_plot(gg, "cond-curves", w=4.5, h=4)

##########

lpdfs <- tibble(
    pl = c(rep(5, 15), rnorm(100, 10, 25), rep(32, 15)),
    cl = c(rep(5, 12), rnorm(100, 25, 10))
)

gg <- ggplot(pivot_longer(lpdfs, everything()), aes(value, group=name)) +
    geom_density(col=alpha("black", 0.8), size=1.7) +
    geom_vline(xintercept=26.5, col="#d175ab", linetype="dashed", size=2.2) +
    geom_segment(
        aes(x=28, y=0.0135, xend=28, yend=0.034), data.frame(),
        col=alpha("black", 0.8), size=1.7, inherit.aes=FALSE, 
        arrow=arrow(type="closed", angle=15, length=unit(0.07, "in"), ends="both")
    ) +
    # facet_wrap(~name, scales="free") +
    theme_svg +
    labs(y="{\\fontfamily{cmss}\\selectfont\\Large density}", x="{\\fontfamily{cmss}\\selectfont\\Large data}")

tikz_plot(gg, "marg-lik", w=4.5, h=4)

########

p3dat <- tibble(
    # S0samp = S0_samp
    betasamp = beta_samp
    # alphasamp = alpha_samp
)

gg <- p3dat |>
    # pivot_longer(everything()) |>
    ggplot(aes(betasamp)) +
    geom_density(col=alpha("black", 0.8), size=1.7) +
    stat_function(fun=~dunif(.x, 0.3, 3)-0.1, alpha=0.8, size=1.3) +
    geom_vline(
        # aes(xintercept=p), tibble(p=c(0.6, 1.25, 0.2), name=c("S0samp", "betasamp", "alphasamp")),
        aes(xintercept=p), tibble(p=c(1.25), name=c("betasamp")),
        col="#ffa24b", linetype="dashed", size=2.2
    ) +
    geom_segment(
        aes(x=1.3, y=dunif(1.3, 0.3, 3)-0.09, xend=1.3, yend=0.41), data.frame(),
        col=alpha("black", 0.8), size=1.7, inherit.aes=FALSE, 
        arrow=arrow(type="closed", angle=15, length=unit(0.07, "in"), ends="both")
    ) +
    # facet_wrap(~name, scales="free") +
    theme_svg +
    labs(y="{\\fontfamily{cmss}\\selectfont\\Large density}", x=NULL)

tikz_plot(gg, "post-dist", w=4.5, h=4)

## phases of an outbreak

phase <- tibble(t=seq(0, 14, 0.1), inf=xtrue[seq(1, 141, 1)])

gg <- ggplot(phase, aes(t, inf)) +
    geom_vline(xintercept=3, col="black", linetype="dotted", size=2.5) +
    geom_vline(xintercept=8, col="black", linetype="dotted", size=2.5) +
    stat_function(fun=~0.01 * exp((1.25*0.6 - 0.2) * .x), col="gray20", alpha=0.8) +
    geom_line(col="#ffa24b", size=2.2) +
    theme_svg +
    lims(y=c(0, 0.5)) +
    labs(y="{\\fontfamily{cmss}\\selectfont\\Large prevalence}", x="{\\fontfamily{cmss}\\selectfont\\Large time}")

tikz_plot(gg, "outbreak-phases", w=5.7, h=2.55)