library(cmdstanr)
library(posterior)
library(tidyverse)
library(VGAM)

rep_number <- function(alpha, beta) {
    beta / alpha
}

outbreak_size <- function(alpha, beta, S0) {
    R <- rep_number(alpha, beta)
    Rinit <- 1 - S0 - 0.01
    r_inf <- 1 + 1/R * lambertW(-S0 * R * exp(-R*(1-Rinit)))
    r_inf - Rinit
}

peak_intensity <- function(alpha, beta, S0) {
    R <- rep_number(alpha, beta)
    ifelse(
        R * S0 < 1,
        0.01,
        0.01 + S0 - 1/R*log(S0) - 1/R*(1-log(1/R))
    )
}

growth_rate <- function(alpha, beta, S0) beta*S0 - alpha

plot_transformation <- function(varname, priors, true) {
    ggplot(priors, aes(.data[[varname]])) +
        geom_density(size=1.3) +
        geom_vline(aes(xintercept=.data[[varname]]), data=true, linetype="dashed", color="orange", size=1.3) +
        labs(x=names(which(texlab == varname))[1], y=NULL) +
        theme_dens
}

priors <- tibble(
    α=runif(10000, 0.05, 0.85),
    β=runif(10000, 0.3, 1.5),
    `S₀`=runif(10000, 0.1, 0.99),
    `rep-number`=rep_number(α, β),
    `outbreak-size`=outbreak_size(α, β, `S₀`),
    `peak-intensity`=peak_intensity(α, β, `S₀`),
    `peak-timing`=read_csv("_research/tpeak-pri.csv")$V,
    `growth-rate`=growth_rate(α, β, `S₀`)
) |> 
    pivot_longer(everything(), "var") |>
    mutate(var=fct_relevel(fct_recode(var, !!!texlab), !!!names(texlab)))

# priors <- tibble(
#     `rep-number`=rep_number(alpha_samp, beta_samp),
#     `outbreak-size`=outbreak_size(alpha_samp, beta_samp, S0_samp),
#     `peak-intensity`=peak_intensity(alpha_samp, beta_samp, S0_samp),
#     `growth-rate`=growth_rate(alpha_samp, beta_samp, S0_samp)
# ) |> 
#     pivot_longer(everything(), "var") |>
#     bind_rows(
#         tibble(var="α", value=dunif(seq(0, 0.9, 0.01), 0.05, 0.85)),
#         tibble(var="β", value=dunif(seq(0.25, 1.55, 0.01), 0.3, 1.5)),
#         tibble(var="S₀", value=dunif(seq(0.05, 1.05, 0.01), 0.1, 0.99))
#     ) |> 
#     mutate(var=fct_relevel(fct_recode(var, !!!texlab), !!!names(texlab)))

alpha_true <- 0.2
beta_true <- 1.25
S0_true <- 0.6
true_vals <- tibble(
    α=alpha_true,
    β=beta_true,
    `S₀`=S0_true,
    `rep-number`=rep_number(α, β),
    `outbreak-size`=outbreak_size(α, β, `S₀`),
    `peak-intensity`=peak_intensity(α, β, `S₀`),
    `peak-timing`=9.29,
    `growth-rate`=growth_rate(α, β, `S₀`)
) |> 
    pivot_longer(everything(), "var") |>
    mutate(var=fct_relevel(fct_recode(var, !!!texlab), !!!names(texlab)))

###

stan_dat <- list(
    max_t=30,
    max_obs_t=8,
    ts = 0:8,
    y=c(7, 22, 38, 50, 59, 115, 163, 183, 242, 239, 244, 230, 213, 183, 136, 148, 116, 109, 78, 69, 52, 52, 34, 44, 26, 23, 16, 13, 17, 13, 12),
    I0= 0.01
)

exec <- cmdstan_model("scripts/fit-sir.stan")
fit <- exec$sample(data=stan_dat)

post8 <- as_draws_df(fit$draws()) |> 
    as_tibble() |> 
    select(`α`=alpha, `β`=beta, `S₀`=S0) |> 
    mutate(
        `rep-number`=rep_number(α, β),
        `outbreak-size`=outbreak_size(α, β, `S₀`),
        `peak-intensity`=peak_intensity(α, β, `S₀`),
        `peak-timing`=read_csv("_research/tpeak-post8.csv")$V,
        `growth-rate`=growth_rate(α, β, `S₀`),
    ) |> 
    pivot_longer(everything(), "var") |>
    mutate(var=fct_relevel(fct_recode(var, !!!texlab), !!!names(texlab)))

all_samps <- bind_rows(
    mutate(priors, dist="prior"),
    mutate(post3, dist="post3"),
    mutate(post8, dist="post8")
) |> 
    mutate(dist=fct_inorder(dist))

density_plot_var <- function(var, idx) {
    p <- all_samps |> 
        filter(var == !!var) |> 
        ggplot(aes(value, group=dist, col=dist)) +
        geom_density(size=1.1, adjust=1.3, alpha=0.7) +
        geom_vline(xintercept=filter(true_vals, var == !!var)$value, col="orange", linetype="dashed", size=1.45) +
        scale_color_manual(values=c("gray30", "#f53db5", "purple")) +
        labs(x=str_extract(var, "\\$.+\\$"), y=NULL) +
        theme_dens
    
    # if (idx %% 2 == 0)
    #     p <- p + theme(plot.background=element_rect(fill="lightblue", color="lightblue"))
    return(p)
}

theme_dens <- theme_bw() +
    theme(
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        legend.position="none",
        plot.margin=unit(c(0, 0, 0, 0), "in")
    )

plot_dens <- imap(true_vals$var, density_plot_var)
plot_dens[[1]] <- plot_dens[[1]] + scale_x_continuous(breaks=c(0.1, 0.45, 0.8))
plot_dens[[3]] <- plot_dens[[3]] + scale_x_continuous(breaks=c(0.1, 0.5, 0.9))
plot_dens[[5]] <- plot_dens[[5]] + scale_x_continuous(breaks=c(0, 0.5, 1))
plot_dens[[6]] <- plot_dens[[6]] + scale_x_continuous(breaks=c(0, 0.4, 0.8))
plot_dens[[8]] <- plot_dens[[8]] + scale_x_continuous(breaks=c(-0.5, 0.25, 1))

# plot_dens <- ggplot(all_samps, aes(value, col=var, group=dist)) +
#     geom_density(size=1.3, adjust=1.5, alpha=0.7) +
#     geom_vline(aes(xintercept=value), col="orange", data=true_vals, linetype="dashed", size=1.3) +
#     facet_wrap(~var, scales="free", nrow=1) +
#     # coord_cartesian(expand=FALSE) +
#     labs(x=NULL, y="Density") +
#     theme_dens
