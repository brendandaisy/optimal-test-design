library(tidyverse)
library(ggthemes)

true_inf <- c(
    0.01, 0.017234986028274994, 0.029226916688214286, 0.0482619354074302, 0.07644676001950913, 0.11397801673090434, 0.15694203765403592,
    0.197068585361762, 0.22556384374277325, 0.23791293991026669, 0.2351898235396991, 0.22162345300496464, 0.2018696959400936, 
    0.1795447345733238, 0.15704515578805922, 0.13574849445435141, 0.11635024064541434, 0.09910377645846799, 0.08402199574081326, 0.0709903078468306, 0.059817465836561626, 0.05029848904490269, 0.04222657248401512, 0.03540215277236312, 
    0.029650387465228304, 0.024814410446680653, 0.02075294687536625, 0.017345265065778624, 0.014491436294653753, 0.012103619476697134, 0.010105741535796938
)
peak <- which.max(true_inf)
tsteps <- seq(2, length(true_inf)-1, by=2)

(inct_org <- read_csv("_research/tmp/res.csv", na="missing"))
(inct_org <- read_csv("data/sims/increasing-tspan/results-03-18.csv"))

inct <- inct_org |>
    select(sig, obs_model, obs_params, known, true=θtrue) |>
    mutate(sig=str_replace_all(sig, "Any|\\[|\\]", "")) |>
    separate(sig, str_c("u", tsteps), sep=",", convert=TRUE, fill="right") |>
    pivot_longer(matches("u\\d+"), names_to="t", values_to="SIG") |>
    mutate(t=as.double(str_extract(t, "\\d+")))

mylab <- function(x) paste0("Unknown: ", x)

inct |>
    mutate(
        rate=ifelse(str_detect(obs_params, "r"), str_extract(obs_params, "r = \\d+"), "r = Inf (Poisson)"),
        ntest=str_extract(obs_params, "n = \\d+")
    ) |>
    # filter(obs_model == "neg_binom") |>
    ggplot(aes(t, SIG, col=rate)) +
    geom_vline(xintercept=peak, col="orange", linetype="dashed") +
    geom_line() +
    geom_point() +
    facet_grid(ntest~known, scales="free_y", labeller=labeller(param_comb=as_labeller(mylab))) +
    labs(x="Days of observation", y="Shannon Information Gain", col="Dispersion")

## Marginal plots ##

known_lab <- list(
    "Set([:α])" = "Unknown: beta, S0",
    "Set([:β])" = "Unknown: alpha, S0",
    "Set([:S₀])" = "Unknown: alpha, beta",
    "Set{Symbol}()" = "Unknown: alpha, beta, S0"
)

recover_sig <- function(df, var) {
    df |> mutate("{{var}}":=str_replace_all({{var}}, "Any|\\[|\\]", "")) |>
        separate({{var}}, str_c("u", as_label(enquo(var)), tsteps), sep=",", convert=TRUE, fill="right")
        # pivot_longer(matches("u\\d+"), values_to=as_label(enquo(var))) |>
        # mutate(t=as.double(str_extract(name, "\\d+"))) |>
        # select(-name)
}

(marg_org <- read_csv("_research/tmp/res.csv", na="missing"))
(marg_org <- read_csv("data/sims/increasing-tspan-marg/results-03-25.csv", na="missing"))

marg <- marg_org |>
    mutate(known=map_chr(known, ~known_lab[[.x]])) |>
    select(-path, true=θtrue, sig_beta=`sig-β`, sig_alpha=`sig-α`, sig_S0=`sig-S₀`)

wmarg <- marg |>
    recover_sig(sig_beta) |>
    recover_sig(sig_alpha) |>
    recover_sig(sig_S0) |>
    pivot_longer(contains("usig"), values_to="SIG") |>
    mutate(var=str_extract(name, "sig_[a-zS]+"), t=as.double(str_extract(name, "\\d+")))

wmarg |>
    # filter(obs_model == "poisson") |>
    drop_na() |>
    mutate(
        rate=ifelse(str_detect(obs_params, "r"), str_extract(obs_params, "r = \\d+"), "r = Inf (Poisson)"),
        ntest=str_extract(obs_params, "n = \\d+")
    ) |>
    # filter(ntest != "n = 10000") |>
    ggplot(aes(t, SIG, col=var, group=interaction(var, rate, gitcommit))) +
    geom_vline(xintercept=peak, col="orange", linetype="dashed") +
    # geom_line(aes(linetype=rate)) +
    geom_point(alpha=0.5, shape=1) +
    facet_grid(ntest~known, scales="free_y") +
    labs(x="Days of observation", y="Shannon Information Gain", col="Dispersion")

