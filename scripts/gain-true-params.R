library(tidyverse)
library(ggthemes)
library(tikzDevice)
library(ggprism)
library(cowplot)

res1 <- read_csv("data/sims/gain-var-beta.csv") |> 
    rename(alpha=α, beta=β, S0=`S₀`) |> 
    mutate(Reff=beta*S0/alpha)

res <- read_csv("data/sims/gain-var-reff.csv") |> 
    rename(alpha=α, beta=β, S0=`S₀`) |> 
    mutate(Reff=beta*S0/alpha)

res <- bind_rows(
    mutate(res1, run="1"),
    mutate(res, run="2")
)

res |> 
    filter(Reff > 0.95) |> 
    ggplot(aes(Reff, md, col=beta)) +
    geom_point() +
    facet_wrap(~lab, nrow=2)
    
ggsave("plots/gain-var-beta1.pdf", width=6.6, height=4.5)
