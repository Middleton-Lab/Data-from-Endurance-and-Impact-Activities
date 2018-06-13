library(readxl)
library(tidyverse)
library(stringr)
library(forcats)
library(cowplot)
library(rethinking)
library(reshape2)
library(ggsci)
library(DescTools)
source("https://raw.githubusercontent.com/kmiddleton/kmmisc/master/R/HPDI_test.R")

set.seed(237486)

iter <- 2 * 2.1e5
warmup <- 2 * 1e4
chains <- 1
cores <- 1

# Function to clean up map2stan tmp saved files
cleanup_map2stan <- function(fm) {
  model_code <- fm@stanfit@stanmodel@model_code
  model_name2 <- attr(model_code, "model_name2")
  unlink(paste0(tempdir(), "/", model_name2, c(".rds", ".stan")))
  unlink(paste0(tempdir(), "/", model_name2, c(".cpp", ".o", ".so")))
}

M <- read_excel("femur sectors_obs2.xlsx")%>%
  filter(Died != 1) %>%
  filter(!is.na(sector_1_P)) %>%
  dplyr::select(-Died, -Green, -Red) %>% 
  as.data.frame()

N <- read_excel("femur sectors_obs1.xlsx")%>%
  filter(Died != 1) %>%
  filter(!is.na(sector_1_P)) %>%
  dplyr::select(-Died, -Green, -Red) %>% 
  as.data.frame()

Obs1 <- M %>% dplyr::select(-Mouse_ID, -Group) %>%
  as.matrix() %>%
  as.numeric()
Obs2 <- N %>% dplyr::select(-Mouse_ID, -Group) %>%
  as.matrix() %>%
  as.numeric()

CohenKappa(Obs1, Obs2, conf.level = 0.95)
##CI= 0.84; .819-.872, "excellent"
 
cor.test(Obs1, Obs2, method= "spearman", exact = FALSE)
## p < 2.2e-16, rho= 0.913

sum(Obs1 != Obs2) / length(Obs1) * 100
## percent disagreement= 9.46%

# Analyze Obs 2
N <- M

## Convert to factors
N <- N %>%
  mutate(Mouse_ID = factor(Mouse_ID),
         Group = fct_relevel(Group, "Sedentary")) %>%
  mutate_all(factor)

## Melt dataframe
N_long <- melt(N, id.vars = c("Mouse_ID", "Group"),
               variable.name = "Sector",
               value.name = "Score")

## Split sector column
N_long$SectorNum <- gsub("sector_", "", N_long$Sector)
N_long$SectorNum <- gsub("_E", "", N_long$SectorNum)
N_long$SectorNum <- gsub("_P", "", N_long$SectorNum)

N_long$SectorSurf <- substr(as.character(N_long$Sector), 10, 10)

## Make graphs
N_long %>%
  ggplot(aes(x = Group, y = Score)) +
  geom_point(position = position_jitter(width = 0.1, height = 0.1)) +
  facet_grid(SectorSurf ~ SectorNum) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Dots sized by number in each group
N_long %>% group_by(SectorNum, SectorSurf, Group, Score) %>%
  tally() %>%
  ggplot(aes(x = Group, y = Score, size = n)) +
  geom_point() +
  facet_grid(SectorSurf ~ SectorNum) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Table of counts
N_long_E <- N_long %>% filter(SectorSurf == "E")
N_long_P <- N_long %>% filter(SectorSurf == "P")

ftable(xtabs(~ Score + Group + SectorNum, data = N_long_P))
ftable(xtabs(~ Score + Group + SectorNum, data = N_long_E))

## Bayesian model ##################################################

# Variables to loop through
vars <- names(N)[!(names(N) %in% c("Mouse_ID", "Group"))]

# Empty list to hold models
fms <- list()

# Open a plot for diagnostics
pdf("ord_logistic_diagnostics.pdf")

for (ii in 1:length(vars)) {  ### SLOW LOOP ###

  # Variable for this iteration
  v <- vars[ii]

  # Setup model data
  Nb <- N %>%
    select_(v, "Group", "Mouse_ID") %>%
    rename_("sector" = v) %>%
    mutate(sector = coerce_index(sector),
           GroupRun = if_else(Group == "Run", 1L, 0L),
           GroupImpact = if_else(Group == "Impact", 1L, 0L),
           GroupSwim = if_else(Group == "Swim", 1L, 0L),
           Mouse_ID = coerce_index(factor(Mouse_ID))) %>%
    dplyr::select(-Group) %>%
    drop_na(sector) %>% 
    as.data.frame()

  print(simplehist(Nb$sector, xlim = c(1, 3) , xlab = v))

  fm <- rethinking::map2stan(
    alist(
      sector ~ dordlogit(phi, cutpoints),
      phi <-
        aRun * GroupRun +
        aImpact * GroupImpact +
        aSwim * GroupSwim,
      cutpoints ~ normal(0, 10),
      c(aRun, aImpact, aSwim) ~ normal(0, 10)
    ),
    data = Nb,
    start = list(cutpoints = c(0, 1)),
    WAIC = FALSE,
    iter = iter,
    warmup = warmup,
    chains = chains,
    cores = cores,
    control = list(adapt_delta = 0.95),
    rng_seed = 324876
  )

  print(traceplot(fm@stanfit, pars = c("aRun", "aImpact", "aSwim")))

  fms[[ii]] <- list(sector = v,
                    fm = fm)

  cleanup_map2stan(fm)
}
dev.off()

save(fms, file = "Sector_Scoring_Models_No_Intercept.Rda")

####################################################################

load("Sector_Scoring_Models_No_Intercept.Rda")

prob <- 0.89

# Data frame to hold summaries
D <- data_frame(
  Sector = character(length = length(fms)),
  Run_L = numeric(length = length(fms)),
  Run_M = numeric(length = length(fms)),
  Run_U = numeric(length = length(fms)),
  Run_Q = numeric(length = length(fms)),
  Impact_L = numeric(length = length(fms)),
  Impact_M = numeric(length = length(fms)),
  Impact_U = numeric(length = length(fms)),
  Impact_Q = numeric(length = length(fms)),
  Swim_L = numeric(length = length(fms)),
  Swim_M = numeric(length = length(fms)),
  Swim_U = numeric(length = length(fms)),
  Swim_Q = numeric(length = length(fms))
)

pdf("ord_logistic_posteriors.pdf")

for (ii in 1:length(fms)) {
  v <- fms[[ii]]$sector
  print(v)
  print(precis(fms[[ii]]$fm))
  post <- extract.samples(fms[[ii]]$fm) %>% as.data.frame()

  p <- post %>%
    dplyr::select(aRun, aImpact, aSwim) %>%
    gather(Group, value) %>%
    ggplot(aes(value, color = Group)) +
    geom_density() +
    labs(title = v)
  print(p)

  D[ii, 1] <- v
  D[ii, 2:ncol(D)] <-
    c(HPDI(post$aRun, prob = prob)[1],
      median(post$aRun),
      HPDI(post$aRun, prob = prob)[2],
      HPDI_test(post$aRun),
      HPDI(post$aImpact, prob = prob)[1],
      median(post$aImpact),
      HPDI(post$aImpact, prob = prob)[2],
      HPDI_test(post$aImpact),
      HPDI(post$aSwim, prob = prob)[1],
      median(post$aSwim),
      HPDI(post$aSwim, prob = prob)[2],
      HPDI_test(post$aSwim))
}
dev.off()

write_csv(D, path = "Sector_Scoring_Output.csv")

####################################################################
## Plotting

D <- read_csv("Sector_Scoring_Output.csv")

dodge_width <- 0.5
font_size <- 10

my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size, face = "bold"),
  legend.text = element_text(size = font_size - 1),
  legend.title = element_text(size = font_size, face = "bold"),
  plot.title = element_text(size = font_size + 1),
  strip.text = element_text(size = font_size - 1)
)

my_scale <-
  scale_y_continuous(limits = c(-4.5, 4.5),
                     breaks = pretty(seq(-4, 4, by = 1)))


Dl <- D %>%
  mutate(Sector = str_replace(Sector, "sector_", ""),
         Sector = str_replace(Sector, "_", ""),
         Sector = fct_inorder(Sector)) %>%
  gather(variable, value, -Sector) %>%
  mutate(Group = str_split(variable, "_", simplify = TRUE)[, 1],
         wh = str_split(variable, "_", simplify = TRUE)[, 2]) %>%
  dplyr::select(-variable)

Dl %>%
  spread(wh, value) %>%
  ggplot(aes(x = Sector, shape = Group, color = Group)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_errorbar(aes(ymin = L, ymax = U), position = "dodge",
                width = dodge_width, size = 1) +
  geom_point(aes(y = M),
             position = position_dodge(width = dodge_width),
             size = 3) +
  scale_shape_discrete(solid = TRUE) +
  theme(legend.position = c(0.07, 0.08),
        plot.title = element_text(size = font_size + 1),
        legend.title = element_blank(),
        legend.key = element_rect(size = 4),
        legend.key.size = unit(1.2, 'lines')) +
  scale_color_nejm() +
  labs(y = "log Odds",
       title = "Log Odds by Sector") +
  my_theme +
  my_scale

# JEB MS Fig
D %>% 
  gather(key, value, -Sector) %>% 
  mutate(Sector = str_replace(Sector, "sector", ""),
         Sector = str_replace_all(Sector, "_", "")) %>% 
  mutate(Sector_num = str_sub(Sector, 1, 1),
         PeriEnd = str_sub(Sector, 2, 2),
         PeriEnd = if_else(PeriEnd == "E", "Endosteal", "Periosteal"),
         Group = str_split(key, "_", simplify = TRUE)[, 1],
         bound = str_split(key, "_", simplify = TRUE)[, 2]) %>% 
  dplyr::select(-Sector, -key) %>% 
  dcast(Sector_num + PeriEnd + Group ~ bound) %>% 
  ggplot(aes(Sector_num, color = Group, group = Group)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = L, ymax = U),
                position = position_dodge(width = 0.5),
                width = 0) +
  geom_point(aes(y = M), position = position_dodge(width = 0.5)) +
  facet_grid(PeriEnd ~ .) +
  scale_color_nejm() +
  labs(y = "log Odds",
       x = "Sector") +
  my_theme
ggsave(last_plot(), file = "../../Figures and tables/Figure_5.pdf",
       height = 3, width = 181/25.4)
