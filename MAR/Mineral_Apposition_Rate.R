library(tidyverse)
library(cowplot)
library(readxl)
library(rethinking)
library(forcats)
library(stringr)
library(reshape2)
library(ggsci)

rerun_stan <- TRUE

options(warn = 1)

font_size <- 10

smallest_HPDI <- function(x) {
  min_interval <- 0.01
  for (ii in seq(0.01, 0.99, 0.01)) {
    Int <- HPDI(x, prob = ii)
    if (abs(sum(sign(Int))) == 2) {
      min_interval <- ii
    }
  }
  return(min_interval)
}

M <- read_excel("MAR measurements.xlsx",
                na = "NA")

# Tally of # of measurements
M %>% group_by(MouseID, Treatment) %>% tally() %>% as.data.frame()

Txs <- M %>% select(MouseID, Position, Treatment) %>% 
  filter(Position == 1) %>% 
  dplyr::select(-Position)

# Aggregate mean MAR per sector
mean_std <- function(x) {
  return(mean(x / 5))
}

M_agr <- M %>%
  dplyr::select(-Position, -Treatment) %>% 
  group_by(MouseID) %>% 
  summarise_all(mean_std) %>% 
  ungroup()

# Plot histogram of MARs
M_agr %>%
  dplyr::select(-MouseID) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_grid(. ~ key)

# Join treatments back in, convert to factor & numeric representation
N <- left_join(M_agr, Txs)

N <- N %>% 
  dplyr::select(-MouseID) %>% 
  gather(-Treatment, key = "Sector", value = "MAR_mean") %>% 
  drop_na()

N$Treatment <- as.factor(N$Treatment)
N <- N %>% 
  mutate(Treatment = fct_relevel(Treatment, "Sedentary", after = 0),
         Treatment_num = as.integer(Treatment),
         MAR_mean_c = MAR_mean - mean(MAR_mean)) %>% 
  as.data.frame()

head(N)

###

if (rerun_stan) {
  
  Sectors <- unique(N$Sector)
  
  n <- length(Sectors)
  
  out <- tibble(
    Sector = Sectors,
    dImp_l = numeric(n),
    dImp_m = numeric(n),
    dImp_u = numeric(n),
    dImp_interval = numeric(n),
    dRun_l = numeric(n),
    dRun_m = numeric(n),
    dRun_u = numeric(n),
    dRun_interval = numeric(n),
    dSw_l = numeric(n),
    dSw_m = numeric(n),
    dSw_u = numeric(n),
    dSw_interval = numeric(n))
  
  pdf("Mineral_Apposition_Rate.pdf")
  
  for (ii in 1:n) {
    S <- Sectors[ii]
    message(paste("\n", S, "\n"))
    N_s <- N %>% filter(Sector == S)
    
    fm <- map2stan(
      alist(
        MAR_mean_c ~ dnorm(mu, sigma),
        mu <- aTx[Treatment_num],
        aTx[Treatment_num] ~ dnorm(0, 5),
        sigma ~ dcauchy(0, 2)
      ),
      data = N_s,
      chains = 1,
      iter = 420000,
      warmup = 20000,
      WAIC = FALSE,
      control = list(max_treedepth = 15)
    )
    print(plot(fm))
    post <- extract.samples(fm)
    dImp <- post$aTx[, 2] - post$aTx[, 1]
    dRun <- post$aTx[, 3] - post$aTx[, 1]
    dSw <- post$aTx[, 4] - post$aTx[, 1]
    
    out[ii, "dImp_l"] <- HPDI(dImp)[1]
    out[ii, "dImp_m"] <- median(dImp)
    out[ii, "dImp_u"] <- HPDI(dImp)[2]
    out[ii, "dImp_interval"] <- smallest_HPDI(dImp)
    out[ii, "dRun_l"] <- HPDI(dRun)[1]
    out[ii, "dRun_m"] <- median(dRun)
    out[ii, "dRun_u"] <- HPDI(dRun)[2]
    out[ii, "dRun_interval"] <- smallest_HPDI(dRun)
    out[ii, "dSw_l"] <- HPDI(dSw)[1]
    out[ii, "dSw_m"] <- median(dSw)
    out[ii, "dSw_u"] <- HPDI(dSw)[2]
    out[ii, "dSw_interval"] <- smallest_HPDI(dSw)
  }
  dev.off()
  
  write_csv(path = "MAR_Analysis.csv", out)
}

#### Plot

out <- read_csv("MAR_Analysis.csv",
                col_types = "cdddddddddddd")

my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size, face = "bold"),
  legend.text = element_text(size = font_size - 1),
  legend.title = element_text(size = font_size, face = "bold"),
  plot.title = element_text(size = font_size + 1),
  strip.text = element_text(size = font_size - 1)
)


out %>% 
  gather(key, value, -Sector) %>% 
  mutate(Sector_num = str_sub(Sector, 1, 1),
         PeriEnd = str_sub(Sector, 2, 2),
         PeriEnd = if_else(PeriEnd == "E", "Endosteal", "Periosteal"),
         Group = str_split(key, "_", simplify = TRUE)[, 1],
         bound = str_split(key, "_", simplify = TRUE)[, 2],
         Group = case_when(
           Group == "dImp" ~ "Impact",
           Group == "dRun" ~ "Run",
           Group == "dSw" ~ "Swim"
         )) %>% 
  dplyr::select(-Sector, -key) %>% 
  dcast(Sector_num + PeriEnd + Group ~ bound) %>% 
  ggplot(aes(Sector_num, color = Group, group = Group)) +
  geom_hline(yintercept = 0, size= .25) +
  geom_errorbar(aes(ymin = l, ymax = u),
                position = position_dodge(width = 0.5),
                width = 0) +
  geom_point(aes(y = m), position = position_dodge(width = 0.5)) +
  facet_grid(PeriEnd ~ .) +
  scale_color_nejm() +
  labs(y = "Difference from Control (µm/day)",
       x = "Sector") +
  my_theme

ggsave(last_plot(), file = "../../Figures and tables/Figure_4.pdf",
       height = 3, width = 181/25.4)

## Absolute MAR means per sector and treatment
MAR_Mean_mean <- N %>% group_by(Treatment, Sector) %>%
  summarise(MAR_Mean_mean = mean(MAR_mean))
write_csv(MAR_Mean_mean, path = "MAR_Mean_mean.csv")
