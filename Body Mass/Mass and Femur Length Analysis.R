library(readxl)
library(ggplot2)
library(dplyr)
library(reshape2)
library(nlme)
library(car)
library(ggsci)
library(rstanarm)
library(rethinking)
library(tidyr)
library(tidyverse)
library(latex2exp)
library(multcomp)
library(lsmeans)
library(cowplot)

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

font_size <- 10

my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size, face = "bold"),
  legend.text = element_text(size = font_size - 1),
  legend.title = element_text(size = font_size, face = "bold"),
  plot.title = element_text(size = font_size + 1),
  strip.text = element_text(size = font_size - 1)
)

M <- read_excel("MDII_ICR measurements.xlsx", na = "NA")

# Drop mice 7 & 34, which were euthanized prior to the end of the
# experiment
M <- M %>% 
  filter(MouseID != 7) %>% 
  filter(MouseID != 34)

# Convert to factors
M$Group <- factor(M$Group)
M$MouseID <- factor(M$MouseID)
head(M)

# Convert from wide to long
## First, drop all the columns that don't have to do with mass
N <- M %>% 
  dplyr::select(-day13mass, -day18mass, -femur_length,
                -ML_diam, -AP_diam, -distal_width,
                -proximal_width, -head_height, -head_depth)

glimpse(N)

## Convert to Long format
N_long <- melt(N, id.vars = c("MouseID", "Group"),
               variable.name = "Day",
               value.name = "Mass")

N_long$Day_Num <- (gsub("mass", "", N_long$Day))
N_long$Day_Num <- as.numeric(gsub("day", "", N_long$Day_Num))
N_long$Day_Num[is.na(N_long$Day_Num)] <- -7

N_long$Day <- as.character(N_long$Day)

str(N_long)

# Drop the one with broken condyles for using femur length
O <- M %>% filter(MouseID != 49)

# Relevel so Control is last.
N_long$Group <- fct_relevel(N_long$Group, "Control", after = Inf)

## Wean mass
ggplot(M, aes(x= wean_mass)) + 
  geom_histogram(binwidth = 1) +
  facet_grid(Group ~ .)

fm_wn <- map2stan(
  alist(
    wean_mass ~ dnorm(mu, sigma),
    mu <- a[Group],
    a[Group] ~ dnorm(25, 10),
    sigma ~ dcauchy(0, 2)
  ),
  data = N,
  iter = iter,
  warmup = warmup,
  control = list(max_treedepth = 15),
  WAIC = FALSE
)

wn <- extract.samples(fm_wn)[['a']] %>% as_tibble()
names(wn) <- c("Control", "Drop", "Run", "Swim")

smallest_HPDI(wn$Swim - wn$Control)
smallest_HPDI(wn$Run - wn$Control)
smallest_HPDI(wn$Drop - wn$Control)
median(wn$Swim - wn$Control)
median(wn$Run - wn$Control)
median(wn$Drop - wn$Control)

wn <- wn %>% gather()
wn$day <- -6

## Mass over time
massplot <- ggplot(N_long, aes(x = Day_Num, y = Mass,
                               group = MouseID,
                               color = Group)) +
  geom_line(alpha = 0.25) +
  geom_smooth(aes(group= Group), method = "loess", se = FALSE, lwd = 2) +
  scale_y_continuous(limits = c(15, 43)) + 
  scale_x_continuous(limits = c(-7, 21),
                     breaks = c(unique(N_long$Day_Num))) + 
  xlab("Treatment Day") +
  ylab("Body Mass (g)") +
  theme(legend.position = c(0.82, 0.22),
        legend.text = element_text(size = 10),
        legend.background = element_rect(linetype = "solid",
                                         size = 0.25,
                                         color = "black"))+
  scale_color_nejm(labels = c("Impact", "Run", "Swim", "Control"),
                       name = element_blank()) +
  my_theme

massplot

# Save MS SI Fig 1
save_plot(file = "../../Figures and tables/Figure_SI1.pdf",
          massplot, base_height = 3.5, base_width = 180/25.4)

## ------------------------------------------------------------------------
mass21plot <- ggplot(M, aes(x = Group, y = day21mass)) +
  geom_boxplot() +
  xlab("Treatment Group") +
  ylab("Body Mass on Day 21")
mass21plot

## ------------------------------------------------------------------------
M_b <- M %>% 
  mutate(Group_num = coerce_index(Group),
         mass_c = day21mass - mean(day21mass)) %>%
  dplyr::select(day21mass, mass_c, Group_num) %>% 
  as.data.frame()

fm_mass_stan <-
  map2stan(
    alist(
      mass_c ~ dnorm(mu, sigma),
      mu <- aGroup[Group_num],
      aGroup[Group_num] ~ dnorm(0, 2),
      sigma ~ dcauchy(0, 2)
    ),
    data = M_b,
    iter = 1.1e4,
    warmup = 1e3,
    WAIC = FALSE
  )

precis(fm_mass_stan, depth = 2, digits = 4)

post <- extract.samples(fm_mass_stan)
Control <- post$aGroup[, 1]
Impact <- post$aGroup[, 2]
Run <- post$aGroup[, 3]
Swim <- post$aGroup[, 4]

data_frame(Control, Impact, Run, Swim) %>% 
  gather(Group, `day21mass`) %>% 
  ggplot(aes(`day21mass`, color = Group)) +
  geom_line(stat = "density")

median(Control - Impact)
median(Control - Run)
median(Control - Swim)
smallest_HPDI(Swim - Control)
smallest_HPDI(Run - Control)
smallest_HPDI(Impact - Control)

## ------------------------------------------------------------------------

# Femur length boxplot
ggplot(O, aes(x = Group, y = femur_length)) +
  geom_boxplot() +
  xlab("Treatment Group") + 
  ylab("Femur Length (mm)")

## Data points and error bars
ggplot(O, aes(x = Group, y = femur_length)) +
  xlab("Treatment Group") +
  ylab("Femur Length (mm)") +
  geom_point(position = position_jitter(width = 0.1)) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 0.7)

## ------------------------------------------------------------------------

M_f <- O %>% 
  mutate(Group_num = coerce_index(Group),
         fem_c = femur_length - mean(femur_length)) %>%
  dplyr::select(femur_length, fem_c, Group_num) %>% 
  as.data.frame()

fm_fem_stan <-
  map2stan(
    alist(
      fem_c ~ dnorm(mu, sigma),
      mu <- aGroup[Group_num],
      aGroup[Group_num] ~ dnorm(0, 2),
      sigma ~ dcauchy(0, 2)
    ),
    data = M_f,
    iter = 1.1e4,
    warmup = 1e3,
    WAIC = FALSE
  )

precis(fm_fem_stan, depth = 2, digits = 4)

post <- extract.samples(fm_fem_stan)
Control <- post$aGroup[, 1]
Impact <- post$aGroup[, 2]
Run <- post$aGroup[, 3]
Swim <- post$aGroup[, 4]

data_frame(Control, Impact, Run, Swim) %>% 
  gather(Group, `femur_length`) %>% 
  ggplot(aes(`femur_length`, color = Group)) +
  geom_line(stat = "density")

median(Control - Impact)
median(Control - Run)
median(Control - Swim)
HPDI(Control - Impact)
HPDI(Control - Run)
HPDI(Control - Swim)
