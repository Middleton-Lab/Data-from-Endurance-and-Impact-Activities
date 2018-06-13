library(readxl)
library(tidyverse)
library(dplyr)
library(broom)
library(cowplot)
library(latex2exp)
library(rethinking)
source("https://raw.githubusercontent.com/kmiddleton/kmmisc/master/R/HPDI_test.R")
library(knitr)

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

rerun_stan <- FALSE
plot_diagnostics <- FALSE

iter <- 2.5e4
warmup <- 5e3
chains <- 4
cores <- 4

## ----load_XS_data--------------------------------------------------------
XSData <- read_excel("Plasticity Mice BoneJ.xlsx") %>% 
  dplyr::select(-Tx, -Scan_ID, -Contrast, -Order, -Label, -Bone_Code, -Slice)

# Take mean for 10 measurements for each mouse.
XSData <- XSData %>% 
  group_by(MouseID) %>% 
  summarise_all(mean) %>% 
  as.data.frame()

## ----load_dissection_data------------------------------------------------
WeanData <- read_excel("MDII_ICR measurements.xlsx", na = "NA")
M <- left_join(WeanData, XSData)

M <- M %>% mutate(Group = factor(Group)) %>% 
  as.data.frame()

M$J <- M$Imin + M$Imax
M$Elip <- (M$Imin/ M$Imax)

glimpse(M)

# Drop mice 7 & 34, which were euthanized prior to the end of the
# experiment
N <- M %>% 
  filter(MouseID != 7) %>% 
  filter(MouseID != 34)

# Variables to analyze in addition to mass and femur length
vars <- c("ML_diam", "AP_diam", "distal_width",
          "proximal_width", "head_height", "head_depth",
          "CSA", "Imin", "Imax", "Zpol",
          "Zmin", "Zmax", "J", "Elip")

# +2 for mass and femur length
nrows <- 2 + length(vars)

# Set up data frame to hold output.
# For each, keep lower 89%, median, upper 89%, and quantile for test vs. 0
out <- tibble(Trait = character(nrows),
              centered_value = numeric(nrows),
              slope_lower_89 = numeric(nrows),
              slope_median = numeric(nrows),
              slope_upper_89 = numeric(nrows),
              slope_diff_0 = numeric(nrows),
              d_lower_89 = numeric(nrows),
              d_median = numeric(nrows),
              d_upper_89 = numeric(nrows),
              d_diff_0 = numeric(nrows),
              d_min_interval = numeric(nrows),
              r_lower_89 = numeric(nrows),
              r_median = numeric(nrows),
              r_upper_89 = numeric(nrows),
              r_diff_0 = numeric(nrows),
              r_min_interval = numeric(nrows),
              s_lower_89 = numeric(nrows),
              s_median = numeric(nrows),
              s_upper_89 = numeric(nrows),
              s_diff_0 = numeric(nrows),
              s_min_interval = numeric(nrows))

# Function to extract samples, plot distributions, and append output data frame.
save_post_est <- function(trait, centered_value, fm, out, r = 1) {
  post <- extract.samples(fm)
  
  Control <- post$aGroup[, 1]
  Drop <- post$aGroup[, 2]
  Run <- post$aGroup[, 3]
  Swim <- post$aGroup[, 4]
  
  p1 <- data_frame(Control, Drop, Run, Swim) %>% 
    gather(Group) %>% 
    ggplot(aes(value, color = Group)) +
    geom_line(stat = "density") +
    labs(title = trait, x = trait)
  print(p1)
  
  dDrop <- Drop - Control
  dRun <- Run - Control
  dSwim <- Swim - Control
  
  p2 <- data_frame(dDrop, dRun, dSwim) %>% 
    gather(Group) %>% 
    ggplot(aes(value, color = Group)) +
    geom_line(stat = "density") +
    labs(title = trait, x = "Difference from Control")
  print(p2)
  
  out$Trait[r] <- trait
  out$centered_value[r] <- centered_value
  out$d_lower_89[r] <- HPDI(dDrop)[1]
  out$d_median[r] <- median(dDrop)
  out$d_upper_89[r] <- HPDI(dDrop)[2]
  out$d_diff_0[r] <- HPDI_test(dDrop)
  out$d_min_interval[r] <- smallest_HPDI(dDrop)
  out$r_lower_89[r] <- HPDI(dRun)[1]
  out$r_median[r] <- median(dRun)
  out$r_upper_89[r] <- HPDI(dRun)[2]
  out$r_diff_0[r] <- HPDI_test(dRun)
  out$r_min_interval[r] <- smallest_HPDI(dRun)
  out$s_lower_89[r] <- HPDI(dSwim)[1]
  out$s_median[r] <- median(dSwim)
  out$s_upper_89[r] <- HPDI(dSwim)[2]
  out$s_diff_0[r] <- HPDI_test(dSwim)
  out$s_min_interval[r] <- smallest_HPDI(dSwim)
  
  # Extract slope and plot if present
  if ("bFemur" %in% names(post)) {
    bFemur <- post$bFemur
    out$slope_lower_89[r] <- HPDI(bFemur)[1]
    out$slope_median[r] <- median(bFemur)
    out$slope_upper_89[r] <- HPDI(bFemur)[2]
    out$slope_diff_0[r] <- HPDI_test(bFemur)
    p <- ggplot(tibble(bFemur), aes(bFemur)) +
      geom_line(stat = "density") +
      labs(title = trait)
    print(p)
  }
  
  return(out)
}

## ----mass----------------------------------------------------------------
ggplot(M, aes(x= Group, y= day21mass)) + 
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1.2) +
  labs(y = TeX("Body Mass (g)"))

M_b <- M %>% 
  drop_na(day21mass) %>% 
  mutate(Group_num = coerce_index(Group),
         mass_c = day21mass - mean(day21mass)) %>%
  dplyr::select(day21mass, mass_c, Group_num) %>% 
  as.data.frame()

if (rerun_stan) {
  fm_mass <-
    map2stan(
      alist(
        mass_c ~ normal(mu, sigma),
        mu <- aGroup[Group_num],
        aGroup[Group_num] ~ normal(0, 2),
        sigma ~ cauchy(0, 2)
      ),
      data = M_b,
      iter = iter,
      warmup = warmup,
      chains = chains,
      cores = cores,
      WAIC = FALSE,
      control = list(adapt_delta = 0.95)
    )
  
  if (plot_diagnostics) plot(fm_mass)
  precis(fm_mass, depth = 2, digits = 4)
  
  out <- save_post_est("Mass", mean(M$day21mass, na.rm = TRUE),
                       fm_mass, out, r = 1)
  save(fm_mass, file = "fm_mass.Rda")
} else {
  load("fm_mass.Rda")
}

## Femur length
ggplot(M, aes(x = Group, y = femur_length)) + 
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1.2) +
  labs(y = TeX("Femur Length (mm)"))

M_b <- M %>% 
  mutate(Group_num = coerce_index(Group),
         femur_length_c = femur_length - mean(femur_length,
                                              na.rm = TRUE)) %>%
  dplyr::select(femur_length, femur_length_c, Group_num) %>% 
  drop_na() %>% 
  as.data.frame()

if (rerun_stan) {
  fm_femur_length <-
    map2stan(
      alist(
        femur_length_c ~ normal(mu, sigma),
        mu <- aGroup[Group_num],
        aGroup[Group_num] ~ normal(0, 2),
        sigma ~ cauchy(0, 2)
      ),
      data = M_b,
      iter = iter,
      warmup = warmup,
      chains = chains,
      cores = cores,
      WAIC = FALSE,
      control = list(adapt_delta = 0.95)
    )
  
  if (plot_diagnostics) plot(fm_femur_length)
  precis(fm_femur_length, depth = 2, digits = 4)
  
  out <- save_post_est("Femur Length",
                       mean(M$femur_length,
                            na.rm = TRUE),
                       fm_femur_length, out, r = 2)
  save(fm_femur_length, file = "fm_femur_length.Rda")
} else {
  load("fm_femur_length.Rda")
}

## ----slow_loop-----------------------------------------------------------
for (ii in 1:length(vars)) {
  y <- vars[ii]
  
  # Simple ANCOVA model for diagnostic plot
  fm <- lm(M[, y] ~ femur_length + Group, M)
  M_fm <- augment(fm, M)
  
  pp <- ggplot(M_fm, aes_string(x = "femur_length", y = y,
                                color = "Group",
                                label = "MouseID")) +
    geom_line(aes_string(x = "femur_length", y = ".fitted")) +
    geom_point(size = 3) +
    labs(x = "Femur Length", y = y) +
    theme(legend.position = c(0.15, 0.8),
          legend.text = element_text(size = 12),
          legend.background = element_rect(linetype = "solid",
                                           size = 0.25,
                                           color = "black")) +
    scale_color_discrete(labels = c("Sedentary", "Impact", "Run", "Swim"),
                         name = element_blank())
  print(pp)
  
  M_tmp <- M
  names(M_tmp)[names(M_tmp) == y] <- "y"
  
  dat_s <- M_tmp %>%
    drop_na(femur_length, y) %>% 
    mutate(Group_num = coerce_index(Group),
           y_c = y - mean(y),
           femur_length_c = femur_length - mean(femur_length)) %>%
    dplyr::select(femur_length_c,
                  y_c,
                  Group_num) %>% 
    as.data.frame()
  if (rerun_stan) {
    fm <- map2stan(
      alist(
        y_c ~ normal(mu, sigma),
        mu <- aGroup[Group_num] + bFemur * femur_length_c,
        bFemur ~ dnorm(0, 5),
        aGroup[Group_num] ~ normal(0, 2),
        sigma ~ cauchy(0, 2)
      ),
      data = dat_s,
      iter = iter,
      warmup = warmup,
      chains = chains,
      cores = cores,
      WAIC = FALSE,
      control = list(adapt_delta = 0.95)
    )
    
    if (plot_diagnostics) plot(fm)
    print(precis(fm, depth = 2, digits = 4))
    
    out <- save_post_est(y, mean(M_tmp$y, na.rm = TRUE), fm, out, r = 2 + ii)
    save(fm, file = paste0(y, ".Rda"))
  } else {
    load(paste0(y, ".Rda"))
  }
}

kable(out)
write_csv(out, path = "Histomorphometry_Output.csv")

## ------------------------------------------------------------------------
font_size <- 10
my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = font_size, face = "bold"),
  legend.text = element_text(size = font_size - 1),
  legend.title = element_text(size = font_size, face = "bold"),
  plot.title = element_text(size = font_size + 1),
  strip.text = element_text(size = font_size - 1),
  axis.text.x = element_text(angle = 45, hjust = 1)
)

pFemLen <- ggplot(M, aes(x = Group, y = femur_length)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("Femur length (mm)")) +
  my_theme


pML <- ggplot(M, aes(x = Group, y = ML_diam)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("ML diameter (mm)")) +
  my_theme


pAP <- ggplot(M, aes(x = Group, y = AP_diam)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("AP diameter (mm)")) +
  my_theme


pCSA <- ggplot(M, aes(x = Group, y = CSA)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("CSA (mm^2)")) +
  my_theme


pImin <- ggplot(M, aes(x = Group, y = Imin)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("I_{min} (mm^4)")) +
  my_theme


pImax <- ggplot(M, aes(x = Group, y = Imax)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("I_{max} (mm^4)")) +
  my_theme


pJ <- ggplot(M, aes(x = Group, y = J)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("J (mm^4)")) +
  my_theme


pElip <- ggplot(M, aes(x = Group, y = Elip)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("Ellipticity")) +
  my_theme


pZmin <- ggplot(M, aes(x = Group, y = Zmin)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("Zmin (mm^4)")) +
  my_theme


pZmax <- ggplot(M, aes(x = Group, y = Zmax)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("Zmax (mm^4)")) +
  my_theme


pZpol <- ggplot(M, aes(x = Group, y = Zpol)) + 
  geom_point(position = position_jitter(width = 0.1), size = 1) +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = TeX("Zpol (mm^4)")) +
  my_theme


## ------------------------------------------------------------------------
p <- plot_grid(pFemLen, pAP, pML, pCSA, pImin, pImax, ncol = 3, nrow = 2)
save_plot(file = "../../Figures and tables/Figure_3.pdf", p,
          base_height = 5, base_width = 181/25.4)

## ------------------------------------------------------------------------

N %>% 
  dplyr::select("Group", "femur_length", vars) %>% 
  group_by(Group) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  write_csv(., path = "Means.csv")
