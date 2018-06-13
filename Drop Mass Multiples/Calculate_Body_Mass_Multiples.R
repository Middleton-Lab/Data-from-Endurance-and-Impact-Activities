library(readxl)
library(tidyverse)
library(ggplot2)
library(cowplot)

font_size <- 10
my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size, face = "bold"),
  legend.text = element_text(size = font_size - 1),
  legend.title = element_text(size = font_size, face = "bold"),
  plot.title = element_text(size = font_size + 1),
  strip.text = element_text(size = font_size - 1)
)

## Make graphs of force traces

LJ_Reader <- function(filename, zero_start = 1, zero_end = 500){
  M <- read.table(filename, skip = 6, header = TRUE)[, 1:2]
  zero <- mean(M$v0[zero_start:zero_end])
  M$v0 <- M$v0 - zero
  M$v0 <- -(M$v0)
  plot(M$v0, type = 'l', main = filename)
  return(M)
}

## Extract data for spreadsheet
C <- read_excel("Drop trace start & end.xlsx")
C <- C %>% filter(is.na(checkme))

file_list <- C$F_name

out <- data.frame(file = file_list,
                  Mouse_ID = numeric(length = length(file_list)),
                  Trial = numeric(length = length(file_list)),
                  Mass = numeric(length = length(file_list)),
                  Max = numeric(length = length(file_list)),
                  Multiples = numeric(length = length(file_list)))
                 
for (i in 1:length(file_list)) {
  infile <- paste0(file_list[i], ".dat")
  start_pt <- C$zero_st[C$F_name == strsplit(infile, ".", fixed = TRUE)[[1]][1]]
  end_pt <- C$end_zero[C$F_name == strsplit(infile, ".", fixed = TRUE)[[1]][1]]
  sit_st <- C$sit_st[C$F_name == strsplit(infile, ".", fixed = TRUE)[[1]][1]]
  sit_end <- C$sit_end[C$F_name == strsplit(infile, ".", fixed = TRUE)[[1]][1]]
  
  M <- LJ_Reader(infile, zero_start = start_pt, zero_end = end_pt)
  out$Mouse_ID[i] <- strsplit(infile, "_")[[1]][1]
  out$Trial[i] <- sub(".dat", "", strsplit(infile, "_")[[1]][2])
  Mass <- mean(M$v0[sit_st:sit_end])
  Max <- max(M$v0)
  out$Mass[i] <- Mass
  out$Max[i] <- Max
  out$Multiples[i] <- Max / Mass
}

write_csv(out, path = "Body_Mass_Multiples.csv")

## Plot distribution of impact forces

out <- read_csv("Body_Mass_Multiples.csv")

ggplot(out, aes(x = Multiples)) +
  geom_histogram(binwidth = 1)

out %>% 
  group_by(Mouse_ID) %>% 
  summarize(mean_Mb = mean(Multiples))

mass_multiples <- ggplot(out, aes(x = Mouse_ID, y = Multiples)) + 
  geom_point(position = position_jitter(width = 0.1), size = 2) +
  geom_hline(yintercept = mean(out$Multiples), linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1,
               color = "red", size = 1) +
  labs(y = "Multiples of Body Mass", x = "Mouse ID") +
  my_theme
ggsave("../../Figures and tables/Figure_SI1.png", plot = mass_multiples,
       width = 180 / 25.4, height = 4)
