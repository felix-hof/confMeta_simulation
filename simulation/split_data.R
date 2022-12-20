############################################################
# This script divides the simulation output (a large list) #
# into three data frames (one for each figure directory)   #
# and saves these as .rds files                            #
############################################################

# load dplyr library
library(dplyr)

# load simulation output
load("RData/simulate_all.RData")

# Prepare data (means) ---
out_meanplots <- out %>%
  bind_rows() %>%
  filter(grepl("_mean$", measure),
         !grepl("^gammaMin", measure),
         !grepl("two sided", method)) %>%
  mutate(method = ifelse(
    grepl("REML", method),
    gsub("REML", "Random Effects, default, REML", method), method),
         measure = gsub("_mean$", "", measure))
saveRDS(out_meanplots, file = "RData/meanplots.rds")

# Prepare data (gammaMin) ----
out_gammamin <- out %>%
  bind_rows() %>%
  filter(grepl("^gammaMin", measure),
         grepl("Harmonic Mean", method),
         !grepl("two sided", method)) %>%
  distinct() %>% # throw out the doubly included gammaMin_mean entries
  mutate(measure = gsub("^gammaMin_", "", measure),
         measure = case_when(measure == "min" ~ "Minimum",
                             measure == "firstQuart" ~ "1. Quartile",
                             measure == "mean" ~ "Mean",
                             measure == "median" ~ "Median",
                             measure == "thirdQuart" ~ "3. Quartile",
                             measure == "max" ~ "Maximum"),
         measure = ordered(measure, levels = c("Maximum", "3. Quartile", "Median", "Mean", "1. Quartile", "Minimum"))) %>%
  arrange(k, dist, bias, large, heterogeneity, I2, method, measure)
saveRDS(out_gammamin, file = "RData/min_pH.rds")

# Prepare data (frequency) ----
out_n <- out %>%
  bind_rows() %>%
  filter(grepl("^n", measure),
         !grepl("_mean$", measure)) %>%
  mutate(measure = gsub("n_", "", measure),
         measure = gsub("gt", "> ", measure),
         measure = ordered(measure, levels = rev(c("> 9", as.character(9:1)))),
         value = value / 1e4)
saveRDS(out_n, file = "RData/rel_freq.rds")
