library(dplyr)

load("RData/simulate_all.RData")

View(out)

# Do we expect NAs in output values?
out %>% 
  filter(is.na(value)) %>% 
  View()
# All NAs happen with method "Harmonic Mean two sided"


# Now Harmonic mean multiplicative looks pretty much the same as Harmonic mean additive
out %>% 
  filter(method %in% c("Harmonic Mean Additive", "Harmonic Mean Multiplicative")) %>% 
  group_by(across(c(-value, -method))) %>% 
  summarise(is_equal = value[method == "Harmonic Mean Additive"] == value[method == "Harmonic Mean Multiplicative"],
            .groups = "drop")

# See whether this is the case also for larger trials
out %>% 
  filter(method %in% c("Harmonic Mean Additive", "Harmonic Mean Multiplicative")) %>% 
  group_by(across(c(-value, -method))) %>% 
  summarise(is_equal = value[method == "Harmonic Mean Additive"] == value[method == "Harmonic Mean Multiplicative"],
            .groups = "drop") %>% 
  filter(measure == "coverage_mean") %>% 
  View()
