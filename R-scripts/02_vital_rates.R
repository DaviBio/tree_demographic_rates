# ==============================================================================
# 02 - Vital rates (recruitment, mortality, turnover and population growth)
# ==============================================================================
# Input:  processed-data/adult_species_level.csv
#         processed-data/juvenile_species_level.csv
# Output: vital-rates/adult_dynamics.csv
#         vital-rates/juvenile_dynamics.csv
#         vital-rates/species_growth_rate.csv
#
# Formulas:
#   annual recruitment rate (Raf): 1 - (Ns/N2)^(1/t)
#   annual mortality rate (Ma):    1 - (Ns/N1)^(1/t)
#   turnover:                      (Raf + Ma) / 2
#   finite population growth rate: (N2/N1)^(1/t)
# where t is the average census interval (in years), N1 the number of
# individuals at the start of the interval, N2 the number at the end, and Ns
# the number of survivors.
# ==============================================================================

library(dplyr)

rm(list = ls()) # clean the environment
dir.create("vital-rates", showWarnings = FALSE)

census_interval <- 5.594278 # average interval between censuses, in years

# dataset ----
species_adult <- read.csv("processed-data/adult_species_level.csv")
species_juvenile <- read.csv("processed-data/juvenile_species_level.csv")

# adults -------
adult_dynamics <- species_adult %>%
  mutate(N1 = Nm + Ns) %>%
  mutate(N2 = Ns + Nr) %>%
  select(Species, N1, Ns, N2) %>%
  mutate(t = census_interval) %>%
  mutate(Raf = 1 - (Ns / N2)^(1 / t)) %>%
  mutate(Ma = 1 - (Ns / N1)^(1 / t)) %>%
  mutate(turnover = (Ma + Raf) / 2) %>%
  mutate(growth_finite = (N2 / N1)^(1 / t))

# rates expressed as percentages
percentual_adult <- adult_dynamics %>%
  mutate(Raf = Raf * 100) %>%
  mutate(Ma = Ma * 100) %>%
  mutate(turnover = turnover * 100)

write.csv(percentual_adult, "vital-rates/adult_dynamics.csv")

# juveniles -----
juvenile_dynamics <- species_juvenile %>%
  mutate(N1 = Nm + Ns) %>%
  mutate(N2 = Ns + Nr) %>%
  select(Species, N1, Ns, N2) %>%
  mutate(t = census_interval) %>%
  mutate(Raf = 1 - (Ns / N2)^(1 / t)) %>%
  mutate(Ma = 1 - (Ns / N1)^(1 / t)) %>%
  mutate(turnover = (Ma + Raf) / 2) %>%
  mutate(growth_finite = (N2 / N1)^(1 / t))

# rates expressed as percentages
percentual_juvenile <- juvenile_dynamics %>%
  mutate(Raf = Raf * 100) %>%
  mutate(Ma = Ma * 100) %>%
  mutate(turnover = turnover * 100)

write.csv(percentual_juvenile, "vital-rates/juvenile_dynamics.csv")

# combined species-level dynamics -----
# (species = juveniles + adults pooled, no distinction by ontogenetic stage)
combined <- merge(adult_dynamics, juvenile_dynamics, by = "Species")

combined_dynamics <- combined %>%
  mutate(N1 = N1.x + N1.y) %>%
  mutate(Ns = Ns.x + Ns.y) %>%
  mutate(N2 = N2.x + N2.y) %>%
  rename(t = t.x) %>%
  select(Species, N1, Ns, N2, t) %>%
  mutate(Raf = 1 - (Ns / N2)^(1 / t)) %>%
  mutate(Ma = 1 - (Ns / N1)^(1 / t)) %>%
  mutate(turnover = (Ma + Raf) / 2) %>%
  mutate(growth_finite = (N2 / N1)^(1 / t))

# rates expressed as percentages
species_growth_rate <- combined_dynamics %>%
  mutate(Raf = Raf * 100) %>%
  mutate(Ma = Ma * 100) %>%
  mutate(turnover = turnover * 100)

write.csv(species_growth_rate, "vital-rates/species_growth_rate.csv")
