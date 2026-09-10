# ==============================================================================
# 03 - Combine vital rates with biogeographic origin and leaf habit
# ==============================================================================
# Input:  vital-rates/adult_dynamics.csv
#         vital-rates/juvenile_dynamics.csv
#         vital-rates/species_growth_rate.csv
#         species-info/biogeo_origin_leaf_habit.csv
# Output: species-info/species_list.csv
#         data-analysis/ontogeny.csv
#         data-analysis/growth_population.csv
# ==============================================================================

library(dplyr)

rm(list = ls()) # clean the environment
dir.create("data-analysis", showWarnings = FALSE)

# datasets ----
adult <- read.csv("vital-rates/adult_dynamics.csv"); adult <- adult[, -1]
juvenile <- read.csv("vital-rates/juvenile_dynamics.csv"); juvenile <- juvenile[, -1]
growth_rate <- read.csv("vital-rates/species_growth_rate.csv"); growth_rate <- growth_rate[, -1]
origin <- read.csv("species-info/biogeo_origin_leaf_habit.csv")

# manual correction of biogeographic origin (07/04/2025)
origin <- origin %>%
  mutate(Origin = if_else(Species == "Myrcia glomerata", "Temperate", Origin))

# adults + origin ----
adult_origin <- merge(adult, origin, by = "Species")
adult_origin <- mutate(adult_origin, Stage = "Adult")

# juveniles + origin ----
juvenile_origin <- merge(juvenile, origin, by = "Species")
juvenile_origin <- mutate(juvenile_origin, Stage = "Juvenile")

# species growth rate + origin ----
growth_origin <- merge(growth_rate, origin, by = "Species")

# keep the same set of species in both ontogenetic stages ----
missing_in_juveniles <- setdiff(adult_origin$Species, juvenile_origin$Species)
missing_in_juveniles

adult_origin <- adult_origin %>%
  filter(Species != "Eugenia oeidocarpa") %>%
  filter(Species != "Laplacea acutifolia") %>%
  filter(Species != "Ocotea pulchella")

# updated species list (used as input by script 01 for the juvenile dataset) ----
species_list <- adult_origin %>%
  select(Species, Origin, Deciduousness)
write.csv(species_list, "species-info/species_list.csv")

# combine adults and juveniles into a single data frame ----
ontogeny <- rbind(juvenile_origin, adult_origin)

# save tables used in the statistical analysis (script 04) ----
write.csv(ontogeny, "data-analysis/ontogeny.csv")
write.csv(growth_origin, "data-analysis/growth_population.csv")
