# ==============================================================================
# 01 - Raw data processing (adult and juvenile census data -> species-level data)
# ==============================================================================
# NOTE ON REPRODUCIBILITY:
# This script documents how the two species-level tables shipped in this
# repository (processed-data/adult_species_level.csv and
# processed-data/juvenile_species_level.csv) were built from the original,
# individual-level forest census data (raw-data/adult_data_name_corrections.csv
# and raw-data/juvenile_data_name_corrections.csv).
#
# The raw census files and the individual/stem-level outputs of this script
# (adult_stems_level.csv, juvenile_stems_level.csv) are NOT distributed in this
# repository. The species-level tables in processed-data/ are the starting
# point for reproducing the analysis - start from script 02.
#
# This script is kept for transparency only and will not run without the
# original raw-data files.
# ==============================================================================

library(dplyr)

# ------------------------------------------------------------------
# Adults ----
# ------------------------------------------------------------------
rm(list = ls()) # clean the environment

adult <- readr::read_csv("raw-data/adult_data_name_corrections.csv")

adult_plots <- adult %>%
  select(TreeID, Tag.No, Stem.Group.ID, Main.Stem.Tag, Plot.Code, Sub.Plot.T1, Sub.Plot.T2,
         Family, SpeciesID, Species, Census.No, Census.Date, F1, F2, D4, POM, Comments) %>%
  filter(Plot.Code %in% c("PRM-01", "PRM-02", "PRM-03", "EEA-01", "EEA-02", "EEA-03",
                           "PNA-01", "PNA-02", "PNA-03")) %>%
  mutate(D.cm = D4 / 10, .keep = c("unused"), .after = D4) %>%
  mutate(basal.area = (pi * (D.cm / 2)^2), .after = D.cm)

# census 1 ----
t1_adult <- adult_plots %>%
  filter(Census.No == 1)

# multi-stemmed individuals
multi_t1_adult <- t1_adult %>%
  filter(!is.na(Stem.Group.ID == TRUE)) %>%
  group_by(Stem.Group.ID, Plot.Code, Sub.Plot.T1, Family, Species) %>%
  summarise(Census.Date_1 = mean(Census.Date),
            D_1 = max(D.cm),
            POM_1 = max(POM),
            basal.area_1 = sum(basal.area),
            n_stems_1 = n()) %>%
  rename(TreeID = Stem.Group.ID)

# single-stemmed individuals
one_1_adult <- t1_adult %>%
  filter(!is.na(Stem.Group.ID) == FALSE) %>%
  select(TreeID, Plot.Code, Sub.Plot.T1, Family, Species, Census.Date, D.cm, POM, basal.area) %>%
  rename(Census.Date_1 = Census.Date,
         D_1 = D.cm,
         POM_1 = POM,
         basal.area_1 = basal.area)

total_t1 <- bind_rows(multi_t1_adult, one_1_adult) # 1551 individuals

# census 2 ----
t2_adult <- adult_plots %>%
  filter(Census.No == 2)

multi_t2_adult <- t2_adult %>%
  filter(!is.na(Stem.Group.ID == TRUE)) %>%
  group_by(Stem.Group.ID, Plot.Code, Sub.Plot.T1, Family, Species) %>%
  summarise(Census.Date_2 = mean(Census.Date),
            D_2 = max(D.cm),
            POM_2 = max(POM),
            basal.area_2 = sum(basal.area),
            n_stems_2 = n()) %>%
  rename(TreeID = Stem.Group.ID)

one_2_adult <- t2_adult %>%
  filter(!is.na(Stem.Group.ID) == FALSE) %>%
  select(TreeID, Plot.Code, Sub.Plot.T1, Family, Species, Census.Date, D.cm, POM, basal.area) %>%
  rename(Census.Date_2 = Census.Date,
         POM_2 = POM,
         D_2 = D.cm,
         basal.area_2 = basal.area)

total_t2 <- bind_rows(multi_t2_adult, one_2_adult) # 1712 individuals

# species with 10 or more individuals ----
species_selection <- total_t1 %>%
  group_by(Species) %>%
  summarise(TreeID = n()) %>%
  filter(TreeID >= 10)

# average census time per plot ----
total_join <- inner_join(total_t1, total_t2)

time <- total_join %>%
  group_by(Plot.Code) %>%
  summarise(t1 = mean(Census.Date_1),
            t2 = mean(Census.Date_2))

time_mean <- time %>%
  mutate(t1_mean = mean(time$t1), t2_mean = mean(time$t2), t = t2_mean - t1_mean) %>%
  mutate(interval = t2 - t1)

# all stems, both censuses ----
full_stems <- full_join(total_t1, total_t2)

# stem-level table (not distributed in this repository)
write.csv(full_stems, "processed-data/adult_stems_level.csv")

# classify surviving, recruited and dead individuals ----
sp.survivals <- full_stems %>%
  filter(!is.na(D_1), !is.na(D_2), D_2 != 0)

sp.recruits <- subset(full_stems, is.na(D_1))

sp.deads <- full_stems %>%
  filter(!is.na(D_1) & (is.na(D_2) | D_2 == 0))

# species-level counts of survivors, recruits and deaths ----
Ns <- sp.survivals %>%
  group_by(Species) %>%
  summarise(Ns = n())

Nr <- sp.recruits %>%
  group_by(Species) %>%
  summarise(Nr = n())

Nm <- sp.deads %>%
  group_by(Species) %>%
  summarise(Nm = n())

# combine counts into a single table
x1 <- left_join(Ns, Nr, by = "Species")
x2 <- left_join(x1, Nm, by = "Species")

# keep species with 10 or more individuals
data_final <- left_join(species_selection, x2, by = "Species")
data_final <- data_final[, -2] # drop the TreeID count column, no longer needed

# replace NA with 0 (species with no recruits/deaths in a given census)
data_final[is.na(data_final)] <- 0

# save the species-level table distributed in processed-data/
write.csv(data_final, "processed-data/adult_species_level.csv")


# ------------------------------------------------------------------
# Juveniles ----
# ------------------------------------------------------------------
rm(list = ls()) # clean the environment

juvenile <- readr::read_csv("raw-data/juvenile_data_name_corrections.csv")

# restrict to the species retained in the adult dataset
# NOTE: species-info/species_list.csv is itself produced by script 03
# (03_biogeographic_origin_and_leaf_habit.R), after the adult vital rates
# have been computed. In the original workflow this script was therefore
# re-run once species_list.csv already existed from a previous iteration.
species <- read.csv("species-info/species_list.csv")
species <- species[, -3]

juvenile_plots <- juvenile %>%
  select(TreeID, Tag.No, Stem.Group.ID, Main.Stem.Tag, Plot.Code, Sub.Plot.T1, Sub.Plot.T2,
         Family, SpeciesID, Species, Census.No, Census.Date, F1, F2, D4, Comments) %>%
  filter(Plot.Code %in% c("PRM-01", "PRM-02", "PRM-03", "EEA-01", "EEA-02", "EEA-03",
                           "PNA-01", "PNA-02", "PNA-03")) %>%
  mutate(D.cm = D4 / 10, .keep = c("unused"), .after = D4) %>%
  mutate(basal.area = (pi * (D.cm / 2)^2), .after = D.cm)

juvenile_selection <- left_join(species, juvenile_plots, by = "Species")

# individuals that recruited into the adult stage
promoted_to_adult <- juvenile_selection %>%
  filter(str_detect(Comments, "(?i)recrutou")) %>%
  mutate(D.cm = 1, basal.area = 1)  # flag D as 1 to signal the individual did not die

# update the original table with the flagged records
juvenile_selection <- juvenile_selection %>%
  rows_update(promoted_to_adult, by = c("TreeID", "Census.No"))

# census 1 ----
t1 <- juvenile_selection %>%
  filter(Census.No == 1)

multi_t1 <- t1 %>%
  filter(Stem.Group.ID != 0) %>%
  group_by(Stem.Group.ID, Plot.Code, Family, Species) %>%
  summarise(Census.Date_1 = mean(Census.Date),
            D_1 = sqrt(sum(D.cm^2)),
            n_stems_1 = n()) %>%
  rename(TreeID = Stem.Group.ID)

stem_t1 <- t1 %>%
  filter(Stem.Group.ID == 0) %>%
  select(TreeID, Plot.Code, Family, Species, Census.Date, D.cm) %>%
  rename(Census.Date_1 = Census.Date,
         D_1 = D.cm)

total_t1 <- bind_rows(multi_t1, stem_t1) # 1011 individuals

# census 2 ----
t2 <- juvenile_selection %>%
  filter(Census.No == 2)

multi_t2 <- t2 %>%
  filter(Stem.Group.ID != 0) %>%
  group_by(Stem.Group.ID, Plot.Code, Family, Species) %>%
  summarise(Census.Date_2 = mean(Census.Date),
            D_2 = sqrt(sum(D.cm^2)),
            n_stems_2 = n()) %>%
  rename(TreeID = Stem.Group.ID)

stem_t2 <- t2 %>%
  filter(Stem.Group.ID == 0) %>%
  select(TreeID, Plot.Code, Family, Species, Census.Date, D.cm) %>%
  rename(Census.Date_2 = Census.Date,
         D_2 = D.cm)

total_t2 <- bind_rows(multi_t2, stem_t2) # 1258 individuals

# combine both censuses ----
all_stems <- full_join(total_t1, total_t2)

# stem-level table (not distributed in this repository)
write.csv(all_stems, "processed-data/juvenile_stems_level.csv")

# classify surviving, recruited and dead individuals ----
sp.survivals <- all_stems %>%
  filter(!is.na(D_1), !is.na(D_2), D_2 != 0)

sp.recruits <- all_stems %>%
  filter(is.na(D_1))

sp.deads <- all_stems %>%
  filter(!is.na(D_1) & (is.na(D_2) | D_2 == 0))

# all_stems = survivors (839) + recruits (255) + deaths (173) = 1267 individuals

# species-level counts of survivors, deaths and recruits ----
Ns <- sp.survivals %>%
  group_by(Species) %>%
  summarise(Ns = n())

Nr <- sp.recruits %>%
  group_by(Species) %>%
  summarise(Nr = n())

Nm <- sp.deads %>%
  group_by(Species) %>%
  summarise(Nm = n())

# combine counts into a single table
x1 <- left_join(Ns, Nm, by = "Species")
data_final <- left_join(x1, Nr, by = "Species")

# replace NA with 0 (species with no recruits/deaths in a given census)
data_final[is.na(data_final)] <- 0

# save the species-level table distributed in processed-data/
write.csv(data_final, "processed-data/juvenile_species_level.csv")
