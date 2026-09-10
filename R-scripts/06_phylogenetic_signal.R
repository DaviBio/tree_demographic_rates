# ==============================================================================
# 06 - Phylogenetic signal in demographic rates
# ==============================================================================
# 1) build the phylogenetic tree for the studied species
# 2) Fritz & Purvis' D statistic for the categorical predictors
# 3) Moran's I test on the residuals of the ordbeta GLM
# ==============================================================================

library(tidyverse)
library(ape)
library(caper)
library(DHARMa)
library(glmmTMB)
library(remotes)
remotes::install_github("jinyizju/V.PhyloMaker2")
library(V.PhyloMaker2)
library(patchwork)

# load the data -----
# the species column must be formatted as "Genus_species".

# convert rates to proportions for the Gamma/ordbeta GLM below
ontogeny <- read.csv("data-analysis/ontogeny.csv") %>%
  mutate(Raf = Raf / 100) %>%
  mutate(Ma = Ma / 100) %>%
  mutate(turnover = turnover / 100) %>%
  rename(species = Species)

# split genus and epithet
species <- ontogeny %>%
  dplyr::select(species, Origin, Deciduousness, Raf, Ma, turnover) %>%
  distinct(species, .keep_all = TRUE) %>%  # keep one row per species
  separate(species, into = c("genus", "epithet"), sep = " ", remove = TRUE) %>%
  unite(species, genus:epithet, sep = "_")

# convert categorical attributes (origin and leaf habit) to factors
df_analysis <- species %>%
  mutate(
    species = as.character(species),
    Origin = as.factor(Origin),
    Deciduousness = as.factor(Deciduousness))

# 1. Build the phylogenetic tree ----
# species list in the format required by V.PhyloMaker2
species_list <- df_analysis %>%
  dplyr::select(species) %>%
  separate(species, into = c("Genus", "Species_epithet"), sep = "_", remove = FALSE) %>%
  mutate(Family = NA) %>% # V.PhyloMaker looks up the family automatically when left NA
  filter(species != "Siphoneugena_reitzii")

# add a placeholder row to force the relationship with Myrciaria floribunda
# (Genus set to "Myrciaria" instead of "Siphoneugena")
siphoneugena_placeholder <- data.frame(
  species = "Siphoneugena_reitzii",
  Genus = "Myrciaria",
  Species_epithet = "reitzii",
  Family = NA
)

species_list_modified <- bind_rows(species_list, siphoneugena_placeholder)

# generate the phylogeny
tree_output <- V.PhyloMaker2::phylo.maker(sp.list = species_list_modified, scenarios = "S3")
tree <- tree_output$scenario.3

# correct the tip label back to the actual species name
tree$tip.label[tree$tip.label == "Myrciaria_reitzii"] <- "Siphoneugena_reitzii"

# quick visual check of the tree
plot(tree, no.margin = TRUE, cex = 0.7)

# 2. Phylogenetic signal test for categorical traits (Fritz & Purvis' D) ----
# the 'species' column must match the tree tip labels exactly.

# remove internal node labels to avoid duplicated names
tree$node.label <- NULL

comparative_data <- comparative.data(phy = tree,
                                      data = df_analysis,
                                      names.col = "species",
                                      vcv = TRUE)

# test for biogeographic origin
phylo_signal_origin <- phylo.d(data = comparative_data, binvar = Origin)
summary(phylo_signal_origin)

# H0: p > 0.05 (0.233) -> cannot reject that Origin is randomly distributed
# H1: p > 0.05 (0.104) -> cannot reject that Origin follows a Brownian model
# conclusion: no strong statistical evidence of phylogenetic signal for biogeographic origin

# test for leaf habit
phylo_signal_leaf_habit <- phylo.d(data = comparative_data, binvar = Deciduousness)
summary(phylo_signal_leaf_habit)

# H0: p > 0.05 (0.133) -> cannot reject that leaf habit is randomly distributed
# H1: p > 0.05 (0.191) -> cannot reject that leaf habit follows a Brownian model
# conclusion: no strong statistical evidence of phylogenetic signal for leaf habit


# 3. Phylogenetic signal in model residuals (Moran's I) ------

# extract randomized quantile residuals from the selected (AIC-best) model for each rate
rec <- glmmTMB(Raf ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = ontogeny)
mort <- glmmTMB(Ma ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = ontogeny)
turn <- glmmTMB(turnover ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = ontogeny)

# extract residuals
dharma_residuals <- simulateResiduals(rec)
res_values <- dharma_residuals$scaledResiduals # 60 residuals

# 60 residuals need to be aligned with 30 species (each species appears twice,
# once per Stage). Solution: average the residuals per species.
residuals_by_species <- tibble(
  species = ontogeny$species,
  Residuals = dharma_residuals$scaledResiduals
) %>%
  group_by(species) %>%
  summarise(Residual_Mean = mean(Residuals))

# species associated with the model observations, in model order
model_species <- df_analysis$species

# build the phylogenetic weight matrix (inverse distance)
phylo_distance <- cophenetic(tree)

# reorder the distance matrix to match the model's species order
phylo_distance_ordered <- phylo_distance[model_species, model_species]

# inverse distance (evolutionary isolation)
inverse_distance_matrix <- 1 / phylo_distance_ordered
diag(inverse_distance_matrix) <- 0 # remove self-distance (Inf)

# replace any Inf values in case of identical/synonymous species in the tree
inverse_distance_matrix[is.infinite(inverse_distance_matrix)] <- 0

# 3.1 Moran's I on the residuals ----
moran_test_residuals <- Moran.I(residuals_by_species$Residual_Mean, weight = inverse_distance_matrix)

print(moran_test_residuals)

# figure and results for the supplementary information ------
# 1. function to extract orthogonal branch coordinates
extract_orthogonal_coordinates <- function(phy) {
  pdf(NULL)
  # force a rectangular layout on the (invisible) plot
  plot(phy, type = "phylogram", plot = FALSE)
  coord <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  dev.off()

  node_coords <- tibble(
    node = 1:(phy$Nnode + length(phy$tip.label)),
    x = coord$xx,
    y = coord$yy,
    isTip = node <= length(phy$tip.label),
    species = c(phy$tip.label, rep(NA, phy$Nnode))
  )

  # build orthogonal (right-angled) segments for each branch
  orthogonal_lines <- tibble()
  for (i in 1:nrow(phy$edge)) {
    parent <- phy$edge[i, 1]
    child <- phy$edge[i, 2]

    x0 <- node_coords$x[parent]
    y0 <- node_coords$y[parent]
    x1 <- node_coords$x[child]
    y1 <- node_coords$y[child]

    # in a rectangular cladogram, the branch leaves the parent (x0,y0),
    # goes vertically to the child's height (x0,y1) and then horizontally (x1,y1)
    branch_segment <- tibble(
      x = c(x0, x0, x1),
      y = c(y0, y1, y1),
      group = i
    )
    orthogonal_lines <- bind_rows(orthogonal_lines, branch_segment)
  }

  return(list(points = node_coords, branches = orthogonal_lines))
}


# extract coordinates for plotting
tree_coords <- extract_orthogonal_coordinates(tree)

tree_plot_data <- tree_coords$points %>%
  filter(isTip) %>%
  left_join(df_analysis, by = "species")


# panel A) biogeographic origin
origin_panel <- ggplot() +
  # geom_path draws the orthogonal connections using the 'group' of each branch
  geom_path(data = tree_coords$branches, aes(x = x, y = y, group = group), color = "#2c3e50", linewidth = 0.6) +
  geom_point(data = tree_plot_data, aes(x = x, y = y, color = Origin), size = 3) +
  geom_text(data = tree_plot_data, aes(x = x, y = y, label = species), hjust = 0, nudge_x = 7, size = 3, fontface = "italic") +
  scale_color_manual(values = c("Tropical" = "#EBC821", "Temperate" = "#692E8C")) +
  labs(title = "A) Biogeographic Origin", color = "Origin") +
  xlim(0, max(tree_coords$points$x) * 1.8) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "bottom")


# panel B) leaf habit
leaf_habit_panel <- ggplot() +
  geom_path(data = tree_coords$branches, aes(x = x, y = y, group = group), color = "#2c3e50", linewidth = 0.6) +
  geom_point(data = tree_plot_data, aes(x = x, y = y, color = Deciduousness), size = 3) +
  geom_text(data = tree_plot_data, aes(x = x, y = y, label = species), hjust = 0, nudge_x = 7, size = 3, fontface = "italic") +
  scale_color_manual(values = c("Deciduous" = "#9E9454", "Evergreen" = "#3F92B5")) +
  labs(title = "B) Leaf Habit", color = "Habit") +
  xlim(0, max(tree_coords$points$x) * 1.8) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_blank(), legend.position = "bottom")

# combined panel
phylo_panel <- origin_panel | leaf_habit_panel
print(phylo_panel)
