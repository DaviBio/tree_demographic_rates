# ==============================================================================
# 04 - Generalized linear models (GLM)
# ==============================================================================
# Tests whether vital rates (recruitment, mortality, turnover) and finite
# population growth differ with biogeographic origin, leaf habit
# (deciduousness) and ontogenetic stage (adult vs. juvenile).
#
# Input:  data-analysis/ontogeny.csv
#         data-analysis/growth_population.csv
# Output: supplementary/model_estimates_table.csv (fixed-effect estimates
#         for the selected models)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(betareg)
library(glmmTMB)
library(performance)
library(DHARMa)

rm(list = ls()) # clean the environment
dir.create("supplementary", showWarnings = FALSE)

plot_theme <- theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        text = element_text(size = 15),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
dodge_position <- position_dodge(0.3)

# dataset ----
ontogeny <- read.csv("data-analysis/ontogeny.csv")
growth <- read.csv("data-analysis/growth_population.csv"); growth <- growth[, -1]

# Are vital rates different between adults and juveniles? ----

# 1) data exploration ----
head(ontogeny)
str(ontogeny)

df <- ontogeny %>%
  select(Species, Stage, Origin, Deciduousness, Raf, Ma, turnover) %>%
  mutate(Stage = as.factor(Stage),
         Origin = as.factor(Origin),
         Deciduousness = as.factor(Deciduousness))
str(df)

means_by_stage <- df %>%
  group_by(Stage) %>%
  summarise(mean = mean(Ma),
            sd = sd(Ma))

raf_density <- df %>%
  ggplot(aes(Raf)) +
  geom_density(colour = "cyan4", fill = "cyan4", alpha = 0.4) +
  theme(legend.position = "none") +
  labs(x = "Recruitment Rate (%)", y = "Density") +
  facet_grid(~Stage) +
  plot_theme

ma_density <- df %>%
  ggplot(aes(Ma)) +
  geom_density(colour = "cyan4", fill = "cyan4", alpha = 0.4) +
  theme(legend.position = "none") +
  labs(x = "Mortality Rate (%)", y = "Density") +
  facet_grid(~Stage) +
  plot_theme

turnover_density <- df %>%
  ggplot(aes(turnover)) +
  geom_density(colour = "cyan4", fill = "cyan4", alpha = 0.4) +
  theme(legend.position = "none") +
  labs(x = "Turnover Rate (%)", y = "Density") +
  facet_grid(~Stage) +
  plot_theme

raf_density / ma_density / turnover_density

# simple comparison
ggplot(df, aes(x = Stage, y = Raf, fill = Stage)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Stage", y = "Recruitment Rate (%)") +
  plot_theme

# zero values ----
df_zero <- df %>%
  filter(if_any(c(Raf, Ma, turnover), ~ . == 0))

juvenile_zero_counts <- df_zero %>%
  filter(Stage == "Juvenile") %>%
  summarise(across(c(Raf, Ma, turnover),
                    ~ sum(.x == 0, na.rm = TRUE),
                    .names = "zeros_{col}"))

adult_zero_counts <- df_zero %>%
  filter(Stage == "Adult") %>%
  summarise(across(c(Raf, Ma, turnover),
                    ~ sum(.x == 0, na.rm = TRUE),
                    .names = "zeros_{col}"))

zero_counts_combined <- bind_rows(juvenile_zero_counts, adult_zero_counts)

stages <- c("Juvenile", "Adult")
zero_counts_by_stage <- data.frame(Stage = stages,
                                    n0_raf = zero_counts_combined$zeros_Raf,
                                    n0_ma = zero_counts_combined$zeros_Ma,
                                    n0_turnover = zero_counts_combined$zeros_turnover)

# juveniles and adults show different vital-rate distribution patterns

# --- Beta GLM ---
# rescale rates to the 0-1 interval
df_proportions <- df %>%
  mutate(Raf = Raf / 100) %>%
  mutate(Ma = Ma / 100) %>%
  mutate(turnover = turnover / 100)

# candidate models ----
m1 <- glmmTMB(Raf ~ Origin * Stage * Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)
m1.1 <- glmmTMB(Ma ~ Origin * Stage * Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)
m1.2 <- glmmTMB(turnover ~ Origin * Stage * Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)
check_model(m1.1)
simulateResiduals(m1.1, plot = TRUE) # no issues detected
r2(m1)
summary(m1)
visreg::visreg(m1)

m2 <- glmmTMB(Raf ~ Origin * Stage + Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)
m2.1 <- glmmTMB(Ma ~ Origin * Stage + Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)
m2.2 <- glmmTMB(turnover ~ Origin * Stage + Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)
check_model(m2)
summary(m2.1)
visreg::visreg(m2)

# most parsimonious model
m3 <- glmmTMB(Raf ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)
m3.1 <- glmmTMB(Ma ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)
m3.2 <- glmmTMB(turnover ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = df_proportions)

check_model(m3.2)
summary(m3)
visreg::visreg(m3.1)

m4 <- glmmTMB(Raf ~ Origin * Deciduousness + Stage, family = ordbeta(link = "logit"), data = df_proportions)
m4.1 <- glmmTMB(Ma ~ Origin * Deciduousness + Stage, family = ordbeta(link = "logit"), data = df_proportions)
m4.2 <- glmmTMB(turnover ~ Origin * Deciduousness + Stage, family = ordbeta(link = "logit"), data = df_proportions)

check_model(m4.1)

# model selection (AIC) ----
# Raf
AIC(m1, m2, m3, m4)
# delta AIC = -0.1261

# Ma
AIC(m1.1, m2.1, m3.1, m4.1)
# delta AIC = -1.0706

# turnover
AIC(m1.2, m2.2, m3.2, m4.2)
# delta AIC = -0.5266

# population growth rate ----
growth_df <- growth

# models with a Gamma distribution
m1_lambda <- glmmTMB(growth_finite ~ Origin * Deciduousness,
                      family = Gamma(link = "log"),
                      data = growth_df)
m2_lambda <- glmmTMB(growth_finite ~ Origin + Deciduousness,
                      family = Gamma(link = "log"),
                      data = growth_df)
m3_lambda <- glmmTMB(growth_finite ~ Origin,
                      family = Gamma(link = "log"),
                      data = growth_df)
m4_lambda <- glmmTMB(growth_finite ~ Deciduousness,
                      family = Gamma(link = "log"),
                      data = growth_df)

# residual diagnostics ----
plot(simulateResiduals(m4_lambda))
check_model(m2_lambda)

summary(m3_lambda)
visreg::visreg(m3.1)

# AIC
AIC(m1_lambda, m2_lambda, m3_lambda, m4_lambda)

# table of fixed-effect estimates for the selected models ----
library(broom.mixed)
library(purrr)

format_publication_table <- function(model_list) {
  map_df(model_list, function(m) {
    tidy(m, effects = "fixed", component = "cond") %>%
      mutate(
        estimate = round(estimate, 4),
        std.error = round(std.error, 4),
        # keep p-value as text to avoid class conflicts when binding rows
        p_formatted = ifelse(p.value < 0.001, "< 0.001", as.character(round(p.value, 3))),
        formula = format(formula(m))
      ) %>%
      select(-p.value) %>%
      rename(p.value = p_formatted)
  }, .id = "Model")
}

final_table <- format_publication_table(list(lambda = m3_lambda, m3 = m3, m3.1, m3.2))

# save the table locally (originally exported to a private Google Sheet)
write.csv(final_table, "supplementary/model_estimates_table.csv", row.names = FALSE)
