# ==============================================================================
# 05 - Figures: predicted vs. observed vital rates and population growth
# ==============================================================================
# Input:  data-analysis/ontogeny.csv
#         data-analysis/growth_population.csv
# Output: figures/vital_rates_predictions.png
#         figures/population_growth_predictions.png
# ==============================================================================

library(glmmTMB)
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggeffects)
library(colorspace)
library(performance)
library(visreg)
library(scales)
library(ggsignif)

rm(list = ls()) # clean the environment
dir.create("figures", showWarnings = FALSE)

# data ----
ontogeny <- read.csv("data-analysis/ontogeny.csv")

df1 <- ontogeny %>%
  select(Species, Stage, Origin, Deciduousness, Raf, Ma, turnover) %>%
  mutate(Stage = as.factor(Stage),
         Origin = as.factor(Origin),
         Deciduousness = as.factor(Deciduousness)) %>%
  mutate(Raf = Raf / 100) %>%
  mutate(Ma = Ma / 100) %>%
  mutate(turnover = turnover / 100)

str(df1)

# 1 - models ----
# most parsimonious model
m_raf <- glmmTMB(Raf ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = df1)
m_ma <- glmmTMB(Ma ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = df1)
m_turnover <- glmmTMB(turnover ~ Origin + Stage + Deciduousness, family = ordbeta(link = "logit"), data = df1)

check_model(m_raf)
summary(m_turnover)
visreg::visreg(m_ma)
p <- visreg(m_ma, "Origin", plot = FALSE)


# 2 - predictions ----
# helper function to generate predictions for all 3 predictors of a model
get_all_preds <- function(model_obj, response_name) {

  # predictors to plot separately
  predictors <- c("Origin", "Stage", "Deciduousness")

  # generate predictions for each predictor
  all_preds <- lapply(predictors, function(predictor) {
    # 'terms' sets the predictor plotted on the x-axis
    # 'type = "fixed"' ignores random effects (none here, but good practice)
    # 'ci.lvl = 0.95' sets the 95% confidence interval
    pred_data <- ggpredict(model_obj, terms = predictor, type = "fixed", ci.lvl = 0.95) %>%
      as_tibble() %>%
      rename(predicted = predicted,
             conf.low = conf.low,
             conf.high = conf.high,
             x_value = x) %>%
      mutate(
        y_variable = response_name,
        predictor_plotted = predictor
      )
    return(pred_data)
  })

  # combine all 3 sets of predictions into a single data frame
  return(bind_rows(all_preds))
}

# predictions for each response variable
preds_raf <- get_all_preds(m_raf, "Raf")
preds_ma <- get_all_preds(m_ma, "Ma")
preds_turnover <- get_all_preds(m_turnover, "turnover")

# combine all predictions into one data frame
all_predictions <- bind_rows(preds_raf, preds_ma, preds_turnover)

# prepare the observed data ----
data_long <- df1 %>%
  pivot_longer(cols = c(Raf, Ma, turnover),
               names_to = "y_variable",
               values_to = "Observed_Value") %>%
  pivot_longer(cols = c(Origin, Deciduousness, Stage),
               names_to = "predictor_plotted",
               values_to = "x_value")

# build the x-axis column used for plotting
plot_data <- data_long %>%
  left_join(all_predictions, by = c("y_variable", "predictor_plotted", "x_value"))

plot_data <- plot_data %>%
  select(-std.error) %>%
  mutate(predictor_plotted = recode(predictor_plotted, "Deciduousness" = "Leaf habit"))

# rename variables for the figure labels
rename_labels <- as_labeller(c(
  Ma = "Mortality",
  Raf = "Recruitment",
  turnover = "Turnover"
))

# significance bar ----
significance_bar_df <- data.frame(
  predictor_plotted = "Stage",    # must match the predictor name exactly
  y_variable = unique(plot_data$y_variable), # applied to every row of that column
  start = 1,
  end = 2,
  y = 0.15,                       # bar height (adjust to the y-axis scale, in %)
  label = "***"
)

# also rename in the predictions data frame
all_predictions <- all_predictions %>%
  mutate(predictor_plotted = recode(predictor_plotted, "Deciduousness" = "Leaf habit"))

# figure ----
prediction_plot <- ggplot(plot_data, aes(x = x_value)) +

  # observed data points
  geom_point(aes(y = Observed_Value),
             fill = "#B3B5A1",
             color = "#B3B5A1",
             size = 2,
             alpha = 0.4) +
  # 95% confidence interval (error bars)
  geom_errorbar(aes(x = x_value, ymin = conf.low, ymax = conf.high),
                data = all_predictions,
                width = 0.2,
                color = "black",
                linewidth = 1) +

  # predicted mean
  geom_point(aes(x = x_value, y = predicted),
             data = all_predictions,
             color = "black",
             size = 3) +

  # significance bar layer
  geom_signif(
    data = significance_bar_df,
    aes(
      xmin = 1,
      xmax = 2,
      y_position = y,
      annotations = label
    ),
    manual = TRUE,
    inherit.aes = FALSE,   # avoids looking for aesthetics from the main plot
    tip_length = 0.03,
    vjust = 0.6,
    textsize = 5
  ) +

  # facets with the renamed labels
  facet_grid(y_variable ~ predictor_plotted,
             scales = "free_x",
             labeller = labeller(y_variable = rename_labels)) +

  labs(
    x = "",
    y = ""
  ) +
  scale_y_continuous(labels = percent_format(),
                      # extra space at the top of the plot
                      expand = expansion(mult = c(0.05, 0.15))
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    strip.text.y = element_text(face = "bold", size = 12),
    strip.text.x = element_text(face = "bold", size = 12),
    strip.placement = "outside",
    plot.margin = unit(c(1, 1, 1, 1), "cm"))


prediction_plot


ggsave(
  filename = "figures/vital_rates_predictions.png",
  plot = prediction_plot,
  width = 20,
  height = 14,
  units = "cm",
  dpi = 300
)

##### -------------- population growth ------------
rm(list = ls()) # clean the environment

# data ----
growth <- read.csv("data-analysis/growth_population.csv")

df <- growth %>%
  select(Species, Origin, Deciduousness, growth_finite) %>%
  mutate(Origin = as.factor(Origin),
         Deciduousness = as.factor(Deciduousness))

# model ----
g1 <- glmmTMB(growth_finite ~ Origin + Deciduousness, family = gaussian(link = "identity"), data = df)
check_model(g1)
simulateResiduals(g1, plot = TRUE) # no issues in residuals or dispersion
summary(g1)
visreg(g1)

# helper function to generate predictions for all predictors of a model
get_all_preds <- function(model_obj, response_name) {

  predictors <- c("Origin", "Deciduousness")

  all_preds <- lapply(predictors, function(predictor) {
    pred_data <- ggpredict(model_obj, terms = predictor, type = "fixed", ci.lvl = 0.95) %>%
      as_tibble() %>%
      rename(predicted = predicted,
             conf.low = conf.low,
             conf.high = conf.high,
             x_value = x) %>%
      mutate(
        y_variable = response_name,
        predictor_plotted = predictor
      )
    return(pred_data)
  })

  return(bind_rows(all_preds))
}

# predictions ----
preds_growth <- get_all_preds(g1, "growth_finite")

# observed data for the boxplot ----
data_long <- df %>%
  pivot_longer(
    cols = c(growth_finite),
    names_to = "y_variable",
    values_to = "Observed_Value"
  ) %>%
  # duplicate rows for each predictor so facet_wrap can display both panels
  expand_grid(predictor_plotted = c("Origin", "Deciduousness"))

# build the x-axis column for the boxplot/predictions
plot_data <- data_long %>%
  left_join(preds_growth, by = c("y_variable", "predictor_plotted")) %>%
  mutate(
    x_value_obs = case_when(
      predictor_plotted == "Origin" ~ Origin,
      predictor_plotted == "Deciduousness" ~ Deciduousness
    ))

rename_labels <- as_labeller(c(
  growth_finite = "Population Growth"
))

# figure ----
prediction_plot <- ggplot(plot_data, aes(x = x_value_obs)) +

  # observed data points
  geom_point(aes(y = Observed_Value),
             fill = "#B3B5A1",
             color = "#B3B5A1",
             size = 2,
             alpha = 0.4) +
  # 95% confidence interval (error bars)
  geom_errorbar(aes(x = x_value, ymin = conf.low, ymax = conf.high),
                data = preds_growth,
                width = 0.2,
                color = "black",
                linewidth = 1) +

  # predicted mean
  geom_point(aes(x = x_value, y = predicted),
             data = preds_growth,
             color = "black",
             size = 3) +

  facet_grid(y_variable ~ predictor_plotted,
             scales = "free_x") +
  theme(strip.text = element_blank()) +

  labs(
    x = "",
    y = "Population Growth"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 11),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    strip.text.y = element_blank(),
    strip.text.x = element_text(face = "bold", size = 12),
    strip.placement = "outside",
    strip.switch.y = "left",
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

prediction_plot


ggsave(
  filename = "figures/population_growth_predictions.png",
  plot = prediction_plot,
  width = 19,
  height = 13,
  units = "cm",
  dpi = 300
)
