# We are looking at a 10 year interval and asking what causes a regional outbreak 
# (defined as 50% of the patches seeing severe defoliation). 
# Weather event that causes increasing levels of mortality. 
# Parasitoid community that has varying levels of dispersal diversity. 
# Mortality caused at a global scale so you can expect some patches hurt harder than others. 
# 10%, 20%, 30%, 40%, 50%, 60%, 70% mortality. 
# Dispersal diversity, so maybe 3 different levels of communities, 
# a low diversity, middle diversity, high diversity. This is a good start. 
# Something like “number of dead patches at end of 10 years.” 

# do 240 timesteps, assuming it might take time to precipitate an outbreak
# number of extinct patches after 100 simulations
low_0 <- pmax(round(rnorm(100, 3, 3)),0) 
low_20 <- pmax(round(rnorm(100, 3, 4)),0) 
low_40 <- pmin(pmax(round(rnorm(100, 10, 10)),0),22) 
low_60 <- pmax(round(rnorm(100, 15, 3)),0) 
low_80 <- pmin(round(rnorm(100, 23, 3)),22) 
mid_0 <- pmax(round(rnorm(100, 3, 2)),0) 
mid_20 <- pmax(round(rnorm(100, 3, 3)),0) 
mid_40 <- pmax(round(rnorm(100, 3, 6)),0) 
mid_60 <- pmax(round(rnorm(100, 14, 2)),0) 
mid_80 <- pmin(round(rnorm(100, 21, 3)),22) 
high_0 <- pmax(round(rnorm(100, 1, 2)),0) 
high_20 <- pmax(round(rnorm(100, 2, 3)),0) 
high_40 <- pmax(round(rnorm(100, 2, 4)),0) 
high_60 <- pmax(round(rnorm(100, 4, 6)),0) 
high_80 <- pmax(round(rnorm(100, 6, 6)),0) 

boxplot(low_0, low_20,low_40,low_60,low_80,mid_0, mid_20,mid_40,mid_60,mid_80,high_0, high_20,high_40,high_60,high_80)

library(ggplot2)
library(viridis)

# Combine all data into a data.frame with categories
data <- data.frame(
  value = c(low_0, low_20, low_40, low_60, low_80,
            mid_0, mid_20, mid_40, mid_60, mid_80,
            high_0, high_20, high_40, high_60, high_80),
  group = rep(c("low", "mid", "high"), each = 500), # 100 per variable, 5 variables per group
  category = rep(c("0", "20", "40", "60", "80"), times = 3, each = 100) # 5 categories per group
)

# Combine group and category into a single factor for x-axis ordering
data$group_category <- factor(paste(data$group, data$category, sep = "_"),
                              levels = c(
                                paste("low", c("0", "20", "40", "60", "80"), sep = "_"),
                                paste("mid", c("0", "20", "40", "60", "80"), sep = "_"),
                                paste("high", c("0", "20", "40", "60", "80"), sep = "_")
                              ))

# Create a cleaner x-axis grouping and ticks
data$tick_label <- rep(c("0%", "20%", "40%", "60%", "80%"), times = 3) # For each group
data$group <- factor(data$group, levels = c("low", "mid", "high")) # Ensure correct order

# Ensure the x-axis maintains the correct ordering
data$group_category <- factor(
  data$group_category,
  levels = c(
    paste("low", c("0", "20", "40", "60", "80"), sep = "_"),
    paste("mid", c("0", "20", "40", "60", "80"), sep = "_"),
    paste("high", c("0", "20", "40", "60", "80"), sep = "_")
  )
)

# Define a pastel color palette for the groups
pastel_colors <- c("low" = "#FDBF6F", "mid" = "#CAB2D6", "high" = "#FF9999")

ggplot(data, aes(x = group_category, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.5) + 
  scale_fill_manual(values = pastel_colors) + # Apply pastel color palette
  scale_x_discrete(
    labels = rep(c("0%", "20%", "40%", "60%", "80%"), 3), # Repeat percentage labels
    breaks = levels(data$group_category)                 # Ensure correct alignment
  ) +
  labs(
    x = "% mortality of FTC and parasitoid caused by disturbance", # Remove the x-axis title
    y = "# of dead forest patches after 240 time steps",
    title = "Number of dead forest patches following disturbances"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for readability
    legend.position = "none",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20) # Adjust margins for clarity
  ) +
  annotate("text", x = 3, y = max(data$value) * 1.05, label = "Low", size = 4, fontface = "bold") +
  annotate("text", x = 8, y = max(data$value) * 1.05, label = "Mid", size = 4, fontface = "bold") +
  annotate("text", x = 13, y = max(data$value) * 1.05, label = "High", size = 4, fontface = "bold")

# Load necessary libraries
library(dplyr)
library(gt)

# Calculate summary statistics
summary_table <- data %>%
  group_by(group) %>%
  summarise(
    Median = median(value),
    Mean = mean(value),
    Max = max(value),
    SD = sd(value)
  )

# Create a publication-ready table with gt
summary_table %>%
  gt() %>%
  tab_header(
    title = "Summary Statistics of Values by Group",
    subtitle = "Descriptive statistics for Low, Mid, and High groups"
  ) %>%
  fmt_number(
    columns = c(Median, Mean, Max, SD),
    decimals = 2
  ) %>%
  cols_label(
    group = "Group",
    Median = "Median",
    Mean = "Mean",
    Max = "Maximum",
    SD = "Standard Deviation"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.subtitle.font.size = 12,
    column_labels.font.weight = "bold",
    data_row.padding = px(5)
  )

# Load necessary libraries
library(dplyr)
library(gt)

# Calculate summary statistics by `% mortality`
mortality_summary <- data %>%
  group_by(category) %>%
  summarise(
    Median = median(value),
    Mean = mean(value),
    Max = max(value),
    SD = sd(value),
    .groups = "drop"
  )

# Create a publication-ready table with gt
mortality_summary %>%
  gt() %>%
  tab_header(
    title = "Summary Statistics of Values by % Mortality",
    subtitle = "Descriptive statistics for each mortality level"
  ) %>%
  fmt_number(
    columns = c(Median, Mean,Max, SD),
    decimals = 2
  ) %>%
  cols_label(
    category = "% Mortality",
    Median = "Median",
    Mean = "Mean",
    Max = "Maximum",
    SD = "Standard Deviation"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.subtitle.font.size = 12,
    column_labels.font.weight = "bold",
    data_row.padding = px(5)
  )


# Load necessary libraries
library(dplyr)
library(gt)

# Example: Perform ANOVA
# Assuming `data` contains the variables:
# `value` (dependent variable), `group` (factor: low, mid, high), and/or `category` (numeric: 0, 20, 40, etc.)
anova_model <- aov(value ~ group * category, data = data)

# Get the ANOVA table
anova_results <- summary(anova_model)

# Convert ANOVA results to a tidy data frame
anova_df <- as.data.frame(anova_results[[1]]) %>%
  tibble::rownames_to_column(var = "Term") %>%
  mutate(
    `F value` = round(`F value`, 2),
    `Pr(>F)` = signif(`Pr(>F)`, 3)  # Significant digits for p-values
  )

# Publication-quality table using `gt`
anova_df %>%
  gt() %>%
  tab_header(
    title = "ANOVA Results",
    subtitle = "Testing the effects of group, category, and their interaction on value"
  ) %>%
  fmt_number(
    columns = vars(`Sum Sq`, `Mean Sq`),
    decimals = 2
  ) %>%
  fmt_number(
    columns = vars(`F value`),
    decimals = 2
  ) %>%
  cols_label(
    Term = "Term",
    `Df` = "Degrees of Freedom",
    `Sum Sq` = "Sum of Squares",
    `Mean Sq` = "Mean Square",
    `F value` = "F-Statistic",
    `Pr(>F)` = "P-Value"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.subtitle.font.size = 12,
    column_labels.font.weight = "bold"
  )

# Perform Tukey HSD for pairwise comparisons
tukey_results <- TukeyHSD(anova_model)

# Print Tukey HSD results
print(tukey_results)

# Convert Tukey HSD results to a tidy data frame
tukey_df <- as.data.frame(tukey_results$`group:category`) %>%
  tibble::rownames_to_column(var = "Comparison") %>%
  mutate(
    `Adjusted P-Value` = signif(`p adj`, 3),  # Format p-values
    Significant = ifelse(`p adj` < 0.05, "Yes", "No")  # Mark significance
  )

# Create a publication-quality table using gt
library(gt)
tukey_df %>%
  gt() %>%
  tab_header(
    title = "Tukey HSD Post-Hoc Test Results",
    subtitle = "Significant Differences Between Group:Category Combinations"
  ) %>%
  fmt_number(
    columns = vars(diff, lwr, upr),
    decimals = 2
  ) %>%
  fmt_number(
    columns = vars(`Adjusted P-Value`),
    decimals = 3
  ) %>%
  cols_label(
    Comparison = "Comparison",
    diff = "Mean Difference",
    lwr = "Lower CI",
    upr = "Upper CI",
    `Adjusted P-Value` = "P-Value",
    Significant = "Significant?"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.subtitle.font.size = 12,
    column_labels.font.weight = "bold"
  )

library(emmeans)

# Compute estimated marginal means
em_means <- emmeans(anova_model, ~ group * category)

# Pairwise comparisons
pairs(em_means)

# Plot estimated marginal means
plot(em_means, comparisons = TRUE, type = "response")




# Reshape to wide format and back to long for ggplot
heatmap_long <- tukey_df %>%
  dplyr::select(Group1, Group2, P_Value) %>%
  pivot_wider(names_from = Group2, values_from = P_Value, values_fill = NA) %>%
  pivot_longer(cols = -Group1, names_to = "Group2", values_to = "P_Value")

# Clean up for ggplot: remove rows with NA (e.g., diagonal)
heatmap_long <- heatmap_long %>%
  filter(!is.na(P_Value))
# Plot the heatmap of p-values
ggplot(heatmap_long, aes(x = Group1, y = Group2, fill = P_Value)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C", trans = "log10", na.value = "gray90") +
  labs(
    title = "Heatmap of Tukey HSD Pairwise Differences",
    x = "Group:Category Combination",
    y = "Group:Category Combination",
    fill = "P-Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )
