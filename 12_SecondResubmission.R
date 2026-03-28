library(mgcv)
library(itsadug)
library(dplyr)
library(gratia)
library(ggplot2)
library(gridExtra)
library(kableExtra)

dat_harmonized = read.csv('dat_harmonized_lr.csv',stringsAsFactors = TRUE)
dat_harmonized$Family_ID=factor(dat_harmonized$Family_ID)
dat_harmonized = dat_harmonized%>%filter(GA>=36 & GA<=41)

# Check how many scans we have at 36 and 41 weeks of age:
table(dat_harmonized$GA)
dat_harmonized%>%count(GA)
# 9 scans from children born at 36 weeks and 9 scans from children born at 41 weeks

# how many children:
unique_children = dat_harmonized[!duplicated(dat_harmonized$Unique_ID),c('Unique_ID','GA')]

# subsetting to 36 to 41:
unique_children%>%count(GA)

#########################
# Show that we need a GAMM:
library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)
library(gridExtra)

metrics <- c("fa", "md", "ad", "rd")
comparison_results <- data.frame()
plot_list <- list()

for (m in metrics) {
  target <- paste0("WB_", m, ".combat")
  
  # 1. Fit Models
  mod_gamm <- gam(as.formula(paste0(target, " ~ s(Mean_RMS) + s(Corr_age) + s(Unique_ID, bs='re') + s(Family_ID, bs='re')")), 
                  method = 'ML', data = dat_harmonized)
  
  mod_quad <- gam(as.formula(paste0(target, " ~ s(Mean_RMS) + poly(Corr_age, 2, raw=TRUE) + s(Unique_ID, bs='re') + s(Family_ID, bs='re')")), 
                  method = 'ML', data = dat_harmonized)
  
  # 2. Comparison Table
  comparison_results <- rbind(comparison_results, data.frame(
    Metric = toupper(m),
    GAMM_AIC = AIC(mod_gamm),
    GAMM_EDF = summary(mod_gamm)$s.table["s(Corr_age)", "edf"],
    Quad_AIC = AIC(mod_quad),
    Quad_DF = 2,
    Delta_AIC = AIC(mod_quad) - AIC(mod_gamm)
  ))
  
  # 3. Extract GAMM Derivative (Handling the dotted column names)
  der_gamm <- derivatives(mod_gamm, select = "s(Corr_age)", n = 200) %>%
    mutate(
      derivative = .derivative,
      lower = .lower_ci,
      upper = .upper_ci,
      Model = "GAMM (Flexible)"
    ) %>%
    rename(age = Corr_age) # gratia uses the variable name from the model
  
  # 4. Calculate Quadratic Derivative Manually
  beta <- coef(mod_quad)
  V <- vcov(mod_quad)
  name1 <- "poly(Corr_age, 2, raw = TRUE)1"
  name2 <- "poly(Corr_age, 2, raw = TRUE)2"
  
  age_seq <- seq(min(dat_harmonized$Corr_age), max(dat_harmonized$Corr_age), length.out = 200)
  
  der_quad <- data.frame(age = age_seq) %>%
    mutate(
      derivative = beta[name1] + 2 * beta[name2] * age,
      se = sqrt(V[name1, name1] + (2 * age)^2 * V[name2, name2] + 2 * (2 * age) * V[name1, name2]),
      lower = derivative - 1.96 * se,
      upper = derivative + 1.96 * se,
      Model = "Quadratic (Linear Deriv)"
    )
  
  # 5. Build the Plot
  # Combine only the necessary columns for plotting
  p_data <- bind_rows(
    der_gamm[, c("age", "derivative", "lower", "upper", "Model")],
    der_quad[, c("age", "derivative", "lower", "upper", "Model")]
  )
  
  plot_list[[m]] <- ggplot(p_data, aes(x = age, y = derivative, group = Model)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model), alpha = 0.2) +
    geom_line(aes(color = Model), size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = toupper(m), x = "Corrected Age", y = "Growth Rate (Deriv)") +
    theme_minimal() + 
    theme(legend.position = "none")
}

# --- Output ---
# Table
library(kableExtra)

comparison_results %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(Delta_AIC = cell_spec(Delta_AIC, color = "white", bold = T,
                               background = spec_color(Delta_AIC, end = 0.9, option = "D", direction = -1))) %>%
  kbl(escape = F, align = "c", caption = "Model Selection: GAMM (Splines) vs. Parametric Quadratic") %>%
  kable_paper("hover", full_width = F) %>%
  add_header_above(c(" " = 1, "GAMM (Flexible)" = 2, "Quadratic (Rigid)" = 2, "Comparison" = 1)) %>%
  column_spec(1, bold = T, border_right = T) %>%
  column_spec(6, bold = T)


# Grid (Shared legend logic)
# Note: uses gridExtra and assumes common color/fill mappings
pdf(file='../Figures/SecondRevision_RtoR_QuadraticGrowthRates.pdf')
grid.arrange(grobs = plot_list, ncol = 2, top = "Comparison of Developmental Growth Rates")
dev.off()

###################
library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)
library(patchwork) # Excellent for combined plots and shared legends

metrics <- c("fa", "md", "ad", "rd")
plot_list <- list()

# 1. Prepare the restricted dataset
dat_restricted <- dat_harmonized %>% filter(GA > 36 & GA < 41)

for (m in metrics) {
  target <- paste0("WB_", m, ".combat")
  form <- as.formula(paste0(target, " ~ s(Mean_RMS) + s(Corr_age) + s(Unique_ID, bs='re') + s(Family_ID, bs='re')"))
  
  # 2. Fit both models
  mod_full <- gam(form, method = 'REML', data = dat_harmonized)
  mod_rest <- gam(form, method = 'REML', data = dat_restricted)
  
  # 3. Create prediction data (Corr_age 0 to 200)
  p_dat <- data.frame(
    Corr_age = seq(0, 200, length.out = 300),
    Mean_RMS = mean(dat_harmonized$Mean_RMS),
    Unique_ID = dat_harmonized$Unique_ID[1],
    Family_ID = dat_harmonized$Family_ID[1]
  )
  
  # 4. Get Simultaneous Confidence Intervals
  pred_full <- fitted_values(mod_full, data = p_dat, exclude = c("s(Unique_ID)", "s(Family_ID)"), 
                             intervals = "simultaneous") %>%
    mutate(Model = "Original (GA 36-41)")
  
  pred_rest <- fitted_values(mod_rest, data = p_dat, exclude = c("s(Unique_ID)", "s(Family_ID)"), 
                             intervals = "simultaneous") %>%
    mutate(Model = "Restricted (GA 37-40)")
  
  plot_data <- bind_rows(pred_full, pred_rest)
  
  # 5. Build Individual Plot
  plot_list[[m]] <- ggplot(plot_data, aes(x = Corr_age, y = .fitted, group = Model)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = Model), alpha = 0.2) +
    geom_line(aes(color = Model), linewidth = 1) +
    labs(title = toupper(m), 
         x = "Corrected Age", 
         y = "Predicted Value") +
    scale_color_manual(values = c("Original (GA 36-41)" = "#2C3E50", "Restricted (GA 37-40)" = "#E67E22")) +
    scale_fill_manual(values = c("Original (GA 36-41)" = "#2C3E50", "Restricted (GA 37-40)" = "#E67E22")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# 6. Combine plots with a shared legend using patchwork
# This automatically handles the alignment and legend extraction
combined_plot <- (plot_list$fa + plot_list$md) / (plot_list$ad + plot_list$rd) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Render the plot
combined_plot + plot_annotation(
  title = "Sensitivity Analysis: Impact of GA Inclusion Range",
  subtitle = "Comparing trajectories when excluding Gestational Ages 36 & 41",
  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5))
)

