#########################################
# Benjamin Risk
# 7 October 2025
# Examine whether there are significant 
# differences between L and R tracts
#########################################

## ---- setup libs & data ----
library(mgcv)
library(gratia)     # derivatives()
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(stringr)
# library(itsadug)  ## not needed now, but fine to keep if you like

dat_harmonized <- read.csv('dat_harmonized_lr.csv', stringsAsFactors = TRUE) %>%
  filter(GA >= 36, GA <= 41)
dat_harmonized$Family_ID=factor(dat_harmonized$Family_ID)

par(mfrow=c(1,2))
boxplot(dat_harmonized$ILF_L_fa.combat)
boxplot(dat_harmonized$ILF_R_fa.combat)
diff = dat_harmonized$ILF_L_fa.combat-dat_harmonized$ILF_R_fa.combat
boxplot(diff)

dat_harmonized$diff=diff
library(lmerTest)
model.lmm <- lmer(diff ~ Mean_RMS + Sex + (1|Unique_ID) + (1|Family_ID),data = dat_harmonized)
summary(model.lmm)

model.diff.gam <- gam(diff ~ s(Corr_age) + s(Mean_RMS) + Sex +
    s(Unique_ID, bs = "re") + s(Family_ID, bs = "re"),
  method = "REML",
  data = dat_harmonized
)
plot(model.diff.gam)
summary(model.diff.gam)


## ---- function: make_long_tract_data ----
make_long_tract_data <- function(tract, metric, data = dat_harmonized) {
  # Build expected column names (L/R versions)
  left_col  <- paste0(tract, "_L_", metric, ".combat")
  right_col <- paste0(tract, "_R_", metric, ".combat")
  
  # Check which exist (some tracts might not have both hemispheres)
  cols_exist <- c(left_col, right_col)[c(left_col, right_col) %in% names(data)]
  if (length(cols_exist) == 0) {
    stop(paste0("No columns found for tract ", tract, " and metric ", metric))
  }
  
  # Select relevant variables
  base_vars <- c("Unique_ID", "Family_ID", "Mean_RMS", "Sex", "Corr_age", "chor_age")
  dat_sub <- data[, c(base_vars, cols_exist), drop = FALSE]
  
  # Reshape to long
  dat_long <- tidyr::pivot_longer(
    dat_sub,
    cols = all_of(cols_exist),
    names_to = "variable",
    values_to = "value"
  )
  
  # Extract hemisphere (L or R) from the variable name
  dat_long <- dat_long %>%
    dplyr::mutate(
      hemisphere = dplyr::case_when(
        stringr::str_detect(variable, "_L_") ~ "L",
        stringr::str_detect(variable, "_R_") ~ "R",
        TRUE ~ NA_character_
      ),
      tract = tract,
      metric = metric
    ) %>%
    dplyr::select(Unique_ID, Family_ID, Mean_RMS, Sex, Corr_age, chor_age,
                  tract, metric, hemisphere, value)
  
  dat_long$hemisphere = factor(dat_long$hemisphere)
  return(dat_long)
}

cor(dat_harmonized$UF_L_fa.combat,dat_harmonized$UF_R_fa.combat)
cor(dat_harmonized$AF_L_fa.combat,dat_harmonized$AF_R_fa.combat)
cor(dat_harmonized$ILF_L_fa.combat,dat_harmonized$ILF_R_fa.combat)
cor(dat_harmonized$IFOF_L_fa.combat,dat_harmonized$IFOF_R_fa.combat)
cor(dat_harmonized$Ci_L_fa.combat,dat_harmonized$Ci_R_fa.combat)
cor(dat_harmonized$Fx_L_fa.combat,dat_harmonized$Fx_R_fa.combat)
cor(dat_harmonized$PT_L_fa.combat,dat_harmonized$PT_R_fa.combat)
cor(dat_harmonized$ATR_L_fa.combat,dat_harmonized$ATR_R_fa.combat)


# correlations should be high in general:
cor(dat_harmonized$AF_L_fa.combat,dat_harmonized$ATR_R_fa.combat)
cor(dat_harmonized$Fx_L_fa.combat,dat_harmonized$ATR_R_fa.combat)

plot(dat_harmonized$AF_L_fa.combat,dat_harmonized$AF_R_fa.combat)

# ---- bilateral panel GAMMs per metric (2x4 each) ----
library(itsadug)

bilateral_tracts <- c("AF","ATR","Ci","Fx","IFOF","ILF","PT","UF")
metrics <- c("fa","md","ad","rd")

# optional: fixed y-lims per metric (comment out to let gg pick)
metric_ylim <- list(
  fa = c(0.10, 0.30),
  md = c(0.90, 1.50),
  ad = c(1.10, 1.80),
  rd = c(0.80, 1.30)
)

# hemisphere colors
# helper: get the left/right hex colors for a tract
table_colors = readxl::read_xlsx("~/Dropbox/WhiteMatterDevelopment/nifti/diffeo_warped/Table_Included_tracts reorderedCol.xlsx")

# overwrite the UF with pink color:
table_colors[1,7]=193
table_colors[2,7]=193


table_colors <- table_colors %>%
  dplyr::mutate(
    hex = grDevices::rgb(`Color RGB R`, `Color RGB G`, `Color RGB B`, maxColorValue = 255)
  )

get_hemi_colors <- function(tract) {
  left_name  <- paste0(tract, "_L")
  right_name <- paste0(tract, "_R")
  left_hex   <- table_colors %>% filter(Abbreviations == left_name)  %>% pull(hex)
  right_hex  <- table_colors %>% filter(Abbreviations == right_name) %>% pull(hex)
  
  c(L = left_hex, R = right_hex)
}

# Adjust AF_L and AF_R colors to more visible goldenrod shades
table_colors <- table_colors %>%
  mutate(
    `Color RGB R` = case_when(
      Abbreviations == "AF_L" ~ 218,
      Abbreviations == "AF_R" ~ 205,
      TRUE ~ `Color RGB R`
    ),
    `Color RGB G` = case_when(
      Abbreviations == "AF_L" ~ 165,
      Abbreviations == "AF_R" ~ 149,
      TRUE ~ `Color RGB G`
    ),
    `Color RGB B` = case_when(
      Abbreviations == "AF_L" ~ 32,
      Abbreviations == "AF_R" ~ 12,
      TRUE ~ `Color RGB B`
    )
  ) %>%
  mutate(hex = grDevices::rgb(`Color RGB R`, `Color RGB G`, `Color RGB B`, maxColorValue = 255))


fit_and_plot_metric <- function(metric) {
  #tiff(sprintf("../Figures/GAMM_%s_2x4_panels.tiff", toupper(metric)),    width = 1200, height = 700, pointsize = 14)
  tiff(sprintf("../Figures/GAMM_%s_2x4_panels_tableColors.tiff", toupper(metric)),
       width = 1200, height = 700, pointsize = 14)  # can also increase pointsize here
  
  par(mfrow = c(2,4), mar = c(5,5,3,2), cex.lab = 1.8, cex.axis = 1.4, cex.main = 1.8)
  
  
  
  
  for (tr in bilateral_tracts) {
    # long data for this tract+metric
    dt <- make_long_tract_data(tr, metric) |>
      dplyr::mutate(hemisphere = factor(hemisphere, levels = c("L","R")))
    
    # fit GAMM
    model.temp <- gam(
      value ~ s(Corr_age, by = hemisphere) + hemisphere + s(Mean_RMS) + Sex +
        s(Unique_ID, bs = "re") + s(Family_ID, bs = "re"),
      method = "REML",
      data = dt
    )
    
    # pick ylim if provided
    ylim_use <- metric_ylim[[metric]]
    # title = tract name (e.g., "ILF"), ylab = metric label
    ylab_use <- toupper(metric)
    
    # colors for table_colors
    hemi_cols = get_hemi_colors(tr)
    
    # Left hemisphere with simultaneous CI
    itsadug::plot_smooth(
      model.temp,
      view = "Corr_age",
      cond = list(hemisphere = "L"),
      print.summary = FALSE,
      hide.label = TRUE,
      ylab = ylab_use,
      xlab = "Corrected age (days)",
      sim.ci = TRUE,
      col = hemi_cols["L"],
      main = tr,
      ylim = ylim_use,
      lwd = 4,
      lty = 1,
      cex=2
    )
    # Right hemisphere overlay
    itsadug::plot_smooth(
      model.temp,
      view = "Corr_age",
      cond = list(hemisphere = "R"),
      print.summary = FALSE,
      hide.label = TRUE,
      sim.ci = TRUE,
      col = hemi_cols["R"],
      add = TRUE,
      lwd = 4,
      lty = 2,
      cex=2
    )
    
    # legend only on the first subplot
    if (tr == bilateral_tracts[1]) {
      legend("topright", legend = c("Left","Right"),
             lty = c(1,2), col = unname(hemi_cols), bty = "n", cex = 2, lwd = 4)
    }
  }
  
  dev.off()
}

# Run for each metric and create four files
invisible(lapply(metrics, fit_and_plot_metric))


## Small differences in IFOF, ILF, and PT
IFOF = make_long_tract_data('IFOF', 'fa') |>
   dplyr::mutate(hemisphere = factor(hemisphere, levels = c("L","R")))

model.temp <- gam(
   value ~ s(Corr_age, by = hemisphere) + hemisphere + s(Mean_RMS) + Sex +
     s(Unique_ID, bs = "re") + s(Family_ID, bs = "re"),
   method = "REML",
   data = IFOF
 )

plot_diff(model.temp,view='Corr_age',comp=list(hemisphere=c('L','R')),ylab='diff in FA',sim.ci=TRUE)
###############

ILF = make_long_tract_data('ILF', 'fa') |>
  dplyr::mutate(hemisphere = factor(hemisphere, levels = c("L","R")))

model.temp <- gam(
  value ~ s(Corr_age, by = hemisphere) + hemisphere + s(Mean_RMS) + Sex +
    s(Unique_ID, bs = "re") + s(Family_ID, bs = "re"),
  method = "REML",
  data = ILF
)

plot_diff(model.temp,view='Corr_age',comp=list(hemisphere=c('L','R')),ylab='diff in FA',sim.ci=TRUE)

###############
PT = make_long_tract_data('PT', 'fa') |>
  dplyr::mutate(hemisphere = factor(hemisphere, levels = c("L","R")))

model.temp <- gam(
  value ~ s(Corr_age, by = hemisphere) + hemisphere + s(Mean_RMS) + Sex +
    s(Unique_ID, bs = "re") + s(Family_ID, bs = "re"),
  method = "REML",
  data = PT
)

plot_diff(model.temp,view='Corr_age',comp=list(hemisphere=c('L','R')),ylab='diff in FA',sim.ci=TRUE)


