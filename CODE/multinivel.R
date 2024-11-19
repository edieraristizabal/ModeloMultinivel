#importar librerias
library(sf) #para importar datos geoespaciales
library(lme4) # para el modelo
library(pscl) # para calcular los R2
library(MuMIn) #para calcular los R2 en el modelo multinivel 
library(ggplot2) # Plotting library
library(ggspatial) # For adding north arrow and scale to maps
library(broom.mixed) # para transformar tabla
library(kableExtra) #para crear tabla en LaTEX
library(sjPlot) #para graficar los effectos
library(knitr) #para LaTEX table
library(broom)
library(dplyr) # for data manipulation
library(grid) #to export the table
library(pROC) # for the ROC curve

#Cargar datos
data = st_read("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_Multiniveles/DATA/df_catchments_spatial.gpkg", quiet = TRUE)

#Para guardar datos
savefile <- "G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_Multiniveles/FIGURAS/"

# Create Y_bin as 1 if lands_rec is 1 or greater, otherwise 0
data$Y_bin <- ifelse(data$lands_rec >= 1, 1, 0)

# Convert 'cuenca' to a factor (categorical variable)
data$cuenca <- as.factor(data$cuenca)
data$kmeans <- as.factor(data$kmeans)

# Ensure predictor variables are standardized
data$elev_mean_std <- scale(data$elev_mean)
data$rel_mean_std <- scale(data$rel_mean)
data$area_std <- scale(data$area)
data$rainfallAnnual_mean_std <- scale(data$rainfallAnnual_mean)

#########################################

# Fit the logistic regression model without random effects
m1 <- glm(
  Y_bin ~ elev_mean_std + rel_mean_std + area_std + rainfallAnnual_mean_std,
  data = data, 
  family = binomial(link = "logit")
)

# Calculate pseudo-R^2 values
pseudo_r2 <- pR2(m1)

# Fitted values for model m3
data$m1fitted <- fitted(m1) 

# Extract model coefficients as a tidy data frame
model_summary <- tidy(m1) %>%
  select(term, estimate, std.error, statistic, p.value) %>%
  rename(
    Variables = term,
    Estimate = estimate,
    Std_Error = std.error,
    Z_value = statistic,
    P_value = p.value
  )

# Rename the variables
model_summary$Variables <- recode(
  model_summary$Variables, 
  "(Intercept)" = "(Intercept)", 
  "elev_mean_std" = "Elevación", 
  "rel_mean_std" = "Relieve", 
  "area_std" = "Área", 
  "rainfallAnnual_mean_std" = "Lluvia anual"
)

# Convert p-values to scientific notation for clarity and reduce precision of all values to 3 decimal places
model_summary <- model_summary %>%
  mutate(
    Estimate = round(Estimate, 3),
    Std_Error = round(Std_Error, 3),
    Z_value = round(Z_value, 3),
    P_value = format.pval(P_value, digits = 3, eps = 0.001)
  )

# Create a separate data frame for pseudo-R^2 values and AIC with the same column structure
pseudo_r2_df <- data.frame(
  Variables = c("McFadden R2", "Cox & Snell R2", "Nagelkerke R2", "AIC"),
  Estimate = round(c(as.numeric(unlist(pseudo_r2[1:3])), AIC(m1)), 3),
  Std_Error = NA,
  Z_value = NA,
  P_value = NA
)

# Combine the model summary with the pseudo-R2 values
final_table <- bind_rows(model_summary, pseudo_r2_df)

# Create a LaTeX table using kable
kable(final_table, format = "latex", booktabs = TRUE, 
      caption = "Logistic Regression Model Summary with Pseudo-R^2 and AIC", 
      align = "c", escape = TRUE) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

##############################################

# Fit the multilevel logistic regression model with random slopes for each predictor
m2 <- glmer(
  Y_bin ~ elev_mean_std + rel_mean_std + area_std + rainfallAnnual_mean_std + 
    (1 + elev_mean_std + rel_mean_std + area_std + rainfallAnnual_mean_std | cuenca),
  data = data,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# Print the summary of the model to see the results
summary(m2)

# Fitted values for model m3
data$m2fitted <- fitted(m2) 

# Extract random effects from the model
ranef_data <- ranef(m2)$cuenca

# Convert to data frame for display
ranef_df <- as.data.frame(ranef_data)
colnames(ranef_df) <- c("Intercept", "Elevación Media Estándar", "Relieve Medio Estándar", "Área Estándar", "Lluvia Anual Media Estándar")
rownames(ranef_df) <- rownames(ranef_data) # Set rownames

# Add a column for cuenca names
ranef_df <- cbind(Cuenca = rownames(ranef_df), ranef_df)
rownames(ranef_df) <- NULL  # Remove rownames from data frame

# Print the data frame as a LaTeX table
kable(
  ranef_df,
  format = "latex",
  booktabs = TRUE,
  caption = "Random Effects for Cuenca",
  align = "c"
)

# Capture the model summary output
summary_text <- capture.output(summary(m2))

# Create a PNG file with higher resolution
png("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_Multiniveles/FIGURAS/m2_summary.png", width = 4000, height = 5300, res = 500)

# Plot the summary text as an image using grid graphics
grid.newpage()
text_grob <- textGrob(
  paste(summary_text, collapse = "\n"), 
  x = 0, y = 1, just = c("left", "top"), 
  gp = gpar(fontsize = 13, fontfamily = "sans")  # Adjusted fontsize for clarity
)
grid.draw(text_grob)

# Close the PNG device
dev.off()

####################################################

# Fit the multilevel logistic regression model with random slopes for each predictor
m3 <- glmer(
  Y_bin ~ elev_mean_std + rel_mean_std + area_std + rainfallAnnual_mean_std + 
    (1 + elev_mean_std + rel_mean_std + area_std + rainfallAnnual_mean_std | kmeans),
  data = data,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# Print the summary of the model to see the results
summary(m3)

# Fitted values for model m3
data$m3fitted <- fitted(m3)  

# Extract random effects from the model
ranef_data <- ranef(m3)$kmeans

# Convert to data frame for display
ranef_df <- as.data.frame(ranef_data)
colnames(ranef_df) <- c("Intercept", "Elevación Media Estándar", "Relieve Medio Estándar", "Área Estándar", "Lluvia Anual Media Estándar")
rownames(ranef_df) <- rownames(ranef_data) # Set rownames

# Add a column for cuenca names
ranef_df <- cbind(KNN5 = rownames(ranef_df), ranef_df)
rownames(ranef_df) <- NULL  # Remove rownames from data frame

# Print the data frame as a LaTeX table
kable(
  ranef_df,
  format = "latex",
  booktabs = TRUE,
  caption = "Random Effects for Cuenca",
  align = "c"
)

###################################################

# Create summary data frames for m2 and m3 with both fixed and random effects
model_summaries <- data.frame(
  Estadístico = c(
    "AIC", "BIC", "Log Verosimilitud", "Deviance", 
    "Efectos fijos: Intercepto", "Efectos fijos: Elevación", 
    "Efectos fijos: Relieve", "Efectos fijos: Área", 
    "Efectos fijos: Lluvia Anual",
    "Efectos aleatorios (Varianza): Intercepto", "Efectos aleatorios (Varianza): Elevación", 
    "Efectos aleatorios (Varianza): Relieve", 
    "Efectos aleatorios (Varianza): Área", 
    "Efectos aleatorios (Varianza): Lluvia Anual"
  ),
  Modelo_m2 = c(
    405.5, 490.8, -182.7, 365.5, 
    2.0011, 0.8566, 1.0958, 0.7033, -0.3884,
    0.000000, 0.003281, 0.829756, 0.007460, 0.278068
  ),
  Modelo_m3 = c(
    420.5, 505.8, -190.2, 380.5, 
    2.0580, 0.7834, 1.2863, 1.2102, -0.4966,
    0.0000000, 0.0006873, 0.1015634, 0.3551130, 0.0007044
  )
)

# Round values for better display
model_summaries$Modelo_m2 <- round(model_summaries$Modelo_m2, 4)
model_summaries$Modelo_m3 <- round(model_summaries$Modelo_m3, 4)

# Create a LaTeX table
kable(
  model_summaries,
  format = "latex",
  booktabs = TRUE,
  caption = "Comparación de Modelos de Regresión Logística Multinivel (m2 y m3) con Efectos Fijos y Aleatorios",
  col.names = c("Estadístico", "Modelo m2", "Modelo m3")
) %>%
  kable_styling(latex_options = c("striped", "hold_position"))
###################################################

# Predicted probabilities for the models
predicted_probs_m1 <- predict(m1, type = "response")
predicted_probs_m2 <- predict(m2, type = "response")
predicted_probs_m3 <- predict(m3, type = "response")

# Generate the ROC curve and calculate AUC for model m2
roc_m1 <- roc(data$Y_bin, predicted_probs_m1)
auc_m1 <- auc(roc_m1)

# Generate the ROC curve and calculate AUC for model m2
roc_m2 <- roc(data$Y_bin, predicted_probs_m2)
auc_m2 <- auc(roc_m2)

# Generate the ROC curve and calculate AUC for model m3
roc_m3 <- roc(data$Y_bin, predicted_probs_m3)
auc_m3 <- auc(roc_m3)

# Create data frames for plotting
roc_data_m1 <- data.frame(
  specificity = 1 - roc_m1$specificities,
  sensitivity = roc_m1$sensitivities,
  model = paste("Modelo M1 (AUC =", round(auc_m1, 3), ")")
)

# Create data frames for plotting
roc_data_m2 <- data.frame(
  specificity = 1 - roc_m2$specificities,
  sensitivity = roc_m2$sensitivities,
  model = paste("Modelo cuencas (AUC =", round(auc_m2, 3), ")")
)

roc_data_m3 <- data.frame(
  specificity = 1 - roc_m3$specificities,
  sensitivity = roc_m3$sensitivities,
  model = paste("Modelo clusters (AUC =", round(auc_m3, 3), ")")
)

roc_data <- rbind(roc_data_m1, roc_data_m2, roc_data_m3)

# Plotting the ROC curves with AUC annotations
roc_plot <- ggplot(roc_data, aes(x = specificity, y = sensitivity, color = model)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  theme_minimal() +
  theme(
    legend.position = c(0.8, 0.2),  # Position the legend inside the plot
    panel.border = element_rect(color = "black", fill = NA, size = 1), # Add border
    legend.title = element_blank()  # Remove legend title
  ) +
  labs(
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  )

# Display the plot
print(roc_plot)

# Save the plot to a file
ggsave("G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_Multiniveles/FIGURAS/roc_curve_plot.png", plot = roc_plot, width = 8, height = 6, dpi = 300)

##################################################

g <- ggplot() +
  geom_sf(data = data, aes(fill = m2_fitted), color = "black") +
  annotation_scale(location = "br", style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.7, "cm"), width = unit(0.6, "cm")) +
  scale_fill_gradientn(colors = c("blueviolet","blue","cornflowerblue","cyan2","green","greenyellow", "yellow", "orange","orangered","red"), name ="Susceptibilidad") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.text.x = element_text(size = 8, margin = margin(t = 1, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size = 8, margin = margin(t = 0, r = 1, b = 0, l = 0)),
    legend.text = element_text(size = 6),
    legend.title.align = 0,
    legend.position = c(0.3, 0.9),
    legend.key.size = unit(0.5, 'cm'),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.title = element_text(size = 10, vjust = .8, hjust = .5)
  )

# Plot the figure before saving
print(g)

ggsave(filename = "G:/My Drive/INVESTIGACION/PAPERS/ELABORACION/Modelo_Multiniveles/FIGURAS/Mapafinal.png", plot = g, width = 7, height = 8.5, units = "in", dpi = 500)
