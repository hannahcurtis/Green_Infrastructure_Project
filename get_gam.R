# Full GAM Process

get_gam <- function(df, num_rf_vars, p_threshold) {
  
# Random Forest

#make this example reproducible
set.seed(1)

# Define response variable
response <- df$Average.infiltration.rate

# Define predictor variables
# MANUALLY CHANGE HERE if interested in different predictors
predictors <- df[, c("Drainage.Area.Ratio", "Age", "Ponding.Depth", "Underdrain", "Percent.Vegetation.Shrub", "Percent.Vegetation.Grass", "Percent.Vegetation.Tree", 
                     "Percent.Vegetation.Prairie", "Percent.Vegetation.Other", "Percent.of.Surface.with.Poor.Infiltration", "Vegetation.Condition", "Engineered.Soil.Depth", 
                     "Engineered.Soil.Percent.Sand", "Engineered.Soil.Percent.Compost", "Soil_PC_1", "Soil_PC_2", "Percent.Nonresidential")]

rf_df <- data.frame(response, predictors)

#fit the random forest model
model <- randomForest(
  formula = response ~ .,
  data = rf_df
)

#display fitted model
model

# Predict the response on the training data
predictions <- predict(model, newdata = df_res_removed)

# Calculate R-squared
actual_values <- response
rss <- sum((actual_values - predictions)^2)  # Residual sum of squares (SSR)
tss <- sum((actual_values - mean(actual_values))^2)  # Total sum of squares (TSS)

r_squared <- 1 - (rss / tss)

# Print the R-squared value
print(r_squared)

# plot the test MSE by number of trees
plot(model)

# produce variable importance plot
varImpPlot(model, main = "")

# Find top 12 variables from model
# MANUALLY CHANGE HERE if interested in taking different number of predictors from random forest model
rf_predictors

# GAM
# TO-DO: 
  # figure out how to smooth non-categorical variables and leave categorical variables as linear
  # figure out what to do if GAM doesn't run due to low edf
  # figure out how to change vars from smoothed to linear based on edf and rerun model
  # figure out how to iteratively remove variables whose p-values are greater than 0.2

# run GAM
gam_model <- gam(response ~ rf_predictors, select = TRUE, data = df)

# Determine if the k values are sufficiently high for the smoothed variables
# If p-values are significant, k is too low
gam.check(gam_model)

#Output of model performance (R2, edf, p-values)
gam_summary <- summary(gam_model)
gam_summary

# Residuals plots, predicted vs measured plots, plots of each predictor vs response
plot(gam_model, all.terms = TRUE, residuals=TRUE, pch=1, cex=1, shade=TRUE, shade.col="lightblue")

# Get list of variables and associated p-values
# Extract p-values from the smooth terms and parametric terms
smooth_pvalues <- gam_summary$s.pv   # p-values for smoothed terms
param_pvalues <- gam_summary$p.pv    # p-values for parametric terms

# Combine smooth and parametric terms
pvalues <- c(smooth_pvalues, param_pvalues)

# Extract the names of the predictor variables in the same order
predictors <- c(rownames(gam_summary$s.table), rownames(gam_summary$p.table))

# Order the predictors by p-value
ordered_indices <- order(pvalues)
ordered_predictors <- predictors[ordered_indices]
ordered_pvalues <- pvalues[ordered_indices]

# Save the results in a list
results_list <- list(
  predictors = ordered_predictors,
  p_values = ordered_pvalues
)

# If edf is 1 or below, change variable from smoothed to parameter 

# Rerun GAM until all edf values are above 1

# Remove variable with the largest p-value and rerun
# Continue removing variables until all variables have p-value less than 0.2

# Create data frame with final predictors
final_predictors

# Check and plot VIF
#run lm on dataframe
vif_model <- lm(response ~ final_predictors, data = df)

# create vector of VIF values
vif_values <- vif(vif_model)

# Increase margin size so variable names all show up
par(mar=c(4,20,4,4))

# create horizontal bar chart to display each VIF value
barplot(vif_values, horiz = TRUE, cex.names = 1.5, cex.axis = 1.5, las=1, col = "steelblue", main="")

# add vertical lines at 5 and 10
abline(v = 5, lwd = 2, lty = 2)
abline(v = 10, lwd = 2, lty = 2)

# Check concurvity
concurvity(gam_model, full=TRUE)

}
