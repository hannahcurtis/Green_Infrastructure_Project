# Full GAM Process

library(mgcv)
library(randomForest)

get_gam <- function(df, num_rf_vars, p_threshold) {
  
# Random Forest

#make this example reproducible
set.seed(1)

# Define response variable
response <- df_res_removed$Average.infiltration.rate
  
# Define predictor variables
# MANUALLY CHANGE HERE if interested in different predictors
predictors <- df_res_removed[, c("Drainage.Area.Ratio", "Age", "Ponding.Depth", "Underdrain", "Percent.Vegetation.Shrub", "Percent.Vegetation.Grass", "Percent.Vegetation.Tree", 
                                  "Percent.Vegetation.Prairie", "Percent.Vegetation.Other", "Percent.of.Surface.with.Poor.Infiltration", "Vegetation.Condition", "Engineered.Soil.Depth", 
                                  "Engineered.Soil.Percent.Sand", "Engineered.Soil.Percent.Compost", "Soil_PC_1", "Soil_PC_2", "Percent.Nonresidential")]
  
rf_df <- data.frame(response, predictors)
  
#fit the random forest model
rf_model <- randomForest(
  formula = response ~ .,
  data = rf_df
)
  
#display fitted model
rf_model
  
# Predict the response on the training data
predictions <- predict(rf_model, newdata = df_res_removed)
  
# Calculate R-squared
actual_values <- response
rss <- sum((actual_values - predictions)^2)  # Residual sum of squares (SSR)
tss <- sum((actual_values - mean(actual_values))^2)  # Total sum of squares (TSS)
  
r_squared <- 1 - (rss / tss)
  
# Print the R-squared value
print(r_squared)
  
# plot the test MSE by number of trees
plot(rf_model)
  
# produce variable importance plot
varImpPlot(rf_model, main = "")
  
# Assume you have already trained your random forest model called 'rf_model'
# Extract variable importance
importance_values <- importance(rf_model)
  
# Convert to a data frame for easier manipulation
importance_df <- as.data.frame(importance_values)
  
# Get the variable names and their importance scores
importance_df$Variable <- rownames(importance_df)
  
# Sort by importance (you may want to use a specific metric, e.g., Mean Decrease Accuracy)
importance_df <- importance_df[order(importance_df[, "IncNodePurity"], decreasing = TRUE), ]
  
# Select the top 12 variables
# MANUALLY CHANGE HERE if interested in taking different number of predictors from random forest model
top_vars <- head(importance_df, 12)
  
# Create the dataframe with only the variable names
rf_predictors <- data.frame(Variable = top_vars$Variable)


# GAM
# TO-DO: 
  # figure out what to do if GAM doesn't run due to low edf

# Full list of predictor variables (smooth terms only)
all_predictors <- c("Drainage.Area.Ratio", "Age", "Ponding.Depth", "Percent.Nonresidential", "Percent.Vegetation.Shrub", 
                    "Percent.Vegetation.Tree", "Percent.Vegetation.Grass", "Percent.Vegetation.Prairie", "Percent.Vegetation.Other", 
                    "Engineered.Soil.Depth", "Engineered.Soil.Percent.Sand", "Engineered.Soil.Percent.Compost", "Soil_PC_1", "Soil_PC_2")

# Subsetted predictor variables you want to include in the model
subsetted_predictors <- rf_predictors

# Parametric variables
initial_parametric_vars <- c("Underdrain", "Vegetation.Condition", "Percent.of.Surface.with.Poor.Infiltration")  # Your fixed parametric terms

# Initialize smooth and parametric terms
smooth_terms <- setdiff(intersect(all_predictors, subsetted_predictors), initial_parametric_vars)
parametric_vars <- initial_parametric_vars

# Set up the response variable and threshold for p-value
response_var <- "Average.infiltration.rate"  # Replace with your response variable

# Define the threshold for p-values and the tolerance for EDF comparison
p_value_threshold <- 0.2
tolerance <- 1.000001  # Adjust tolerance if necessary

# Loop until all terms have p-values below the threshold
iteration <- 1  # Initialize iteration counter
while (TRUE) {
  
  # Construct the formula with `s()` for smooth terms and linear for parametric terms
  formula <- as.formula(
    paste(response_var, "~", 
          paste(c(paste("s(", smooth_terms, ")", sep=""), parametric_vars), collapse = " + ")
    )
  )
  
  # Fit the GAM model
  model <- gam(formula, data = df_res_removed)
  
  # Print the current model formula for diagnostics
  cat("Iteration", iteration, "Model Formula:\n")
  print(formula)
  
  # Get p-values and EDFs for each term
  summary_model <- summary(model)
  smooth_p_values <- summary_model$s.table[, "p-value", drop = FALSE]  # p-values for smooth terms
  smooth_edfs <- summary_model$s.table[, "edf"]         # EDFs for smooth terms
  parametric_p_values <- summary_model$p.table[, "Pr(>|t|)", drop = FALSE]  # p-values for parametric terms
  
  # Display EDF values for debugging
  cat("EDF values for smooth terms:\n")
  print(smooth_edfs)
  
  # Combine p-values for both smooth and parametric terms
  all_p_values <- c(smooth_p_values, parametric_p_values)
  names(all_p_values) <- c(rownames(smooth_p_values), rownames(parametric_p_values))
  
  # Check if any p-values are above the threshold
  if (all(all_p_values < p_value_threshold, na.rm = TRUE)) {
    cat("All p-values are below the threshold. Ending loop.\n")
    break  # Exit loop if all terms have p-values below threshold
  }
  
  # Identify smooth terms with EDF ≤ 1
  low_edf_terms <- names(smooth_edfs[smooth_edfs <=  1.000001])
  
  if (length(low_edf_terms) > 0) {
    cat("Converting the following smooth terms to parametric due to low EDFs:\n", 
        paste(low_edf_terms, collapse = ", "), "\n")
    
    for (term in low_edf_terms) {
      # Print the term being processed for debugging
      cat("Processing term for conversion:", term, "\n")
      
      # Remove 's()' to get the clean term
      clean_term <- gsub("^s\\((.*)\\)$", "\\1", term)
      
      # Remove the smooth term directly using logical indexing
      smooth_terms <- smooth_terms[smooth_terms != clean_term]  # Exclude the smooth term
      
      # Add to parametric terms
      parametric_vars <- union(parametric_vars, clean_term)  
      
      # Print updated terms for debugging
      cat("Updated smooth terms:\n", paste(smooth_terms, collapse = ", "), "\n")
      cat("Updated parametric variables:\n", paste(parametric_vars, collapse = ", "), "\n")
    }
  } else {
    cat("No smooth terms with EDF ≤ 1 found for conversion.\n")
  }
  
  # Find the term with the highest p-value above the threshold
  max_p_var <- names(all_p_values)[which.max(all_p_values)]
  max_p_value <- max(all_p_values, na.rm = TRUE)  # Ensure no NA issues
  
  cat("Highest p-value term:", max_p_var, "with p-value:", max_p_value, "\n")
  
  # If the term with the highest p-value is a smooth term, remove it
  if (max_p_var %in% rownames(smooth_p_values)) {
    cat("Removing", max_p_var, "from smooth terms due to high p-value.\n")
    smooth_terms <- setdiff(smooth_terms, max_p_var)
  } else {
    # If it's a parametric term and its p-value is high, remove it
    cat("Removing", max_p_var, "from parametric terms due to high p-value.\n")
    parametric_vars <- setdiff(parametric_vars, max_p_var)
  }
  
  # Increment iteration counter for diagnostics
  iteration <- iteration + 1
}

# View final model summary
summary(model)

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
