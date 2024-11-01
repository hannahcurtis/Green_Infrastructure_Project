# Full GAM Process

library(mgcv)
library(randomForest)
library(visreg)

get_gam <- function(df, x_vars, x_linear, y_var, pca_columns, num_rf_vars, p_threshold) {

# PCA 
  
# Run PCA on columns of interest
pca_results <- prcomp(pca_columns, scale=TRUE, center=TRUE)
  
# Summary of proportion of variance explained for each component 
print(summary(pca_results))
  
# Save first five principal components as variables
# Can save more components as needed
first_pc_soil_all <- pca_results$x[, 1]
second_pc_soil_all <- pca_results$x[, 2]
third_pc_soil_all <- pca_results$x[, 3]
fourth_pc_soil_all <- pca_results$x[, 4]
fifth_pc_soil_all <- pca_results$x[, 5]

# Create a dataframe with first five principal components
pca_results_df <- data.frame(Soil_PC_1 = first_pc_soil_all, Soil_PC_2 = second_pc_soil_all, 
                              Soil_PC_3 = third_pc_soil_all, Soil_PC_4 = fourth_pc_soil_all, Soil_PC_5 = fifth_pc_soil_all)
  
# Bind the first two principal components to the original data frame
df <- cbind(df, pca_results_df)

# Randfom Forest

#make this example reproducible
set.seed(1)

# Define response variable
response <- y_var

# Define smooth predictor variables
smooth_predictors <- df[, x_vars]

# Define linear predictor variables
linear_predictors <- df[, x_linear]

# Combine smooth and linear variables
predictors <- cbind(smooth_predictors, linear_predictors)

# create data frame of predictor and response variables of interest
rf_df <- data.frame(response, predictors)

#fit the random forest model
rf_model <- randomForest(
  formula = response ~ .,
  data = rf_df
)

#display fitted model
rf_model

# plot the test MSE by number of trees
plot(rf_model)

# produce variable importance plot
varImpPlot(rf_model, main = "")

# Determine most important x variables
importance_values <- importance(rf_model)  # Get importance values
top_vars <- head(order(importance_values[, "IncNodePurity"], decreasing = TRUE), num_rf_vars)  # Get indices of top x variables

# Get the names of the top x variables
rf_predictors <- rownames(importance_values)[top_vars]

# Convert variables to a list
rf_predictors <- as.list(rf_predictors)

# GAM

# Subsetted predictor variables to include in the model from random forest
subsetted_predictors <- rf_predictors

# Full list of predictor smoothed predictor variables 
all_predictors <- x_vars

# Parametric variables
initial_parametric_vars <- x_linear  # Your fixed parametric terms

# Initialize smooth and parametric terms from overlap of rf vars and all vars of interest
smooth_terms <- intersect(all_predictors, subsetted_predictors)
parametric_vars <- intersect(initial_parametric_vars, subsetted_predictors)

# Define the response variable
response_var <- "y_var" 

# Define the threshold for p-values and the tolerance for EDF comparison
p_value_threshold <- p_threshold
tolerance <- 1.000001  # Adjust tolerance if necessary

run_GAM <- function(smooth_terms, parametric_vars, response_var, data) {
  # Construct the formula with `s()` for smooth terms and linear for parametric terms
  formula <- as.formula(
    paste(response_var, "~", 
          paste(c(paste("s(", smooth_terms, ")", sep=""), parametric_vars), collapse = " + ")
    )
  )
  # Print model formula
  print(formula)
  
  return(gam(formula, select=TRUE, data = df))
}

# Loop until all terms have p-values below the threshold
iteration <- 1  # Initialize iteration counter
while (TRUE) {
  cat("Iteration:", iteration) # display iteration number
  model <- run_GAM(smooth_terms, parametric_vars, response_var, df) # run GAM
  
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
      
      # Remove 's()' to convert term from smooth to linear
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
    clean_term <- gsub("^s\\((.*)\\)$", "\\1", max_p_var)
    smooth_terms <- setdiff(smooth_terms, clean_term)
  } else {
    # If it's a parametric term and its p-value is high, remove it
    cat("Removing", max_p_var, "from parametric terms due to high p-value.\n")
    parametric_vars <- setdiff(parametric_vars, max_p_var)
  }
  
  # Increment iteration counter for diagnostics
  iteration <- iteration + 1
}

# After your loop
y_var <- df$Average.infiltration.rate

# view model check
print(gam.check(model))

# view final model summary
print(summary(model))

# Now use y_var in your analysis or plotting functions
#visreg(model, "parametric_variable_name", y_var = y_var)

# Save the final smooth and parametric terms as a character vector
final_terms <- c(smooth_terms, parametric_vars)

# Create a formula for the VIF model
formula_string <- paste("response ~", paste(final_terms, collapse = " + "))

# Run the lm model using the constructed formula
vif_model <- lm(as.formula(formula_string), data = df)

# Create vector of VIF values
vif_values <- vif(vif_model)

# widen plot so variable names aren't cut off
par(mar=c(4,25,4,4))

#create horizontal bar chart to display each VIF value
barplot(vif_values, horiz = TRUE, cex.names = 1.5, cex.axis = 1.5, las=1, col = "steelblue", main="")

#add vertical lines at 5 and 10 (considered high VIF score if above)
abline(v = 5, lwd = 2, lty = 2)
abline(v = 10, lwd = 2, lty = 2)

# Check concurvity
print(concurvity(model, full=TRUE))

}
