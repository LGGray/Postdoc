# Power Calculation for CD8 and CD4 IL2RG Analysis using merged_data
# Using old data columns with .x suffix

library(pwr)

# ============================================================================
# Check available data and columns
# ============================================================================

cat("Objects in environment:\n")
print(ls())

cat("\n\nChecking for merged_data...\n")
if(exists("merged_data")) {
  cat("merged_data found!\n")
  cat("Dimensions:", dim(merged_data), "\n")
  cat("Column names:\n")
  print(colnames(merged_data))
  
  # Look for columns with .x suffix (old data)
  x_cols <- grep("\\.x$", colnames(merged_data), value=TRUE)
  cat("\nColumns with .x suffix (old data):\n")
  print(x_cols)
  
} else {
  cat("merged_data not found in environment\n")
}

# ============================================================================
# POWER ANALYSIS FUNCTIONS
# ============================================================================

power_analysis_wilcox <- function(data, outcome_var, group_var, var_name_display) {
  # Extract groups
  groups <- unique(na.omit(data[[group_var]]))
  
  group1 <- data[data[[group_var]] == groups[1], outcome_var]
  group2 <- data[data[[group_var]] == groups[2], outcome_var]
  
  # Remove NA values
  group1 <- na.omit(group1)
  group2 <- na.omit(group2)
  
  # Sample sizes
  n1 <- length(group1)
  n2 <- length(group2)
  
  # Skip if insufficient data
  if(n1 < 2 || n2 < 2) {
    cat("Insufficient data for", var_name_display, "\n")
    return(NULL)
  }
  
  # Descriptive statistics
  mean1 <- mean(group1)
  mean2 <- mean(group2)
  sd1 <- sd(group1)
  sd2 <- sd(group2)
  
  # Effect size (Cohen's d)
  cohens_d <- (mean1 - mean2) / sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2))
  
  # Perform Wilcox test for actual p-value
  wilcox_result <- wilcox.test(group1, group2)
  
  # Results
  results <- list(
    variable = var_name_display,
    column_name = outcome_var,
    n1 = n1,
    n2 = n2,
    total_n = n1 + n2,
    group1_name = as.character(groups[1]),
    group2_name = as.character(groups[2]),
    group1_mean = mean1,
    group2_mean = mean2,
    group1_sd = sd1,
    group2_sd = sd2,
    mean_diff = mean2 - mean1,
    cohens_d = cohens_d,
    wilcox_p_value = wilcox_result$p.value,
    wilcox_U_statistic = wilcox_result$statistic
  )
  
  return(results)
}

# ============================================================================
# POST-HOC POWER CALCULATION
# ============================================================================

post_hoc_power_ttest <- function(cohens_d, n1, n2, alpha = 0.05) {
  if(n1 == n2) {
    power_result <- pwr.t.test(n = n1, d = abs(cohens_d), sig.level = alpha, 
                               type = "two.sample", alternative = "two.sided")
  } else {
    n_avg <- (n1 + n2) / 2
    power_result <- pwr.t.test(n = n_avg, d = abs(cohens_d), sig.level = alpha,
                               type = "two.sample", alternative = "two.sided")
  }
  
  return(power_result$power)
}

# ============================================================================
# SAMPLE SIZE NEEDED FOR SPECIFIED POWER
# ============================================================================

sample_size_for_power <- function(cohens_d, target_power = 0.80, alpha = 0.05) {
  if(abs(cohens_d) < 0.01) {
    return(NA)  # Effect size too small
  }
  
  sample_size_result <- pwr.t.test(d = abs(cohens_d), power = target_power, 
                                    sig.level = alpha, type = "two.sample",
                                    alternative = "two.sided")
  
  return(ceiling(sample_size_result$n))
}

# ============================================================================
# PERFORM ANALYSIS ON MERGED DATA
# ============================================================================

if(exists("merged_data") && exists("CONDITION")) {
  
  # Find CD8 and CD4 columns with .x suffix
  cd8_col <- grep("CD8.*\\.x$", colnames(merged_data), value=TRUE, ignore.case=TRUE)
  cd4_col <- grep("CD4.*\\.x$", colnames(merged_data), value=TRUE, ignore.case=TRUE)
  
  cat("\n================================================\n")
  cat("POWER ANALYSIS - OLD DATA (with .x suffix)\n")
  cat("================================================\n")
  
  # CD8 Analysis
  if(length(cd8_col) > 0) {
    cat("\nSearching for CD8 column...\n")
    cd8_col <- cd8_col[1]
    cat("Using column:", cd8_col, "\n")
    
    cd8_power <- power_analysis_wilcox(merged_data, cd8_col, "CONDITION", "CD8 IL2RG (OLD DATA)")
    
    if(!is.null(cd8_power)) {
      cat("\n============================================\n")
      cat("CD8 IL2RG Analysis - OLD DATA\n")
      cat("============================================\n")
      cat("Sample sizes -", cd8_power$group1_name, ":", cd8_power$n1, "|", 
          cd8_power$group2_name, ":", cd8_power$n2, "\n")
      cat("Total N:", cd8_power$total_n, "\n")
      cat(cd8_power$group1_name, "Mean ± SD:", 
          round(cd8_power$group1_mean, 2), "±", round(cd8_power$group1_sd, 2), "\n")
      cat(cd8_power$group2_name, "Mean ± SD:", 
          round(cd8_power$group2_mean, 2), "±", round(cd8_power$group2_sd, 2), "\n")
      cat("Mean difference:", round(cd8_power$mean_diff, 2), "\n")
      cat("Cohen's d:", round(cd8_power$cohens_d, 3), "\n")
      cat("Wilcox p-value:", formatC(cd8_power$wilcox_p_value, format="e", digits=3), "\n")
      
      # Post-hoc power
      cd8_post_hoc_power <- post_hoc_power_ttest(cd8_power$cohens_d, cd8_power$n1, cd8_power$n2)
      cat("Post-hoc power (achieved with current N):", round(cd8_post_hoc_power, 3), "\n")
      
      # Sample size for 80% power
      cd8_n_needed <- sample_size_for_power(cd8_power$cohens_d, target_power = 0.80)
      if(!is.na(cd8_n_needed)) {
        cat("Sample size needed for 80% power:", cd8_n_needed, "per group\n")
        if(cd8_n_needed > cd8_power$n1) {
          cat("RECOMMENDATION: Current sample size may be UNDERPOWERED\n")
        } else if(cd8_n_needed < cd8_power$n1) {
          cat("RECOMMENDATION: Current sample size is ADEQUATE\n")
        }
      }
    }
  } else {
    cat("\nNo CD8 column with .x suffix found\n")
  }
  
  # CD4 Analysis
  if(length(cd4_col) > 0) {
    cat("\nSearching for CD4 column...\n")
    cd4_col <- cd4_col[1]
    cat("Using column:", cd4_col, "\n")
    
    cd4_power <- power_analysis_wilcox(merged_data, cd4_col, "CONDITION", "CD4 IL2RG (OLD DATA)")
    
    if(!is.null(cd4_power)) {
      cat("\n============================================\n")
      cat("CD4 IL2RG Analysis - OLD DATA\n")
      cat("============================================\n")
      cat("Sample sizes -", cd4_power$group1_name, ":", cd4_power$n1, "|", 
          cd4_power$group2_name, ":", cd4_power$n2, "\n")
      cat("Total N:", cd4_power$total_n, "\n")
      cat(cd4_power$group1_name, "Mean ± SD:", 
          round(cd4_power$group1_mean, 2), "±", round(cd4_power$group1_sd, 2), "\n")
      cat(cd4_power$group2_name, "Mean ± SD:", 
          round(cd4_power$group2_mean, 2), "±", round(cd4_power$group2_sd, 2), "\n")
      cat("Mean difference:", round(cd4_power$mean_diff, 2), "\n")
      cat("Cohen's d:", round(cd4_power$cohens_d, 3), "\n")
      cat("Wilcox p-value:", formatC(cd4_power$wilcox_p_value, format="e", digits=3), "\n")
      
      # Post-hoc power
      cd4_post_hoc_power <- post_hoc_power_ttest(cd4_power$cohens_d, cd4_power$n1, cd4_power$n2)
      cat("Post-hoc power (achieved with current N):", round(cd4_post_hoc_power, 3), "\n")
      
      # Sample size for 80% power
      cd4_n_needed <- sample_size_for_power(cd4_power$cohens_d, target_power = 0.80)
      if(!is.na(cd4_n_needed)) {
        cat("Sample size needed for 80% power:", cd4_n_needed, "per group\n")
        if(cd4_n_needed > cd4_power$n1) {
          cat("RECOMMENDATION: Current sample size may be UNDERPOWERED\n")
        } else if(cd4_n_needed < cd4_power$n1) {
          cat("RECOMMENDATION: Current sample size is ADEQUATE\n")
        }
      }
    }
  } else {
    cat("\nNo CD4 column with .x suffix found\n")
  }
  
  cat("\n================================================\n")
  cat("INTERPRETATION GUIDE\n")
  cat("================================================\n")
  cat("Cohen's d Effect Sizes:\n")
  cat("  0.2 = small\n")
  cat("  0.5 = medium\n")
  cat("  0.8 = large\n\n")
  cat("Power Interpretation:\n")
  cat("  < 0.80 = UNDERPOWERED (insufficient sample size)\n")
  cat("  0.80-0.90 = Adequate to good\n")
  cat("  > 0.90 = Conservative (larger than necessary)\n")
  
} else {
  cat("\nmerged_data or CONDITION not found in environment\n")
  cat("Please load the merged_data first\n")
}
