# Assignment 2: Factor Pricing

# Libraries
library(quantmod)
library(PerformanceAnalytics)
library(TTR)
library(zoo)
library(sandwich)
library(lmtest)
library(matrixStats)
library(plyr)


# Load provided data
load("D:/data_ex2_factormodels_20181129.Rdata")

# Constants
number_of_observations <- nrow(DT)
number_of_portfolios <- 100

# Functions
analize_and_fill_na <- function(){
  
  # Extracting the rows with missing data
  columns_with_NA <- DT[ , colSums(is.na(DT)) > 0]
  cat("\nnumber of columns with NA is: ", ncol(columns_with_NA))
  
  columns_names <- colnames(columns_with_NA)
  # Calculate the percentage of NA
  for(i in 1:ncol(columns_with_NA)){
    number_of_NA <- sum(is.na(columns_with_NA[i]))
    percentage <- number_of_NA/nrow(columns_with_NA)
    cat("\npercentage of NA in column ", columns_names[i]," is: ", percentage, '\n')
  }
  
  # Next we plot 4 graphs to get the structure 
  columns_with_NA[is.na(columns_with_NA)] <- -1
  for(i in 1:ncol(columns_with_NA)){
    plot(columns_with_NA[[i]], x = 1:length(columns_with_NA[[i]]))
  }
  
  DT[is.na(DT)] <<- 0
  cat("\nNA's been filled by 0\n")
}
get_portfolio_excess_ret <- function(){
  
  risk_free <- DT$RF
  portfolios <- cbind(DT[8:107])
  
  for (i in 1:ncol(portfolios)){
    portfolios[[i]] <- portfolios[[i]] - risk_free
  }
  
  return(portfolios)
}

get_factors <- function(){
  FF_Factors <- cbind(DT[2:6])
  FF_Factors <- cbind(FF_Factors, DT["HML_dev"])
  return(FF_Factors)
}

build_cross_sectional_model <- function(coeff_beta,
                                        portfolio_excess_ret,
                                        formula_gamma,
                                        header_gamma){
  
  coeffs_gamma <- matrix(nrow = number_of_observations, ncol = length(header_gamma))
  
  for (i in 1:number_of_observations) {
    formula_gamma <- as.formula(formula_gamma)
    cross_sec_model <- lm(formula_gamma)
    coeffs_gamma[i,] <- cross_sec_model$coefficients
  }
  
  colnames(coeffs_gamma) <- header_gamma
  coeffs_gamma <- as.data.frame(coeffs_gamma)
  
  print(summary(cross_sec_model))
  
  return(coeffs_gamma)
}

build_models <- function(portfolio_excess_ret, intercept, HML, HML_DEV ){
  
  FF_Factors <- get_factors()
  
  # No HML, No HML_DEV, with intercept
  formula_beta <- "portfolio_excess_ret[[j]] ~  FF_Factors$Mkt.RF +
                                                FF_Factors$SMB    +
                                                FF_Factors$RMW    +
                                                FF_Factors$CMA"
  header_beta <- c("beta.Mkt.rf",
                   "beta.SMB",
                   "beta.RMW",
                   "beta.CMA")
  
  # No HML, No HML_DEV, with intercept
  formula_gamma <- "t(portfolio_excess_ret[i,]) ~  coeff_beta$beta.Mkt.rf +
                                              coeff_beta$beta.SMB    +
                                              coeff_beta$beta.RMW    +
                                              coeff_beta$beta.CMA"
  header_gamma <- c("gamma.Mkt.rf",
                    "gamma.SMB",
                    "gamma.RMW",
                    "gamma.CMA")
  
  if(intercept == TRUE){
    header_beta <- c("alpha", header_beta)
    header_gamma <- c("alpha", header_gamma)
  }else{
    formula_beta <- paste(formula_beta,'-1')
    formula_gamma <- paste(formula_gamma,'-1')
  }
  
  if(HML == TRUE){
    formula_beta <- paste(formula_beta,'+',"FF_Factors$HML")
    header_beta <- c(header_beta, "beta.HML")
    formula_gamma <- paste(formula_gamma,'+',"coeff_beta$beta.HML")
    header_gamma <- c(header_gamma, "gamma.HML")
  }else{
    # HML not included
    if(HML_DEV == TRUE){
      formula_beta <- paste(formula_beta,'+',"FF_Factors$HML_dev")
      header_beta <- c(header_beta, "beta.HML_dev")
      
      formula_gamma <- paste(formula_gamma,'+',"coeff_beta$beta.HML_dev")
      header_gamma <- c(header_gamma, "gamma.HML_dev")
    }
  }
  
  coeff_beta <- matrix(nrow = number_of_portfolios, ncol = length(header_beta))
  
  colnames(coeff_beta) <- header_beta

  for (j in 1:number_of_portfolios){
    formula_beta <- as.formula(formula_beta)
    ff_5_model <- lm(formula_beta)
    coeff_beta[j,] <- ff_5_model$coefficients
    
  }

  
  coeff_beta <- as.data.frame(coeff_beta)
  
  coeffs_gamma <- build_cross_sectional_model(coeff_beta,
                                              portfolio_excess_ret,
                                              formula_gamma,
                                              header_gamma)
  
  coeffs <- list(coeff_beta = coeff_beta, coeffs_gamma = coeffs_gamma)
  
  return(coeffs)
}

calc_avg_risk_premia_for_factors <- function(coeff_cross_sectional_model){
  means <- colMeans(coeff_cross_sectional_model)
  means <- t(means)
  means <- data.frame(means)
  return(means)
}
calc_std_of_risk_premia_for_factors <- function(coeff_cross_sectional_model){
  names <- colnames(coeff_cross_sectional_model)
  stds <- colSds(as.matrix(coeff_cross_sectional_model))
  stds <- t(stds)
  stds <- data.frame(stds)
  colnames(stds) <- names

  return(stds)
}
find_95_conf_intervals <- function(means, stds, size){
  
  intervals <- data.frame(left = 0, right = 0)
  
  for(i in 1:ncol(means)){
    
    error <- qt(0.975, df = size - 1) * stds[[i]] / sqrt(size)
    left <- means[[i]] - error
    right <- means[[i]] + error
    
    intervals[i,] <- c(left, right)
  }
  
  names <- colnames(means)
  rownames(intervals) <- names
  
  return(intervals)

}
calculate_stats <- function(param_list, number_of_observations){
  
  element_names <- names(param_list)
  
  means <- calc_avg_risk_premia_for_factors(param_list[[1]]$coeffs_gamma)
  stds <- calc_std_of_risk_premia_for_factors(param_list[[1]]$coeffs_gamma)
  intervals <- find_95_conf_intervals(means, stds, number_of_observations)
  intervals <- round(intervals, 4)
  intervals_list <- list(intervals)
  
  for(i in 2:length(param_list)){
    means_i <- calc_avg_risk_premia_for_factors(param_list[[i]]$coeffs_gamma)
    stds_i <- calc_std_of_risk_premia_for_factors(param_list[[i]]$coeffs_gamma)
    intervals_i <- find_95_conf_intervals(means_i, stds_i, number_of_observations)
    
    means <- rbind.fill(means, means_i)
    stds <- rbind.fill(stds, stds_i)
    intervals <- round(intervals_i, 4)
    intervals_list[[i]] <- intervals_i
  }
  
  means <- round(means, 4)
  stds <- round(stds, 4)
  
  
  rownames(means) <- element_names
  rownames(stds) <- element_names
  names(intervals_list) <- element_names
  
  return_list <- list(means = means, stds = stds, intervals = intervals_list)
  return(return_list)
  
}

value_redundancy_check <- function(){
  # check if the value redundant
  FF_Factors <- get_factors()
  regrHML <- lm(FF_Factors$HML ~ FF_Factors$Mkt.RF + FF_Factors$SMB + FF_Factors$RMW + FF_Factors$CMA)
  return(summary(regrHML))
  
}

# MAIN
analize_and_fill_na()
portfolio_excess_ret <- get_portfolio_excess_ret()
coeffs_ttf <- build_models(portfolio_excess_ret, intercept = TRUE, HML = TRUE, HML_DEV = FALSE)
coeffs_tff <- build_models(portfolio_excess_ret, intercept = TRUE, HML = FALSE, HML_DEV = FALSE)
coeffs_tft <- build_models(portfolio_excess_ret, intercept = TRUE, HML = FALSE, HML_DEV = TRUE)
coeffs_ftf <- build_models(portfolio_excess_ret, intercept = FALSE, HML = TRUE, HML_DEV = FALSE)
coeffs_fff <- build_models(portfolio_excess_ret, intercept = FALSE, HML = FALSE, HML_DEV = FALSE)
coeffs_fft <- build_models(portfolio_excess_ret, intercept = FALSE, HML = FALSE, HML_DEV = TRUE)

param_list <- list(ttf = coeffs_ttf,
                   tff = coeffs_tff,
                   tft = coeffs_tft,
                   ftf = coeffs_ftf,
                   fff = coeffs_fff,
                   fft = coeffs_fft)

return_list <- calculate_stats(param_list, number_of_observations)


value_redundancy_check()










