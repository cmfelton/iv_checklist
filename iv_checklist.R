#Below, we provide example R code for generating the sensitivity plots shown in the main text. 
#We name the dataset sibs. The outcome (private school attendance) is named priv. 
#The treatment (having two siblings) is named treat. The instrument (having the first two children be of the same sex) is 
#named samesex. All other variables are additional covariates included in their model specifications.

##CLEARING ENVIRONMENT
rm(list = ls())

##INSTALLING AND LOADING NECESSARY PACKAGES

#install.packages("tidyverse")
#install.packages("lmtest")
#install.packages("sandwich")
#install.packages("AER")
#install.packages("ivpack")
#install.packages("boot")
#install.packages("assertive")
#install.packages("rlang")
#install.packages("sensemakr")
#install.packages("retrodesign")

library(tidyverse)
library(lmtest)
library(sandwich)
library(AER)
library(ivpack)
library(boot)
library(assertive)
library(rlang)
library(sensemakr)
library(retrodesign)

##SET WORKING DIRECTORY
#here you can use setwd() or set working 
#directory to source file location and have
#conley.rds in the same directory as this script

##READING IN DATA

sibs <- readRDS("cg2006.rds")

##for .csv file:
#df_name <- read_csv("filename")

##for .dta file:
#install.packages(readstata13)
#library(readstata13)
#df_name <- read.dta13("filename")

##ROBUST F-TEST
#requires lmtest package

#first-stage model WITHOUT instrument
mod_no_iv <- lm(treat ~ parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                  age_diff + factor(AGE_kid2), data = sibs)

#first-stage model WITH instrument
mod_iv <- lm(treat ~ samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
               age_diff + factor(AGE_kid2), data = sibs)

#partial F-test using Huber-White standard errors
waldtest(mod_no_iv, mod_iv, vcov = vcovHC(mod_iv, type = "HC0"))

##ROBUST STANDARD ERRORS
#requires ivpack package

iv_mod <- ivreg(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                  factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                  age_diff + factor(AGE_kid2), data = sibs)

#robust SE on treatment effect
se <- robust.se(iv_mod)[2,2]
point <- coefficients(iv_mod)[2]

#95% CIs
#for 90%, use qnorm(.95)
robust_CI <- c(point - qnorm(.975)*se, point + qnorm(.975)*se)

##CLUSTER-ROBUST STANDARD ERRORS
#requires ivpack package

#Conley and Glauber (2006) do not use clustered data, so we 
#simply illustrate the syntax for generating cluster-robust
#SEs

#iv_mod <- ivreg(outcome ~ treatment + covariate | instrument + covariate, data = data)

#se <- cluster.robust.se(iv_mod, clusterid = iv_mod$model$clusterid)[2,2]

#point <- coefficients(iv_mod)[2]

#95% CIs
#for 90%, use qnorm(.95)
#robust_CI <- c(point - qnorm(.975)*se, point + qnorm(.975)*se)

##ANDERSON-RUBIN CONFIDENCE INTERVALS

iv_mod <- ivreg(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                  factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                  age_diff + factor(AGE_kid2), data = sibs)

anderson.rubin.ci(iv_mod, conflevel = 0.95)$confidence.interval

#extracting confidence interval
ci.char <- anderson.rubin.ci(iv_mod, conflevel = 0.95)$confidence.interval

#converting it to numeric
ci.num <- as.numeric(unlist(regmatches(ci.char, gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", ci.char, perl=TRUE))))

##BOOTSTRAPPED CONFIDENCE INTERVALS
#requires boot package

iv_fun <- function(data, i) {
  dat <- data[i,]
  iv_mod <- ivreg(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                    factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                    age_diff + factor(AGE_kid2), data = dat)
  
  return(coef(iv_mod)[2])
}

#WARNING: takes a long time to run
#due to large number of obs
iv_boot <- boot(data = sibs,
                statistic = iv_fun,
                R = 10000)

#type = "perc" uses the 
#percentile bootstrap
boot.ci(iv_boot, type = "perc")

#to extract CI as vector

iv_bs_ci <- boot.ci(iv_boot, type = "perc")
iv_bs_ci$percent[4:5]

##EXPECTED TYPE-M ERROR

ols_fit <- lm(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                factor(AGE_kid2), data = sibs)

hypothetical_effect_size <- abs(ols_fit$coefficients[2])

se <- robust.se(iv_mod)[2,2]

retrodesign(A = hypothetical_effect_size,
            s = se,
            alpha = 0.05)

#type-S probability is very low (0.0001)
#expected type-M error is 1.43

##CINELLI AND HAZLETT BIAS ANALYSIS

#fitting ITT / reduced-form regression
itt_fit <- lm(priv ~ samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                factor(AGE_kid2), data = sibs)

#unmeasured confounder with partial
#R^2 sufficient to reduce point
#estimate to 0
robustness_value(itt_fit, covariates = "samesex")

#unmeasured confounder with partial
#R^2 sufficient to reduce the 
#point estimate enough for 95%
#confidence intervals to overlap 0
robustness_value(itt_fit, covariates = "samesex", alpha = 0.05)

#sensemakr contour plot for IV
sibs_sens_iv <- sensemakr(model = itt_fit,
                          treatment = "samesex",
                          benchmark_covariates = "parent_avg_educ")

#the red diamond shows
#the partial R^2 of
#parent_avg_educ
plot(sibs_sens_iv)

#sensemakr contour plot for OLS
ols_fit <- lm(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                factor(AGE_kid2), data = sibs)

sibs_sens_ols <- sensemakr(model = ols_fit,
                           treatment = "treat",
                           benchmark_covariates = "parent_avg_educ",
                           kd = 5)

#the red diamond shows
#partial R^2 values for
#an unmeasured confounder
#5x as strong as parent_ave_educ
plot(sibs_sens_ols)

##BIAS ANALYSIS PLOTS

#writing ER_sens_plot function
ER_sens_plot <- function(data,
                         formula, #make sure treatment is first variable in right-hand side
                         outcome,
                         instrument,
                         treatment,
                         theta,
                         ci_type = c("anderson-rubin", 
                                     "robust", 
                                     "clustered"),
                         clusterid) {
  
  #stop and return error msg if NULL
  assert_is_data.frame(data)
  assert_is_formula(formula)
  assert_is_vector(theta)
  assert_is_not_null(treatment)
  assert_is_not_null(outcome)
  assert_is_not_null(instrument)
  
  #extract length of theta
  len <- length(theta)
  
  #some error and warning messages
  if (len < 2) {
    stop("theta must be a vector of length > or = 2.")
  }
  
  if (missing(ci_type)) {
    warning("ci_type not specified; defaulting to \"anderson-rubin\"")
  }
  
  ci_type <- match.arg(ci_type)
  
  if (ci_type == "clustered" & missing(clusterid)) {
    stop("ci_type is set to \"clustered\", but no clusterid was specified.")
  }
  
  #updates formula to take y_adj as outcome
  formula <- update.formula(formula, y_adj ~ . )
  
  #remove parentheses from RHS to fit ivreg syntax
  formula[[3]] <- formula[[3]][[2]]
  
  #create empty vectors for sensitivity parameters
  point <- upper <- lower <- rep(NA, len)
  
  for (i in 1:len) {
    
    #adjust outcome using theta
    data <- data %>%
      mutate(y_adj = !!sym(outcome) - !!sym(instrument)*theta[i])
    
    #run ivreg with adjusted outcome
    iv_mod <- ivreg(formula = formula, data = data, x = TRUE)
    
    if (is.null(iv_mod$terms$instruments)) {
      stop("Your formula includes no instruments. Please respecify. See ?ivreg.")
    }
    
    #store adjusted point estimate
    point[i] <- coefficients(iv_mod)[2]
    
    if (ci_type == "anderson-rubin") {
      
      ci.char <- anderson.rubin.ci(iv_mod, conflevel = 0.95)$confidence.interval
      
      ci.num <- as.numeric(unlist(regmatches(ci.char, gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", ci.char, perl=TRUE))))
      
      upper[i] <- ci.num[2]
      lower[i] <- ci.num[1]
      
    } else if (ci_type == "robust") {
      
      se <- robust.se(iv_mod)[2,2]
      
      upper[i] <- point[i] + qnorm(.975)*se
      lower[i] <- point[i] - qnorm(.975)*se
      
    } else if (ci_type == "clustered") {
      
      se <- cluster.robust.se(iv_mod, clusterid = iv_mod$model$clusterid)[2,2]
      
      upper[i] <- point[i] + qnorm(.975)*se
      lower[i] <- point[i] - qnorm(.975)*se
      
    }
    
  }
  
  sens_df <- data.frame(point, upper, lower, theta)
  
  sens_plot <- ggplot(sens_df, aes(x = theta, y = point)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 2, alpha = .1) +
    xlab("Direct Effect of Instrument on Outcome (\u03b8)") +
    ylab("Adjusted Treatment Effect Estimate") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 20))
  
  ggsave("sens_plot.png")
  
  print(sens_plot)
  return(sens_df)
  
}

#calling function
#make sure treatment is first variable in right-hand side
ER_sens_plot(data = sibs,
             formula = priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
               factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
               age_diff + factor(AGE_kid2),
             outcome = "priv",
             instrument = "samesex",
             treatment = "treat",
             theta = seq(from = 0.01, to = -0.01, by = -0.0005),
             ci_type = "anderson-rubin")

##OLS VS. 2SLS COMPARISON PLOT

compare_OLS_IV <- function(data,
                           ols_formula, #make sure treatment is first variable in right-hand side
                           iv_formula, #make sure treatment is first variable in right-hand side
                           ols_ci_type = c("robust",
                                           "classical", 
                                           "clustered"),
                           iv_ci_type = c("anderson-rubin",
                                          "classical", 
                                          "robust", 
                                          "clustered"),
                           ols_cluster_id,
                           iv_cluster_id) {
  
  assert_is_data.frame(data)
  assert_is_formula(ols_formula)
  assert_is_formula(iv_formula)

  if (missing(ols_ci_type)) {
    warning("ols_ci_type not specified; defaulting to \"robust\"")
  }
  
  ols_ci_type <- match.arg(ols_ci_type)
  
  if (missing(iv_ci_type)) {
    warning("iv_ci_type not specified; defaulting to \"anderson-rubin\"")
  }
  
  iv_ci_type <- match.arg(iv_ci_type)
  
  if (ols_ci_type == "clustered" & missing(ols_cluster_id)) {
    stop("ols_ci_type is set to \"clustered\", but no ols_cluster_id was specified.")
  }
  
  if (iv_ci_type == "clustered" & missing(iv_cluster_id)) {
    stop("iv_ci_type is set to \"clustered\", but no iv_cluster_id was specified.")
  }
  
  est <- data.frame(estimator = c("OLS", "2SLS"), 
                    point = c(NA, NA),
                    upper = c(NA, NA),
                    lower = c(NA, NA))
  
  iv_mod <- ivreg(formula = iv_formula, data = data, x = TRUE)
  ols_mod <- lm(formula = ols_formula, data = data)
  
  ols_point <- est$point[1] <- ols_mod$coefficients[2]
  iv_point <- est$point[2] <- iv_mod$coefficients[2]
  
  #OLS CI
  
  if (ols_ci_type == "robust") {
  
    se <- coeftest(ols_mod, vcov = vcovHC(ols_mod, type="HC0"))[2,2]
  
    est$upper[1] <- ols_point + qnorm(.975)*se
    est$lower[1] <- ols_point - qnorm(.975)*se
  
  } else if (ols_ci_type == "classical") {
    
    se <- summary(ols_mod)$coefficients[2,2]
    
    est$upper[1] <- ols_point + qnorm(.975)*se
    est$lower[1] <- ols_point - qnorm(.975)*se
    
  } else if (ols_ci_type == "clustered") {
  
    se <- coeftest(ols_mod, vcov = vcovCL, cluster = ~ols_cluster_id)[2,2]
  
    st$upper[1] <- ols_point + qnorm(.975)*se
    est$lower[1] <- ols_point - qnorm(.975)*se
    
  }  
    
  #2SLS CI
  
  if (iv_ci_type == "anderson-rubin") {
  
    ci.char <- anderson.rubin.ci(iv_mod, conflevel = 0.95)$confidence.interval
    ci.num <- as.numeric(unlist(regmatches(ci.char, gregexpr("(?>-)*[[:digit:]]+\\.*[[:digit:]]*", ci.char, perl=TRUE))))
  
    est$upper[2] <- ci.num[2]
    est$lower[2] <- ci.num[1]
  
  } else if (iv_ci_type == "robust") {
    
    se <- robust.se(iv_mod)[2,2]
    
    upper[2] <- point[2] + qnorm(.975)*se
    lower[2] <- point[2] - qnorm(.975)*se
    
  } else if (iv_ci_type == "clustered") {
    
    se <- cluster.robust.se(iv_mod, clusterid = iv_mod$model$iv_cluster_id)[2,2]
      
    upper[2] <- point[2] + qnorm(.975)*se
    lower[2] <- point[2] - qnorm(.975)*se
      
  } else if (iv_ci_type == "classical") {
   
    se <- summary(iv_mod)$coefficients[2,2]
      
    upper[2] <- point[2] + qnorm(.975)*se
    lower[2] <- point[2] - qnorm(.975)*se
    
  }
  
  est_plot <- ggplot(data = est, aes(estimator, point)) +
    geom_point(size = 5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 2) +
    xlab("Estimator") +
    ylab("Estimated Treatment Effect") +
    ylim(-.10, 0.05) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 20))
  
  ggsave("est_plot.png")
  
  print(est_plot)
  return(est)
  
}

#calling function
#make sure treatment is first variable in right-hand side
compare_OLS_IV(data = sibs,
               ols_formula = priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                 factor(AGE_kid2),
               iv_formula = priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                 factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                 age_diff + factor(AGE_kid2),
               ols_ci_type = "robust",
               iv_ci_type = "anderson-rubin")

