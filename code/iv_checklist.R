# Below, we provide example R code for generating the sensitivity plots shown in the main text. 
# We name the dataset sibs. The outcome (private school attendance) is named priv. 
# The treatment (having two siblings) is named treat. The instrument (having the first two children be of the same sex) is 
# named samesex. All other variables are additional covariates included in their model specifications.

## INSTALLING AND LOADING NECESSARY PACKAGES

#install.packages("lmtest")
#install.packages("sandwich")
#install.packages("ivmodel")
#install.packages("boot")
#install.packages("assertive")
#install.packages("rlang")
#install.packages("sensemakr")
#install.packages("retrodesign")
#install.packages("estimatr")

library(lmtest)
library(sandwich)
library(ivmodel)
library(boot)
library(assertive)
library(rlang)
library(sensemakr)
library(retrodesign)
library(estimatr)

## READING IN DATA

sibs <- readRDS(here("data", "cg2006.rds"))

## ROBUST F-TEST
## Requires lmtest package

# First-stage model WITHOUT instrument
fs_mod_no_iv <- lm(treat ~ parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                  age_diff + factor(AGE_kid2), data = sibs)

# First-stage model WITH instrument
fs_mod <- lm(treat ~ samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
               age_diff + factor(AGE_kid2), data = sibs)

# Fartial F-test using Huber-White standard errors
waldtest(fs_mod_no_iv, fs_mod, vcov = vcovHC(fs_mod, type = "HC2"))

# Other HC options produce extremely similar results

#waldtest(fs_mod_no_iv, fs_mod, vcov = vcovHC(fs_mod, type = "HC0"))
#waldtest(fs_mod_no_iv, fs_mod, vcov = vcovHC(fs_mod, type = "HC1"))
#waldtest(fs_mod_no_iv, fs_mod, vcov = vcovHC(fs_mod, type = "HC3"))

## HETEROSKEDASTICITY-ROBUST CONFIDENCE INTERVALS
iv_mod <- iv_robust(formula = priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                      factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                      age_diff + factor(AGE_kid2), 
                    data = sibs,
                    se_type = "HC2")

# Full regression output
iv_mod

# Point estimate
iv_mod$coefficients[2]

# Confidence interval based on robust Huber-White HC2 SEs
c(iv_mod$conf.low[2], iv_mod$conf.high[2])

## CLUSTER-ROBUST STANDARD ERRORS
## Requires estimatr package

# Conley and Glauber (2006) do not use clustered data, so we 
# simply illustrate the syntax for generating cluster-robust
# SEs

# iv_mod <- iv_robust(outcome ~ treatment + covariate | instrument + covariate, 
#                     data = data,
#                     clusters = cluster_varname)

# iv_mod can be examined for point estimates / CIs

## ANDERSON-RUBIN CONFIDENCE INTERVALS
## Requires ivmodel package

iv_mod <- ivmodelFormula(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                 factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                 age_diff + factor(AGE_kid2),
               data = sibs,
               heteroSE = F,
               deltarange = NULL)

iv_mod$AR$ci[2]

## BOOTSTRAPPED CONFIDENCE INTERVALS
## Requires boot package

iv_fun <- function(data, i) {
  dat <- data[i,]
  iv_mod <- ivreg(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                    factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                    age_diff + factor(AGE_kid2), data = dat)
  
  return(coef(iv_mod)[2])
}

# WARNING: takes a long time to run
# due to large number of obs
iv_boot <- boot(data = sibs,
                statistic = iv_fun,
                R = 1000)

# type = "perc" uses the 
# percentile bootstrap
boot.ci(iv_boot, type = "perc")

# to extract CI as vector
iv_bs_ci <- boot.ci(iv_boot, type = "perc")
iv_bs_ci$percent[4:5]

## CLASSICAL CONFIDENCE INTERVALS

iv_mod <- ivreg(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                  factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                  age_diff + factor(AGE_kid2), data = sibs)

est <- coef(iv_mod)[2]

se <- sqrt(diag(vcov(iv_mod)))[2]

est + qnorm(.975)*se
est - qnorm(.975)*se


## EXPECTED TYPE-M ERROR

# Here we just use the OLS estimate
# for the hypothetical effect size

ols_fit <- lm(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                factor(AGE_kid2), data = sibs)

hypothetical_effect_size <- abs(ols_fit$coefficients[2])

retrodesign(A = hypothetical_effect_size,
            s = se,
            alpha = 0.05)

# Type-S probability is very low (0.0001)
# Expected type-M error is 1.43

## CINELLI AND HAZLETT BIAS ANALYSIS

# Fitting ITT / reduced-form regression
itt_fit <- lm(priv ~ samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                factor(AGE_kid2), data = sibs)

# Unmeasured confounder with partial
# R^2 sufficient to reduce point
# estimate to 0
robustness_value(itt_fit, covariates = "samesex")

# Unmeasured confounder with partial
# R^2 sufficient to reduce the 
# point estimate enough for 95%
# confidence intervals to overlap 0
robustness_value(itt_fit, covariates = "samesex", alpha = 0.05)

# sensemakr contour plot for IV
sibs_sens_iv <- sensemakr(model = itt_fit,
                          treatment = "samesex",
                          benchmark_covariates = "parent_avg_educ",
                          kd = -1)

pdf(here("plots", "ivcontour.pdf"), 
    height = 6.38,
    width = 6.38,
    pointsize = 20)

par(mgp = c(3,.5,0)) # moving tick labels closer to ticks
plot(sibs_sens_iv,
     nlevels = 9, # number of contours
     lim = 0.015, # x- and -axis limits
     lim.y = 0.015,
     cex.axis = .6, # relative size of axis labels
     xlab = "", 
     ylab = "",
     label.text = F, # don't include text for adjusted estimate
     labcex = .5, # change contour text size (0.005, 0.004, etc)
     show.unadjusted = F, # don't show unadjusted estimate
     list.par = list(mar = c(2,2,1,1), pty = "s")) # change plot margins
points(0, 0, pch = 17, col = "black", cex = 1) # manually plotting unadjusted estimate

#saving plot
dev.off() 

## If you set show.unadjusted = T (the default), it will automatically
## run the points() code but include a label. I set it to F and removed
## the label to add my own in LaTeX. 

## To extract info from contour plot to make your own:
#
# iv_sens_plot <- plot(sibs_sens_iv)
#
# contour(iv_sens_plot$r2dz.x, iv_sens_plot$r2yz.dx, iv_sens_plot$value)

## BIAS ANALYSIS PLOTS based on Conley et al. (2012) and Wang et al. (2018)

# Function for getting df for plot
ER_sens_df <- function(data,
                       formula, # make sure treatment is first variable in right-hand side
                       outcome,
                       instrument,
                       treatment,
                       theta, # range of plausible values of exclusion restriction violations
                       ci_type = c("anderson-rubin", 
                                   "robust",
                                   "classical")) {
  
  # stop and return error msg if NULL
  assert_is_data.frame(data)
  assert_is_formula(formula)
  assert_is_vector(theta)
  assert_is_not_null(treatment)
  assert_is_not_null(outcome)
  assert_is_not_null(instrument)
  
  # extract length of theta
  len <- length(theta)
  
  # some error and warning messages
  if (len < 2) {
    stop("theta must be a vector of length > or = 2.")
  }
  
  if (missing(ci_type)) {
    warning("ci_type not specified; defaulting to \"anderson-rubin\"")
  }
  
  ci_type <- match.arg(ci_type)
  
  # updates formula to take y_adj as outcome
  formula <- update.formula(formula, y_adj ~ . )
  
  # remove parentheses from RHS to fit ivreg syntax
  formula[[3]] <- formula[[3]][[2]]
  
  # create empty vectors for sensitivity parameters
  point <- upper <- lower <- rep(NA, len)
  
  for (i in 1:len) {
    
    # adjust outcome using theta
    data <- data %>%
      mutate(y_adj = !!sym(outcome) - !!sym(instrument)*theta[i])
    
    if (ci_type == "anderson-rubin") {
      
      iv_mod <- ivmodelFormula(formula = formula,
                               data = data,
                               heteroSE = F,
                               deltarange = NULL)
      
      point[i] <- coef(iv_mod)[4,2]
      upper[i] <- iv_mod$AR$ci[2]
      lower[i] <- iv_mod$AR$ci[1]
      
    } else if (ci_type == "robust") {
      
      iv_mod <- iv_robust(formula = formula, 
                          data = data,
                          se_type = "HC2")
      
      point[i] <- iv_mod$coefficients[2]
      lower[i] <- iv_mod$conf.low[2]
      upper[i] <- iv_mod$conf.high[2]
      
    } else if (ci_type == "classical") {
      
      iv_mod <- ivreg(formula = formula,
                      data = data)
      
      point[i] <- coef(iv_mod)[2]
      se <- sqrt(diag(vcov(iv_mod)))[2]
      
      upper[i] <- point[i] + qnorm(.975)*se
      lower[i] <- point[i] - qnorm(.975)*se
      
    }
    
  }
  
  sens_df <- data.frame(point, upper, lower, theta)
  return(sens_df)
  
}

# Calling function
# Make sure treatment is first variable in right-hand side
sens_df <- ER_sens_df(data = sibs,
             formula = priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
               factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
               age_diff + factor(AGE_kid2),
             outcome = "priv",
             instrument = "samesex",
             treatment = "treat",
             theta = seq(from = 0.01, to = -0.01, by = -0.0005),
             ci_type = "anderson-rubin")

ITT_est <- itt_fit$coefficients[2]

ggplot(sens_df, aes(x = theta, y = point)) +
  geom_line(linewidth = .7) +
  geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 2, alpha = .1) +
  #xlab("Direct Effect of Instrument on Outcome (\u03b8)") + # Commenting out because we're doing axis labels in TikZ
  #ylab("Adjusted Treatment Effect Estimate") +
  xlab("") +
  ylab("") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = .67) + 
  geom_vline(xintercept = ITT_est, color = "red", linetype = "dashed", linewidth = .67) + 
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20)) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = c(-0.008, -0.004, 0, 0.004, 0.008)) 

ggsave(here("plots", "sens_plot.png"), height = 5, width = 8, unit = "in")

## 2SLS CI COMPARISON PLOT

# Function for getting df for plot
compare_IV_CI <- function(data,
                          formula # Make sure treatment is first variable in right-hand side
                          ) {
  
  assert_is_data.frame(data)
  assert_is_formula(formula)

  ci_df <- data.frame(est = rep(NA, 4),
                      ci_type = c("anderson-rubin", "robust", "bootstrap", "classical"),
                      upper = rep(NA, 4),
                      lower = rep(NA, 4))
  
  iv_mod <- ivreg(formula = formula, data = data, x = TRUE)
  
  # Point estimate
  est <- coef(iv_mod)[2]
  ci_df$est <- rep(est, 4)
  
  # Classical CIs
  
  se <- sqrt(diag(vcov(iv_mod)))[2]
  
  ci_df$upper[4] <- as.numeric(est + qnorm(.975)*se)
  ci_df$lower[4] <- as.numeric(est - qnorm(.975)*se)
    
  # AR CI
  
  iv_mod <- ivmodelFormula(formula = formula,
                           data = data,
                           heteroSE = F,
                           deltarange = NULL)
  
  ci_df$upper[1] <- iv_mod$AR$ci[2]
  ci_df$lower[1] <- iv_mod$AR$ci[1]
  
  # Robust CI
  
  iv_mod <- iv_robust(formula = formula, 
                      data = data,
                      se_type = "HC2")
  
  ci_df$lower[2] <- iv_mod$conf.low[2]
  ci_df$upper[2] <- iv_mod$conf.high[2]
  
  # Bootstrap CI
  
  iv_fun <- function(data, i) {
    dat <- data[i,]
    iv_mod <- ivreg(formula = formula, data = dat)
    
    return(coef(iv_mod)[2])
  }
  
  iv_boot <- boot(data = data,
                  statistic = iv_fun,
                  R = 1000)
  
  iv_bs_ci <- boot.ci(iv_boot, type = "perc")
  ci_df$lower[3] <- iv_bs_ci$percent[4]
  ci_df$upper[3] <- iv_bs_ci$percent[5]
  
  return(ci_df)
  
}

ci_for_plot_df <- compare_IV_CI(data = sibs,
                                formula = priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
                                  factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
                                  age_diff + factor(AGE_kid2))

ggplot(data = ci_for_plot_df, aes(ci_type, est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), size = .7, width = 0) +
  xlab("") +
  ylab("Effect on Probability of Attending Private School") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 20)) +
  scale_x_discrete(labels = c("Anderson\u2013Rubin", "Robust", "Bootstrapped", "Classical")) #\u2013 creates an en-dash
  