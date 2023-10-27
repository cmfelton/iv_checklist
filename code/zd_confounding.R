library(AER)
library(tidyverse)

rm(list = ls())

n.obs <- 5000
n.sims <- 1000

## WITH dz confounding

sims <- data.frame(tau_est = rep(NA, n.sims),
                   f_test = rep(NA, n.sims))

set.seed(111)
for (i in 1:n.sims) {

  uzd <- rbinom(n.obs, 1, .5)
  udy <- rbinom(n.obs, 1, .5)

  dat <- data.frame(uzd, udy) %>%
    mutate(z_prob = case_when(
                            uzd == 1 ~ .75,
                            uzd == 0 ~ .25
                           ),
        z = rbinom(n.obs, 1, prob = z_prob),
        d0_prob = case_when(
                            uzd == 1 & udy == 1 ~ .9,
                            uzd == 0 & udy == 0 ~ .1,
                            uzd == 1 & udy == 0 ~ .7,
                            uzd == 0 & udy == 1 ~ .4
                           ),
        d0 = rbinom(n.obs, 1, prob = d0_prob),
        comp_ind = rbinom(n.obs, 1, prob = .5),
        comp = case_when(
                         d0 == 0 & comp_ind == 1 ~ "C",
                         d0 == 0 & comp_ind == 0 ~ "NT",
                         d0 == 1 ~ "AT"
                        ),
        d1 = case_when(
                       comp == "AT" | comp == "NT" ~ d0,
                       comp == "C" ~ 1
                      ),
        d = case_when(
                      comp == "C" & z == 1 ~ d1,
                      .default = d0
                     ),
        tau = case_when(
                        comp == "C" ~ .9, # LATE
                        .default = 0
                       ),
        y = rnorm(n.obs, udy*.5 + d*tau))

  for_f <- summary(lm(d ~ z, dat))

  iv_est <- ivreg(y ~ d | z, data = dat)

  sims$tau_est_iv[i] <- iv_est$coefficients[2]

  ols_est <- lm(y ~ d, data = dat)
  
  sims$tau_est_ols[i] <- ols_est$coefficients[2]
  
  sims$f_test[i] <- for_f$fstatistic[1]

}

# LATE is .9

mean(sims$tau_est_iv) # 0.3432843
mean(sims$tau_est_ols) # 0.2270462
mean(sims$f_test) # 1310.655, so estimation bias is negligible

## WITHOUT dz confounding

sims <- data.frame(tau_est = rep(NA, n.sims),
                   f_test = rep(NA, n.sims))

set.seed(111)
for (i in 1:n.sims) {
  
  udy <- rbinom(n.obs, 1, .5)
  
  dat <- data.frame(udy) %>%
    mutate(z_prob = case_when(
      uzd == 1 ~ .75,
      uzd == 0 ~ .25
    ),
    z = rbinom(n.obs, 1, prob = z_prob),
    d0_prob = case_when(
      udy == 1 ~ .8,
      udy == 0 ~ .2
    ),
    d0 = rbinom(n.obs, 1, prob = d0_prob),
    comp_ind = rbinom(n.obs, 1, prob = .5),
    comp = case_when(
      d0 == 0 & comp_ind == 1 ~ "C",
      d0 == 0 & comp_ind == 0 ~ "NT",
      d0 == 1 ~ "AT"
    ),
    d1 = case_when(
      comp == "AT" | comp == "NT" ~ d0,
      comp == "C" ~ 1
    ),
    d = case_when(
      comp == "C" & z == 1 ~ d1,
      .default = d0
    ),
    tau = case_when(
      comp == "C" ~ .9, # LATE
      .default = 0
    ),
    y = rnorm(n.obs, udy*.5 + d*tau))
  
  for_f <- summary(lm(d ~ z, dat))
  
  iv_est <- ivreg(y ~ d | z, data = dat)
  
  sims$tau_est_iv[i] <- iv_est$coefficients[2]
  
  ols_est <- lm(y ~ d, data = dat)
  
  sims$tau_est_ols[i] <- ols_est$coefficients[2]
  
  sims$f_test[i] <- for_f$fstatistic[1]
  
}

# LATE is .9

mean(sims$tau_est_iv) # 0.897298
mean(sims$tau_est_ols) # 0.4179464
mean(sims$f_test) # 359.2909, so estimation bias is negligible
range(sims$f_test) # 223.0412 527.0398