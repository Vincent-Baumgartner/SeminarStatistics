rm(list = ls())
setwd("C:/Uni Statistik/Bachelorarbeit")

library(tidyverse)
library(mice)
source("EM_Regression.R")
set.seed(42)
#-------------------------------------------------------------------------------
#                             Simulations
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# (1)                   [Simple Linear Regression]
#                              
#                        y_i = a + b * x_i + eps_i
#-------------------------------------------------------------------------------

n <- 1000       # sample size
a <- 0       # intercept
b <- 0.25       # slope
sigma_y2 <- 50  # Var(eps_i)
mu_x <- 0      # E[x_i]
sigma_x2 <- 500 # Var(x_i)

## complete Data simulation
x <- rnorm(n, mu_x, sqrt(sigma_x2))
y <- a+ x * b + rnorm(n, 0, sqrt(sigma_y2))

## Missing Data Mechanism by logistic function
by <- 1
ay <- -a - b*mu_x

p_miss <- plogis(ay + by * y)
mean(p_miss)

R <- rbinom(n, size = 1, prob = p_miss)
mean(R)
x_mis <- x
x_mis[R == 1] <- NA

df_mis <- data.frame(y = y, x_mis = x_mis)
md.pattern(df_mis)
df <- data.frame(y, x, mis = is.na(x_mis))

## complete case analysis
lm_model <- lm(y ~., data = df_mis)
summary(lm_model)

## em algorithm
em_model <- lm_em(y, as.matrix(x_mis), verbose = TRUE)
em_model


p <- df |>
  ggplot(aes(x = x, y = y)) +
  geom_point(aes(alpha = mis), color = ifelse(df$mis, "grey50", "black"), 
             pch = 20, size = 5) +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.2), guide = "none") +
  geom_abline(aes(intercept = a, slope = b, color = "True Line"), 
              linewidth = 2) +
  geom_abline(aes(intercept = coef(lm_model)[1], slope = coef(lm_model)[2], 
                  color = "CC Regression"), linewidth = 2) +
  geom_abline(aes(intercept = em_model$alpha, slope = em_model$beta,
                  color = "EM Regression"), linewidth = 2) +
  scale_colour_manual(name = "Line type", 
                      values = c("True Line"       = "#09622A",
                                 "CC Regression"   = "#9C0824",
                                 "EM Regression"   = "#2E5A87")) +
  theme_minimal() 
p

# Monte Carlo sim:
n_sim <- 5000
one_sim <- function(n = 100, a = 0, b = 0.25, sy2 = 50, mux = 50 , sx2 = 500, 
                    by = 1) {
  
  x <- rnorm(n, mu_x, sqrt(sx2))
  y <- a + b * x + rnorm(n, 0, sqrt(sy2))
  
  ay <- -a - b*mu_x
  p_miss <- plogis(ay + by * y)
  R <- rbinom(n, 1, p_miss)
  x_mis <- x; x_mis[R == 1] <- NA
  
  df_mis <- data.frame(y = y, x_mis = x_mis)
  
  # CC‐Fit
  cc_fit   <- lm(y ~ x_mis, data = df_mis)
  cc_int   <- coef(cc_fit)[1]
  cc_slope <- coef(cc_fit)[2]
  
  # EM‐Fit
  em_out   <- lm_em(y, as.matrix(x_mis), tol = 1e-6, verbose = TRUE)
  em_int   <- em_out$alpha
  em_slope <- em_out$beta
  
  return(c(cc_int, cc_slope, em_int, em_slope))
}
res_mat <- replicate(n_sim, one_sim(), simplify = TRUE)
res_df <- as_tibble(t(res_mat),
                    .name_repair = "unique") %>%
  set_names(c("cc_int", "cc_slope", "em_int", "em_slope"))

summary_df <- tibble(
  method    = c("CC", "CC", "EM", "EM"),
  param     = c("intercept", "slope", "intercept", "slope"),
  estimates = list(res_df$cc_int, res_df$cc_slope,
                   res_df$em_int, res_df$em_slope)
) %>%
  mutate(
    true = if_else(param == "intercept", a, b),
    bias = map2_dbl(estimates, true, ~ mean(.x - .y)),
    rmse = map2_dbl(estimates, true, ~ sqrt(mean((.x - .y)^2)))
  ) %>%
  select(method, param, true, bias, rmse)

summary_df
write.csv(summary_df,
          file = "summary_df.csv",
          row.names = FALSE)
read.csv("summary_df.csv")
