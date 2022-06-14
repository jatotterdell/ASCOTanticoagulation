# Simulate a simple trial with an ordinal outcome.

library('rmsb')
library('tidyverse')
library(tidyverse)
library(labelled)
library(kableExtra)
library(cmdstanr)
library(posterior)
library(bayestestR)
library(bayesplot)
library(matrixStats)

options(mc.cores = parallel::detectCores())   # use max # CPUs

# ==========================================================================================
## Compute the mean and median x after shifting the probability
## distribution by an odds ratio under the proportional odds model
states <- who_labels <-c("Not hospitalised, no limitations on activities",
                         "Not hospitalised, limitation on activities",
                         "Hospitalised, not requiring supplemental oxygen and no longer requiring ongoing medical care ",
                         "Hospitalised, not requiring supplemental oxygen but requiring ongoing medical care",
                         "Hospitalised, requiring supplemental oxygen",
                         "Hospitalised, on non-invasive ventilation or high flow oxygen devices",
                         "Hospitalised, on invasive mechanical ventilation or ECMO",
                         "Death")
iStates <- seq(states)

# ==========================================================================================
# Multinomial simulation.

simulateTrial <- function()
{
  ## Define outcome probabilities for control group
  p0 <- c(416, 146, 5, 6, 2, 0, 2, 19)
  p0 <- p0/sum(p0)

  ## Define outcome probabilities for the treatment groups
  p1   <- pomodm(p = p0, odds.ratio = 1/1.05)
  p2   <- pomodm(p = p0, odds.ratio = 1/1.05)
  p3   <- pomodm(p = p0, odds.ratio = 1/1.05)

  # Simulate two arm trial data.
  dMulti0  <- rmultinom(1, size = 600, prob = p0)
  dMulti1  <- rmultinom(1, size = 600, prob = p1)
  dMulti2  <- rmultinom(1, size = 300, prob = p2)
  dMulti3  <- rmultinom(1, size = 50,  prob = p3)

  sample0  <- rep(states, dMulti0) %>% factor(ordered = iStates, levels = states)
  sample1  <- rep(states, dMulti1) %>% factor(ordered = iStates, levels = states)
  sample2  <- rep(states, dMulti2) %>% factor(ordered = iStates, levels = states)
  sample3  <- rep(states, dMulti3) %>% factor(ordered = iStates, levels = states)

  # Munge simulated data.
  data           <- rbind(cbind(sample0, 0),
                          cbind(sample1, 1),
                          cbind(sample2, 2),
                          cbind(sample3, 3))
  data           <- data.frame(data)
  colnames(data) <- c("y", "CAssignment")
  data[["CAssignment"]] <- as.factor(data[["CAssignment"]])
  data
}

make_sim_data <- function(dat, vars, outcome)
{
  X <- model.matrix(
    as.formula(paste("~", paste(vars, collapse = " + "))),
    data = dat)

  nXassign <- sum(grepl("Assignment", colnames(X)))

  X <- as.matrix(X[,-1])

  y <- dat[[outcome]]
  N <- dim(X)[1]
  K <- dim(X)[2]
  J <- max(y)

  list(N = N, K = K, X = X, y = y, J = J, p_par = rep(1, J),
       beta_sd = rep(1, nXassign))
}

data <- simulateTrial()

with(data, table(y, CAssignment))

stan_data <- make_sim_data(data, "CAssignment", "y")

run_stan <- function(model, data)
{
  model[["sample"]](data = data,
                    seed = seed,
                    refresh = 0,
                    iter_sampling = 2500,
                    chains = 8,
                    parallel_chains = 8)
}

seed <- 59579
model_simple <- cmdstan_model(file.path("stan", "logistic_cumulative_ordinal.stan"))

results <- run_stan(model_simple, stan_data)

results$summary(variables = c("p",
                              "alpha",
                              "beta"))

# ==========================================================================================
# Save things.

# save(file = file.path("..", "..", "output", "simList.Rdata"), simList)


# End of script.

