library(ggplot2)
library(gtools)
library(dplyr)
library(tidyr)


## We tried 10 differents seed from 2024 to 2034.
set.seed(2031)


## Parameters
## group 0 (g0) corresponds to high genetic risk score
N1 <- 16
N2 <- 16
tau_g0 <- 1
tau_g1 <- 0
tau_seq <- seq(-1.5, 2.5, by = 0.05)


## Data generation
X1 <- c(rep(0, N1/2), rep(1, N1/2))
X2 <- c(rep(0, N2/2), rep(1, N2/2))
Y_0 <- rnorm(N1+N2, 0, 1)
Y_1 <- Y_0 + tau_g0 * (1-c(X1, X2)) + tau_g1 * c(X1, X2)

treat <- function(n) sample(c(rep(0, n/2), rep(1, n/2)))
Z1 <- c(treat(N1/2), treat(N1/2))
Z2 <- c(treat(N2/2), treat(N2/2))
Y1 <- Y_0[1:N1] * (1-Z1) + Y_1[1:N1] * Z1
Y2 <- Y_0[(N1+1):N2] * (1-Z2) + Y_1[(N1+1):N2] * Z2


## Statistics for selection and hypothesis test
SATE <- function(y, z) {
  mean_diff <- mean(y[z==1]) - mean(y[z==0])
  sd_est <- sqrt(var(y[z==1])/sum(z) + var(y[z==0])/sum(1-z))
  mean_diff / sd_est
}

Delta <- function(Y, X, Z) {
  sate_g0 <- SATE(Y[X==0], Z[X==0])
  sate_g1 <- SATE(Y[X==1], Z[X==1])
  (sate_g1 - sate_g0) / sqrt(2)
}


delta_holdout <- Delta(Y1, X1, Z1)
delta_full <- Delta(c(Y1, Y2), c(X1,X2), c(Z1, Z2))


## Only proceed if both selection rules lead to subgroup "high"
if (delta_holdout < qnorm(0.2) && delta_full < qnorm(0.2)) {
  ## Defining short-hands
  Y1_g0 <- Y1[X1==0]
  Y2_g0 <- Y2[X2==0]
  Y_g0 <- c(Y1[X1==0], Y2[X2==0])
  Y1_g1 <- Y1[X1==1]
  Y2_g1 <- Y2[X2==1]
  Z1_g0 <- Z1[X1==0]
  Z2_g0 <- Z2[X2==0]
  Z_g0 <- c(Z1[X1==0], Z2[X2==0])
  Z1_g1 <- Z1[X1==1]
  Z2_g1 <- Z2[X2==1]

  ## observed value of test statistic
  t_obs <- SATE(c(Y1_g0, Y2_g0), c(Z1_g0, Z2_g0))

  ## generating all different combinations of treated people
  comb1 <- combinations(n = N1/2, r = N1/4)
  comb2 <- combinations(n = N2/2, r = N2/4, v=N1/2 + (1:(N2/2)))
  comb <- merge(as.data.frame(comb1), as.data.frame(comb2), by = NULL)


  ## Computing p-value for specific tau
  generate_pval <- function(tau) {

    ## Computing selection and test statistics when people with indices a
    ## are treated
    helper <- function(a) {
      ## treatment assignments and imputed outcomes
      z <- rep(0, (N1+N2)/2) ## treatment assignments in group 0
      z[a] <- 1
      Y_imp <- Y_g0 + (z - Z_g0) * tau
      Y_imp1 <- Y_imp[1:(N1/2)]
      Y_imp2 <- Y_imp[(N1/2) + (1:(N2/2))]
      z1 <- z[1:(N1/2)]
      z2 <- z[(N1/2) + (1:(N2/2))]

      ## Selection statistics for 3 scenarios:
      ## holdout units, no holdout units, second stage only RT
      sel_holdout <- Delta(c(Y_imp1, Y1_g1), X1, c(z1, Z1_g1)) < qnorm(0.2)
      sel_full <- Delta(c(Y_imp1, Y1_g1, Y_imp2, Y2_g1),
                        c(X1, X2),
                        c(z1, Z1_g1, z2, Z2_g1)) < qnorm(0.2)
      sel_second <- all(z1 == Z1_g0)

      ## value of test statistic under treatment z
      test_stat <- SATE(Y_imp, z)

      c(sel_holdout, sel_full, sel_second, test_stat)
    }

    ## assembling values of statistics for all feasible treatment assignments
    res_mat <- t(apply(comb, 1, helper))
    res_df <- as.data.frame(res_mat)
    colnames(res_df) <- c("sel_holdout", "sel_full", "sel_second", "test_stat")
    res_df <- mutate(res_df, smaller = test_stat > t_obs)

    ## Computing different p-values
    pval_holdout <- filter(res_df, sel_holdout == 1) |>
      summarize(pval_holdout = mean(smaller)) |>
      pull()

    pval_full <- filter(res_df, sel_full == 1) |>
      summarize(pval_holdout = mean(smaller)) |>
      pull()

    pval_second <- filter(res_df, sel_second == 1) |>
      summarize(pval_second = mean(smaller)) |>
      pull()

    c(pval_holdout, pval_full, pval_second)
  }


  ## Computing and storing p-value curve over sequences of tau's
  pval_mat <- t(sapply(tau_seq, generate_pval))
  pval_df <- as.data.frame(pval_mat)
  colnames(pval_df) <- c("holdout", "full", "second")


  ## Helper functions for visualization
  compute_mid <- function(a, b, p) {
    x <- lead(a)
    y <- lead(b)
    lambda <- (p-b) / (y-b)
    a * (1-lambda) + x * lambda
  }

  cross_fun <- function(a, b) {
    (a-b) * (lead(a)-b) < 0
  }

  touch_fun <- function(a, b) {
    abs(a - b) < 0.001
  }


  ## Levels for estimation and 90% confidence interval
  cross_val <- 0.5
  confint_val <- 0.1


  ## Preparing data for plots
  est_confint_df <- pval_df |>
    mutate(tau = tau_seq) |>
    mutate(across(1:3, list(tau.midAest = ~compute_mid(tau, ., cross_val),
                            touchAest = ~touch_fun(., cross_val),
                            crossAest = ~cross_fun(., cross_val),
                            tau.midAci = ~compute_mid(tau, ., confint_val),
                            touchAci = ~touch_fun(., confint_val),
                            crossAci = ~cross_fun(., confint_val)),
                  .names = "{.fn}_{.col}")) |>
    rename(pval_holdout = holdout, pval_full = full, pval_second = second) |>
    pivot_longer(cols = -tau,
                 names_to = c(".value", "kind"),
                 names_sep = "_") |>
    pivot_longer(cols = -(1:3),
                 names_to = c(".value", "type"),
                 names_sep = "A") |>
    filter(touch | cross) |>
    mutate(x = if_else(cross, tau.mid, tau),
           y = if_else(type == "est", cross_val, confint_val)) |>
    select(x, y, tau, kind, type)


  pval_df <- pval_df |>
    mutate(tau = tau_seq) |>
    pivot_longer(cols = -tau,
                 names_to = "kind",
                 values_to = "pval")


  join_df2 <- pval_df |>
    left_join(est_confint_df, by = c("tau" = "tau", "kind" = "kind"))

  name_list <- list("full" = "without hold-out units",
                    "holdout" = "with hold-out units",
                    "second" = expression(2^nd ~ stage ~ RT))
  labeller_fun <- function(variable, value){
    name_list[value]
  }


  ## Plotting the results
  ggplot(data = join_df2) +
    geom_hline(yintercept = c(0.1,0.5), col = "black",
               linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = 1, col = "gray",
               linetype = "solid", linewidth = 0.3) +
    geom_segment(aes(x = x, y = y, xend = x, yend = -Inf,
                     col = type),
                 linewidth = 0.5,
                 alpha = 1,
                 na.rm = TRUE
    ) +
    geom_line(aes(x = tau, y = pval),
              linewidth = 0.7,
              col = "#009E73",
              alpha = 1) +
    facet_grid(~factor(kind, levels=c("full", "holdout", "second")),
               labeller = labeller_fun) +
    coord_cartesian(xlim = c(-1.5, 2.5)) +
    labs(x = expression(tau),
         y = expression(P(tau))) +
    scale_colour_discrete(name = "",
                          breaks = c("ci", "est"),
                          labels = c("Lower end of 90%\nconfidence interval",
                                     "Treatment effect\nestimate"),
                          na.translate = F) +
    theme_bw() +
    theme(strip.text.x = element_text(size = 11),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.height = unit(1.2,"cm"))

} else {
  print("Choose another seed.")
}
