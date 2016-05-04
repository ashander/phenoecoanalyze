library(phenoecosim)
context("Analytical functions")


context("Approximate predictions")
param_string <- '{"K0":100000000,
"omegaz2":25,
"omega_Wmax":1.15,
"opt_Wmax":0,
"A":0,
"B":2,
"sigma_xi2":0.5,
"fractgen":5,
"rho_gen":false,
"rho_dev_sel":0.8,
"delta":2.5,
"alpha":0.5,
"var_a":0.5,
"Vb":0.005,
"Cab":0,
"var_ab":0,
"seed":111,
"Ve":0.5,
"Npop0": 10000,
"Nc": 500,
"t_vals": [5, 10, 20, 40, 80, 100, 150, 200],
"Tlim": 200,
"nrep": 500,
"interp_alpha": false,
"predict_alpha": false,
"exact": false}'


test_that("(cross-package) Nt pred approx equal to median Gamma dens", {
try(library(jsonlite, quietly=TRUE, warn.conflicts=FALSE),
      silent = TRUE)
  try(params  <- fromJSON(param_string)) #stream_in needs minified json so use fromJSON
  if(!exists('params'))
    stop("need R package jsonlite")

  invisible(list2env(params, envir=environment()))


  ## @knitr run-sims
  tau <- 1 / fractgen

  if (rho_gen) {
    rho <- rho_gen
  } else {
    rho <- rho_dev_sel
  }

  params_dep <- list(tau = tau,
                     ## reaction norm elevation in the reference environment
                     ## and plsticity must be uncorrelated
                     Gmat=G(var_a, Vb, var_ab),
                     R0 = omega_Wmax,
                     rho_dev_sel = rho,
                     rho_sel = rho ^ (1 / tau),
                     rho_dev = rho ^ (1 / tau),
                     rho_dev_selp = rho ^ ((1 + tau) / tau),
                     rho_sel_devp = rho ^ ((1 - tau) / tau)
                     )
  invisible(list2env(params_dep, envir=environment()))

  sigma_xi <- sqrt(sigma_xi2)
  omegaz <- sqrt(omegaz2)
  t2 <- noisy.evo.demo::Timescale_2(var_a, sigma_xi, delta, Vb)
  times <- 1:ceiling(t2)
  pred_median <- sapply(times, function(t) noisy.evo.demo::qNgamma(t, 0.5, rho, tau, omegaz^2, var_a, Vb, delta, Ve, B, alpha,
                                                   interp_alpha = FALSE, exact = FALSE, .R0=R0, .Npop0=Npop0,
                                                   .sigma_xi = sigma_xi))
  pred_Nt <- Nt(times, env=delta, omegaz=omegaz, rho=rho, alpha=alpha, R0=R0, N0=Npop0, PPpred = FALSE,
               .var_a=var_a, .Vb=Vb, .Ve=Ve, .tau=tau, .sigma_xi=sigma_xi, .B=B)
  expect_equal(pred_Nt, pred_median, tol=0.02)
  expect_equal_to_reference(pred_Nt, 'Nt.rds')

  pred_median <- sapply(times, function(t) noisy.evo.demo::qNgamma(t, 0.5, rho, tau, omegaz^2, var_a, Vb, delta, Ve, B, alpha,
                                                   interp_alpha = FALSE, exact = FALSE, .R0=R0, .Npop0=Npop0,
                                                   .sigma_xi = sigma_xi, predict_alpha = TRUE))
  pred_Nt <- Nt(times, env=delta, omegaz=omegaz, rho=rho, alpha=alpha, R0=R0, N0=Npop0, PPpred = TRUE,
               .var_a=var_a, .Vb=Vb, .Ve=Ve, .tau=tau, .sigma_xi=sigma_xi, .B=B)
  expect_equal(pred_Nt, pred_median, tol=0.02)
  expect_equal_to_reference(pred_median, 'gamma-median-predict.rds', 0.02)
  expect_equal_to_reference(pred_Nt, 'Nt-predict.rds')
})
