
#' @export
Phi <- function(Gaa, Gbb, delta)
    Gbb * delta ^ 2 / (Gaa + Gbb * delta ^ 2)

## porting from .nb file
# setting up for stoch load under plasticity after Michel Chevin et \
#al 2014. Note that delayed evaluation, :=, is necessary to substitute \
#alpha for predicted alpha in these super nested expressions

#' @export
Sigma_psi <- function(B, sigma_xi, alpha, rho_dev_sel)
    sqrt(  B ^ 2 * sigma_xi^2 * (1 + alpha * (alpha - 2 * rho_dev_sel)))

# (* stoch load under plasticity in AC env after Michel Chevin et al \2014
#' @export
Rho_psi <- function(rho_dev_sel, rho_sel, rho_dev, rho_dev_selp, rho_sel_devp, alpha)
    ## TODO --errora unused sel_dev?
    (rho_sel + alpha ^ 2 * rho_dev - alpha *
     (rho_sel_devp + rho_dev_selp)) / (1 + alpha * (alpha - 2 * rho_dev_sel))

## autocorrelated environment
#' @export
Rho_psiAC <- function(rho_dev_sel, alpha, tau)
    rho_dev_sel ^ (1/tau) *
    (1 - (alpha * (rho_dev_sel ^ (-1) -  rho_dev_sel) ) /
                      (1 + alpha * (alpha - 2 * rho_dev_sel)))

# sig_psi already has the B ^ 2 factor in it,
# accounting for \ difference btween these and stochload below

##  also using the white \ noise approximation in cases where alpha >
## [Rho] then the formula for Subscript[\[Rho], \[Psi]AC] implies
## negative autocorrelation at generational scale. This is necessary
## because the current derivation assumes this autocorrelation is
## positive (technically relying on the log of this quantity).
##  1  <= alpha ( \[Rho] ^ -1 -  \[Rho] )) /
##                (1 + alpha (alpha - 2\[Rho])) <->  \[Rho] < alpha

#' @export
Stochload_LS <- function(rho_dev_sel,  tau, B, sigma_xi, alpha, gamma_sh, var_a){
 	gamma_sh/2 * sigma_xi ^ 2 / (gamma_sh * var_a / abs(log(rho_dev_sel)) +  1)
}

#' @export
StochloadMC <- function(rho_dev_sel, alpha, gammash, var_a, rho_psi, sigma_psi) {
    ## if we can't compute it,
    ## return the whitenoise approximate load, otherwise
    ## return the load corrected for autocorrelation
    ## (agrees with lande/shannon in both cases)
    retval <- ifelse (rho_dev_sel == 0 |
                      alpha * (rho_dev_sel ^ (-1) - rho_dev_sel) / (1 + alpha * (alpha - 2 * rho_dev_sel)) >= 1,
                      gammash / 2 * sigma_psi ^ 2 * ( var_a * gammash / 2 + 1),
                      gammash / 2 * sigma_psi ^ 2  / (  gammash  * var_a / abs(log(rho_psi)) + 1)
                      )
    return(retval)
}

#' @export
Malad <- function(alpha, B, delta)
  (1 - alpha) ^ 2 * B ^ 2 * delta ^ 2

#' @export
Wbarfn <- function(t, alpha, B, gammash, delta, rmax, Gaa, Gbb,
                   tau, rho_dev_sel,
                   sigma_xi, PPpred=FALSE, phase2=FALSE)  {
  phi <-  Phi(Gaa, Gbb, delta)
  lrpgrmod <- LRPgr_mod(alpha, B, gammash, delta, rmax, Gaa, Gbb,
                        tau, rho_dev_sel, sigma_xi, PPpred)

  if (phase2) {
    abpred <- ab_pred(t, alpha, B, gammash, delta, Gaa, Gbb, rho_dev_sel, sigma_xi)
    a <- abpred[["a"]]
    b <- abpred[["b"]]
    lag_load <- gammash / 2 * (a - 0 + delta * (b - B)) ^ 2
    lrpgrmod <- LRPgr_mod(b / B, B, gammash, delta, rmax, Gaa, Gbb, tau, rho_dev_sel, sigma_xi, PPpred)
  } else {
    malad <- Malad(alpha, B, delta)
    lag_load <- gammash /2 * malad  * exp(- 2 * gammash * Gaa / (1 - phi) * t) 
  }
  rmax + lrpgrmod - lag_load
}

#' @export
Zpred <- function(t, alpha, B, gammash, delta, Gaa, Gbb) {
      phi <-  Phi(Gaa, Gbb, delta)
      B * delta * (1 - (1 - alpha) * exp( -gammash * Gaa / (1 - phi) * t))
    }

#' @export
ab_pred <- function(t, alpha, B, gammash, delta, Gaa, Gbb, rho_dev_sel, sigma_xi) {
    phi <-  Phi(Gaa, Gbb, delta)
    # to leading order in sigma_xi / delta
    eigenvec1 <- c( delta * (1 - phi), phi)
    eigenvec2 <- c( delta, -1)
    x_inf <- c(0 + B * delta * (1 - rho_dev_sel), B * rho_dev_sel)
    constants <- c(-B * (1 - alpha), -B * (alpha * (1-phi) - rho_dev_sel + phi))
    eigenvals <- c(Gaa / (1 - phi), Gaa * sigma_xi^2 / delta^2 * phi)
    pred_aorb <- function(i){
      x_inf[i] + constants[1] * eigenvec1[i] * exp( - gammash * eigenvals[1] * t) +
        constants[2] *  eigenvec2[i] * exp( - gammash * eigenvals[2] * t)
    }
    list(a = pred_aorb(1), b = pred_aorb(2))
}

#' @export
LRPgr_mod <- function(alpha, B, gammash, delta, rmax, Gaa, Gbb,
                      tau, rho_dev_sel,
                      sigma_xi, PPpred, LS=FALSE) {
    phi <-  Phi(Gaa, Gbb, delta)
    if (PPpred) {
        alphamax <- alpha   +  phi * (1 - alpha)
        lrpgrmod <- - StochloadMC(rho_dev_sel, alphamax, gammash, Gaa,
                                  Rho_psiAC(rho_dev_sel, alphamax, tau),
                                  Sigma_psi(B, sigma_xi, alphamax, rho_dev_sel))
	if(LS)
		return(- Stochload_LS(rho_dev_sel, tau, B, sigma_xi, alphamax, gammash, Gaa))
    }
    else
        lrpgrmod <- - StochloadMC(rho_dev_sel, alpha, gammash, Gaa,
                                  Rho_psiAC(rho_dev_sel, alpha, tau),
                                  Sigma_psi(B, sigma_xi, alpha, rho_dev_sel))
    lrpgrmod
}

#' @export
LRPgr <- function(alpha, B, gammash, delta, rmax, Gaa, Gbb,
                  tau, rho_dev_sel, sigma_xi, PPpred=FALSE, LS=FALSE)  {
    lrpgrmod <- LRPgr_mod(alpha, B, gammash, delta, rmax, Gaa, Gbb,
                          tau, rho_dev_sel, sigma_xi, PPpred, LS)
    return (rmax + lrpgrmod)
}


#' @export
Pgr <- function(t, alpha, B, gammash, delta, rmax, Gaa, Gbb,
                tau, rho_dev_sel,
                sigma_xi, PPpred=FALSE, LS=FALSE)  {
    malad <- Malad(alpha, B, delta)
    phi <-  Phi(Gaa, Gbb, delta)
    lrpgr <- LRPgr(alpha, B, gammash, delta, rmax, Gaa, Gbb,
                   tau, rho_dev_sel, sigma_xi, PPpred, LS)
    return (lrpgr * t - malad * (1 - phi) / (4 * Gaa) * (1 - exp(-2 * gammash * Gaa / (1 - phi) * t)))
}

#' @export
Timescale_2 <- function(var_a, sigma_xi, delta, Vb) {
  eval2 <- var_a * sigma_xi / delta ^ 2 * Phi(var_a, Vb, delta)
  1 / eval2
}


#' @export
Nt <- function(t, env, omegaz, rho, alpha, R0=omega_Wmax, N0=Npop0, PPpred=FALSE,
               .var_a=var_a, .Vb=Vb, .Ve=Ve, .tau=tau, .sigma_xi=sigma_xi, LS=FALSE, .B=B)
    log(N0) + Pgr(t, alpha, .B, gammash=1 / (omegaz ^ 2 +  Vz_(.var_a, .Vb, env, .Ve)),
                  delta=env,
                  rmax=log(R0) - 1 / 2 * log(1 + Vz_(.var_a, .Vb, env, .Ve) / (omegaz ^ 2)),
                  Gaa=.var_a, Gbb=.Vb, tau=.tau,
                  rho_dev_sel=rho, sigma_xi = .sigma_xi,
                  PPpred=PPpred, LS=LS)

#' @export
LogNTrN0 <- function(alpha, B, gammash, delta, rmax,Gaa, Gbb,
                     tau, rho_dev_sel,
                     sigma_xi,
                     PPpred=FALSE)  {
   malad <- Malad(alpha, B, delta)
   phi <-  Phi(Gaa, Gbb, delta)
   lrpgr <- LRPgr(alpha, B, gammash, delta, rmax, Gaa, Gbb,
                  tau, rho_dev_sel,
                  sigma_xi, PPpred)
   undef_res <- ifelse(lrpgr > 0, lrpgr, -Inf) # return value for else

   return(ifelse( malad / ( 2 * lrpgr) >= 1,
                     (1 - phi) / (4 * Gaa * gammash)  * (2 * lrpgr * (1 + log(malad / (2 * lrpgr ))) - malad),
                     undef_res))
   ## TODO CHECK WITH ABOVE RE gammash
}
