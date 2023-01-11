## This script can be used to simulate results from a random effects mode with
## selection publication bias, compute several CIs for the effect estimates and
## assess the CIs using several metrics.
##
## Here we investigate the effect of inbalanced study sizes. To this end, we
## we make one or two studies ten times the size of the others.
##
## It provides the functions:
## - simREbias():   simulates data from a random effects model with selection
##                  bias
##                  for meta analyses, where the effects (theta) describe mean
##                  differences.
## - simRE():       simulates data from a random effects model for meta
##                  analyses, where the effects (theta) describe mean
##                  differences.
## - pAccept():     computes the probability of publishing a study under the
##                  assumption of a 'moderate' and 'strong' publication bias
##                  as mentioned in Henmi & Copas, 2009.
## - sim2CIs():     takes data from simRE() as input and computes several
##                  confidence intervals for the combined effect estimate
## - CI2measures(): takes output from simRE() and sim2CI() and computes several
##                  quality measures of the CIs (width, coverage, score)
## - sim():         run entire simulation: generate data, compute CIs,
##                  assess CIs
##
## Florian Gerber, florian.gerber@uzh.ch, Oct. 14, 2021
rm(list = ls())
library(meta)
devtools::install_github("felix-hof/hMean")
library(hMean)
library(tidyverse)
library(rlang)
library(doParallel)
library(doRNG)
library(RhpcBLASctl)
blas_set_num_threads(1) # multi threading of BLAS
library(tictoc)
library(metafor)
library(sn)

#' Simulate effect estimates and their standard errors using a random effects
#' model
#'
#' Simulate effect estimates and their standard error using a random effects
#' model.
#' @param k number of trials
#' @param sampleSize sample size of the trial
#' @param effect effect size
#' @param I2 Higgin's I^2 heterogeneity measure
#' @param heterogeneity The heterogeneity model, the studies are simulated from.
#' Either "additive" or "multiplicative".
#' @param dist distribution to simulate the study effect. Either "t" or
#' "Gaussian".
#' the sample size as specified by \code{sampleSize}.
#' @return a matrix \code{k} x 2 matrix with columns
#' \code{theta} (effect estimates) and
#' \code{se} (standard errors).
simRE <- function(k, sampleSize, effect, I2,
                  heterogeneity,
                  dist,
                  large) {

    # get args
    n <- rep(sampleSize, k)
    # include large studies
    if (large != 0) n[seq_len(large)] <- n[seq_len(large)] * 10
    # stuff for additive model
    if (heterogeneity == "additive") {
        eps2 <- 1 / k * sum(2 / n)
        tau2 <- eps2 * I2 / (1 - I2)
        if (dist == "t") {
            ## the sn::rst(xi=0, omega, nu) distribution has variance
            ## omega^2 nu/(nu-2) (if nu>2)
            ## where nu is the degrees of freedom (dof).
            ## So if we want the variance to be tau^2, then
            ## omega^2 = tau^2 * (nu-2)/nu
            ## We use nu=4 dof then omega^2 = tau^2/2, so half as
            ## large as the heterogeneity variance under normality.
            delta <- rst(n = k, xi = effect, omega = sqrt(tau2 / 2), nu = 4)
        } else {
            delta <- rnorm(n = k, mean = effect, sd = sqrt(tau2))
        }
        theta <- rnorm(n = k, mean = delta, sd = sqrt(2 / n))

    } else { ## multiplicative model
        phi <- 1 / (1 - I2)
        eps2 <- 1 / k * sum(2 / n)
        tau2 <- eps2 * (phi - 1)
        if (dist == "t") {
            ## the sn::rst(xi=0, omega, nu) distribution has variance
            ## omega^2 nu/(nu-2) (if nu>2)
            ## where nu is the degrees of freedom (dof).
            ## So if we want the variance to be tau^2, then
            ## omega^2 = tau^2 * (nu-2)/nu
            ## We use nu=4 dof then omega^2 = tau^2/2, so half as
            ## large as the heterogeneity variance under normality.
            ## sample sequentially with marginal variance equal to
            ## (phi-1)*2/n + 2/n = phi*2/n
            delta <- rst(n = k, xi = effect, omega = sqrt(tau2 / 2), nu = 4)
        } else {  ## Gaussian, sample directly from marginal
            delta <- rnorm(n = k, mean = effect, sd = sqrt(tau2))
        }
        theta <- rnorm(n = k, mean = delta, sd = sqrt(2 / n))
    }
    se <- sqrt(rchisq(n = k, df = 2 * n - 2) / (n * (n - 1)))
    o <- cbind("theta" = theta, "se" = se, "delta" = delta)
    rownames(o) <- NULL
    return(o)
}


#' computes the probability of publishing a study under the assumption
#' of 'moderate' and 'strong' publication bias as mentioned in
#'
#' Henmi & Copas, Confidence intervals for random effects
#' meta-analysis and robustness to publication bias, 2009, eq. 23
#'
#' @param theta vector of study effects
#' @param se vector of study effects standard errors
#' @param bias either 'strong' or 'moderate'.
#' Indicating the amount of publication bias.
#' @return probabilities of publishing the studies
#' @examples
#' pAccept(theta  = c(0, 0, 1, 1, 2, 2),
#'         sigma2 = c(1, 2, 1, 2, 1, 2), bias = "moderate")
#' pAccept(theta  = c(0, 0, 1, 1, 2, 2),
#'         sigma2 = c(1, 2, 1, 2, 1, 2), bias = "strong")
pAccept <- function(theta, se, bias) {
    ## Begg & Mazumdar, Biometrics, 1994
    ## moderate bias: beta = 4, gamma = 3
    ## strong bias:   beta = 4, gamma = 1.5

    if (bias == "moderate") {
        beta <- 4
        gamma <- 3
    } else {
        beta <- 4
        gamma <- 1.5
    }

    exp(-beta * (dnorm(-theta / se))^gamma)
}

#' Simulate effect estimates and their standard errors using a random effects
#' model under none, moderate, or strong publication bias
#'
#' @param k number of trials
#' @param sampleSize sample size of the trial
#' @param effect effect size
#' @param I2 Higgin's I^2 heterogeneity measure
#' @param heterogeneity The heterogeneity model, the studies are simulated from.
#' Either "additive" or "multiplicative".
#' @param dist distribution to simulate the study effect. Either "t" or
#' "Gaussian".
#' @param large A number in \code{c(0,1,2)} indicating the number of studies
#' that have ten times the sample size as specified by \code{sampleSize}.
#' Publication bias is only applied to the smaller studies with sample size
#' specified by \code{sampleSize}.
#' @param bias either 'none', 'moderate' or 'strong' as used in
#' Henmi & Copas (2010).
#' @references
#' Henmi, M. and Copas, J. B. (2010). Confidence intervals for random effects
#' meta-analysis and robustness to publication bias. Statistics in Medicine,
#' 29(29):2969-2983.
#' @return a matrix \code{k} x 2 matrix with columns
#' \code{theta} (effect estimates) and
#' \code{se} (standard errors).
#' @examples
#' simREbias(4, sampleSize= 50, effect=.2, I2=.3, dist="Gaussian", large=0)
#' simREbias(4, sampleSize= 50, effect=.2, I2=.3, dist="Gaussian", large=1)
#' simREbias(4, sampleSize= 50, effect=.2, I2=.3, dist="Gaussian", large=2)
#' simREbias(4, sampleSize= 50, effect=.2, I2=.3, dist="Gaussian", large=2,
#'           bias = "moderate")
#' simREbias(4, sampleSize= 50, effect=.2, I2=.3, dist="Gaussian",
#'           large=1, bias = "strong")
simREbias <- function(k, sampleSize, effect, I2,
                      heterogeneity = c("additive", "multiplicative"),
                      dist = c("t", "Gaussian"),
                      large,
                      bias = c("none", "moderate", "strong"),
                      verbose = TRUE,
                      check_inputs = TRUE) {
    # input checks
    if (check_inputs) {
        stopifnot(length(k) == 1L,
                  is.numeric(k),
                  is.numeric(sampleSize),
                  is.finite(sampleSize),
                  length(sampleSize) == 1L,
                  is.numeric(effect),
                  length(effect) == 1L,
                  is.numeric(effect),
                  is.finite(effect),
                  length(I2) == 1,
                  is.numeric(I2),
                  0 <= I2, I2 < 1,
                  is.character(heterogeneity),
                  length(heterogeneity) == 1L,
                  !is.na(heterogeneity),
                  is.character(dist),
                  length(dist) == 1L,
                  !is.na(dist),
                  is.numeric(large),
                  length(large) == 1L,
                  is.finite(large),
                  large %in% c(0, 1, 2),
                  is.character(bias),
                  length(bias) == 1L,
                  !is.na(bias),
                  k >= large)
    }

    bias <- match.arg(bias)
    dist <- match.arg(dist)
    heterogeneity <- match.arg(heterogeneity)

    if (bias == "none") {
        o <- simRE(k = k, sampleSize = sampleSize, effect = effect, I2 = I2,
            heterogeneity = heterogeneity, dist = dist, large = large)
        ## add attributes and return
        attr(o, "heterogeneity") <- heterogeneity
        attr(o, which = "effect") <- effect
        return(o)
    }

    ## first ignore the 'large'
    o <- simRE(k = k * 3, sampleSize = sampleSize, effect = effect, I2 = I2,
        heterogeneity = heterogeneity, dist = dist, large = 0)
    pa <- pAccept(theta = o[, "theta"], se = o[, "se"], bias = bias)
    keep <- rbernoulli(n = k * 3, p = pa)
    while (k > sum(keep)) {
        if (verbose) cat(".")
        o2 <- simRE(k = k * 3, sampleSize = sampleSize, effect = effect,
            I2 = I2, heterogeneity = heterogeneity, dist = dist,  large = 0)
        pa2 <- pAccept(theta = o2[, "theta"], se = o2[, "se"], bias = bias)
        keep2 <- rbinom(n = k * 3, size = 1, p = pa2)
        #keep2 <- rbernoulli(n = k * 3, p = pa2)
        o <- rbind(o, o2)
        keep <- c(keep, keep2)
    }
    o <- o[keep, ][1:k, ]

    ## add large studies
    if(large != 0){
        oLarge <- simRE(k = large, sampleSize = sampleSize, effect = effect, I2 = I2, heterogeneity = heterogeneity, dist = dist, large = large)
        o <- rbind(oLarge, o[-seq_len(large),]) 
    }
    
    ## add attributes and return
    attr(o, which = "heterogeneity") <- heterogeneity
    attr(o, which = "effect") <- effect
    o
}


#' Confidence intervals from effect estimates and their standard errors
#'
#' Takes the output of \code{simRE} and returns CIs for the combined effect using the
#' indicated methods.
#' @param x matrix output from \code{simRE}.
#' @return a tibble with columns \code{lower}, \code{upper}, and \code{method}.
sim2CIs <- function(x){
    ## Henmi & Copas confidence Interval
    HC <- metafor::hc(object = metafor::rma(yi = x[, "theta"], sei = x[, "se"], 
                                            control = list(maxiter = 10000, stepadj = 0.25)))
    
    ## standard metagen with REML estimation of tau
    REML <- meta::metagen(TE = x[, "theta"], seTE = x[, "se"], sm = "MD", 
                          method.tau = "REML", 
                          control = list(maxiter = 10000, stepadj = 0.25))
    
    ## Hartung & Knapp
    HK <- meta::metagen(TE = x[, "theta"], seTE = x[, "se"], sm = "MD", 
                        method.tau = "REML", hakn = TRUE,
                        control = list(maxiter = 10000, stepadj = 0.25))
    
    ## HMean2sided
    # if(nrow(x) <= 5) {
    #     HM2_f <- hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"],
    #                           alternative = "two.sided", distr = "f",
    #                           heterogeneity = "additive")
    #     HM2_chisq <- hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"],
    #                               alternative = "two.sided", distr = "chisq",
    #                               heterogeneity = "additive")
    # } else {
    #     HM2_f <- list(CI = cbind(lower = NA_real_, upper = NA_real_))
    #     HM2_chisq <- list(CI = cbind(lower = NA_real_, upper = NA_real_))
    # }
    
    ## Note: these here are actually not really necessary anymore. In an earlier version, the
    ## argument tau2 defaulted to 0 and the heterogeneity was always additive.
    
    ## HMeanNone (additive model with tau2 = 0)
    HM_f <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"],
                         alternative = "none", distr = "f", tau2 = 0,
                         heterogeneity = "additive", check_inputs = TRUE)
    HM_chisq <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"],
                             alternative = "none", distr = "chisq", tau2 = 0,
                             heterogeneity = "additive", check_inputs = TRUE)
    
    ## HMeanNone_tau2 (additive with estimated tau2)
    HM_tau2_f <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                              tau2 = REML$tau2, alternative = "none", distr = "f",
                              heterogeneity = "additive", check_inputs = TRUE)
    HM_tau2_chisq <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                                  tau2 = REML$tau2, alternative = "none", distr = "chisq",
                                  heterogeneity = "additive", check_inputs = TRUE)
    
    ## HMeanNone_phi (multiplicative with estimated phi)
    phi <- hMean::estimatePhi(thetahat = x[, "theta"], se = x[, "se"])
    HM_phi_f <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                             phi = phi, alternative = "none", distr = "f", 
                             heterogeneity = "multiplicative", check_inputs = TRUE)
    HM_phi_chisq <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                                 phi = phi, alternative = "none", distr = "chisq", 
                                 heterogeneity = "multiplicative", check_inputs = TRUE)
    
    ## HMean_None with k-trials (multiplicative with estimated phi)
    HM_ktrial_none <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                                          alternative = "none",
                                          heterogeneity = "none", check_inputs = TRUE,
                                          pValueFUN = hMean::kTRMu)
    HM_ktrial_mult <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                                          phi = phi, alternative = "none",
                                          heterogeneity = "multiplicative", check_inputs = TRUE,
                                          pValueFUN = hMean::kTRMu)
    HM_ktrial_add <- hMean::hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                                         tau2 = REML$tau2, alternative = "none",
                                         heterogeneity = "additive", check_inputs = TRUE,
                                         pValueFUN = hMean::kTRMu)
    
    tib <- tibble(lower = c(HC$ci.lb,
                            REML$lower.random,
                            REML$lower.predict,
                            HK$lower.random,
                            HK$lower.predict,
                            #HM2_chisq$CI[,"lower"],
                            #HM2_f$CI[,"lower"],
                            HM_chisq$CI[,"lower"],
                            HM_f$CI[,"lower"],
                            HM_tau2_chisq$CI[,"lower"],
                            HM_tau2_f$CI[,"lower"],
                            HM_phi_chisq$CI[,"lower"],
                            HM_phi_f$CI[,"lower"],
                            HM_ktrial_none$CI[,"lower"],
                            HM_ktrial_add$CI[,"lower"],
                            HM_ktrial_mult$CI[,"lower"]),
                  upper = c(HC$ci.ub,
                            REML$upper.random,
                            REML$upper.predict,
                            HK$upper.random,
                            HK$upper.predict,
                            #HM2_chisq$CI[,"upper"],
                            #HM2_f$CI[,"upper"],
                            HM_chisq$CI[,"upper"],
                            HM_f$CI[,"upper"],
                            HM_tau2_chisq$CI[,"upper"],
                            HM_tau2_f$CI[,"upper"],
                            HM_phi_chisq$CI[,"upper"],
                            HM_phi_f$CI[,"upper"],
                            HM_ktrial_none$CI[,"upper"],
                            HM_ktrial_add$CI[,"upper"],
                            HM_ktrial_mult$CI[,"upper"]),
                  method = c("Henmy & Copas CI",
                             "REML CI",
                             "REML PI", 
                             "Hartung & Knapp CI",
                             "Hartung & Knapp PI",
                             #"Harmonic Mean two sided CI (chisq)",
                             #"Harmonic Mean two sided CI (f)",
                             rep("Harmonic Mean CI (chisq)", nrow(HM_chisq$CI)),
                             rep("Harmonic Mean CI (f)", nrow(HM_f$CI)),
                             rep("Harmonic Mean Additive CI (chisq)", nrow(HM_tau2_chisq$CI)),
                             rep("Harmonic Mean Additive CI (f)", nrow(HM_tau2_f$CI)),
                             rep("Harmonic Mean Multiplicative CI (chisq)", nrow(HM_phi_chisq$CI)),
                             rep("Harmonic Mean Multiplicative CI (f)", nrow(HM_phi_f$CI)),
                             rep("k-Trials CI", nrow(HM_ktrial_none$CI)),
                             rep("k-Trials Additive CI", nrow(HM_ktrial_add$CI)),
                             rep("k-Trials Multiplicative CI", nrow(HM_ktrial_mult$CI))
                             ))
    out <- list(CIs = tib,
                model = attributes(x)$heterogeneity,
                gamma = tibble("method" = c("Harmonic Mean CI (chisq)", "Harmonic Mean CI (f)", 
                                            "Harmonic Mean Additive CI (chisq)", "Harmonic Mean Additive CI (f)", 
                                            "Harmonic Mean Multiplicative CI (chisq)", "Harmonic Mean Multiplicative CI (f)",
                                            "k-Trials CI", "k-Trials Additive CI", "k-Trials Multiplicative CI"),
                               "gamma_min" = c(min(HM_chisq$gamma[,2]), min(HM_f$gamma[,2]), 
                                               min(HM_tau2_chisq$gamma[,2]), min(HM_tau2_f$gamma[,2]),
                                               min(HM_phi_chisq$gamma[,2]), min(HM_phi_f$gamma[,2]),
                                               min(HM_ktrial_none$gamma[,2]), min(HM_ktrial_add$gamma[,2]), min(HM_ktrial_mult$gamma[,2])),
                               "x_gamma_min" = c(HM_chisq$gamma[which.min(HM_chisq$gamma[,2]),1], HM_f$gamma[which.min(HM_f$gamma[,2]),1],
                                                 HM_tau2_chisq$gamma[which.min(HM_tau2_chisq$gamma[,2]),1], HM_tau2_f$gamma[which.min(HM_tau2_f$gamma[,2]),1],
                                                 HM_phi_chisq$gamma[which.min(HM_phi_chisq$gamma[,2]),1], HM_phi_f$gamma[which.min(HM_phi_f$gamma[,2]),1],
                                                 HM_ktrial_none$gamma[which.min(HM_ktrial_none$gamma[,2]),1],
                                                 HM_ktrial_add$gamma[which.min(HM_ktrial_add$gamma[,2]),1], HM_ktrial_mult$gamma[which.min(HM_ktrial_mult$gamma[,2]),1])),
                theta = x[, "theta"],
                delta = x[, "delta"],
                effect = attributes(x)$effect)
    out
}



#' Computes quality measures for CIs
#'
#' @param x a list with elements \code{CIs}, \code{model}, \code{gamma}, \code{theta}, 
#' \code{delta} and \code{effect} as obtained from \code{sim2CIs}.
#' @param pars a \code{data.frame} with one row. It should have column names \code{k}, \code{sampleSize}, \code{effect}, 
#' \code{I2}, \code{heterogeneity}, \code{dist} and \code{large}. These are passed as parameters to
#' \link{simREbias()}. In order to simulate another study in order to assess prediction interval coverage.
#' @return a tibble with columns
#' \item{\code{method}}{method}
#' \item{\code{width}}{with of the intervals}
#' \item{\code{coverage}}{covarage of the true value 0}
#' \item{\code{score}}{interval score as defined in Gneiting and Raftery (2007)}
#' \item{\code{coverage_effects}}{Proportion of study effects covered by the interval(s).}
#' \item{\code{n}}{Number of intervals}
CI2measures <- function(x, pars) {
    
    methods <- unique(x$CIs$method)
    
    foreach(i = seq_along(methods), .combine = rbind) %do% {
        
        # Subset by method
        x_sub <- x$CIs %>% filter(method == methods[i]) %>%
            select(lower, upper) %>%
            as.matrix()
        
        # calculate proportion of deltas covered by the interval
        coverage_effects <- vapply(x$delta, function(delta){
            any(x_sub[,"lower"] <= delta & delta <= x_sub[,"upper"])
        }, logical(1L)) %>%
            mean()
        
        # calculate how many times at least one study-specific effect is covered
        found <- FALSE
        for(z in x$delta){
          if(any(x_sub[,"lower"] <= z & z <= x_sub[,"upper"])){
            found <- TRUE
            break
          }
        }
        coverage_effects_min1 <- as.numeric(found)
        
        # calculate whether all deltas covered by CI
        coverage_effects_all <- as.numeric(all(vapply(x$delta, function(delta) any(x_sub[,"lower"] <= delta & delta <= x_sub[,"upper"]), logical(1L))))
        
        # calculate whether interval covers future study (only for hMean, hMean_additive, hMean_multiplicative)
        new_study <- simREbias(k = 1, sampleSize = pars$sampleSize, effect = pars$effect, 
                               I2 = pars$I2, heterogeneity = pars$heterogeneity, 
                               dist = pars$dist, large = 0, bias = "none")[, "delta"]
        coverage_prediction <- as.numeric(any(x_sub[,"lower"] <= new_study & new_study <= x_sub[,"upper"]))
        
        # get gamma_min to later attach it to the function output
        if(grepl("Harmonic Mean.*CI|k-Trials.*CI", methods[i])){
            gamma_min <- x$gamma %>% filter(method == methods[i]) %>% pull(gamma_min)
        } else {
            gamma_min <- NA_real_
        }
        
        # calculate total width of the interval(s)
        width <- {x_sub[,"upper"] - x_sub[,"lower"]} %>% sum(.)
        
        # calculate whether interval covers true effect
        coverage_true <- as.numeric(any(x_sub[,"lower"] <= x$effect & x$effect <= x_sub[,"upper"]))
        
        # Calculate score for CI methods
        if(grepl("CI", methods[i])){
            score <- width + (2/0.05) * min(abs(x_sub[,"lower"]), abs(x_sub[,"upper"])) * (1 - coverage_true)
        } else {
            score <- NA_real_
        }
        
        # count number of intervals
        if(grepl("Harmonic Mean.*CI|k-Trials.*CI", methods[i])){
            n <- nrow(x_sub)
        } else {
            n <- NA_real_
        }
        
       out <- tibble(method = methods[i], 
                     coverage_true = coverage_true, 
                     coverage_effects = coverage_effects,
                     coverage_effects_min1 = coverage_effects_min1,
                     coverage_effects_all = coverage_effects_all,
                     coverage_prediction = coverage_prediction,
                     gammaMin = gamma_min, 
                     n = n,
                     width = width,
                     score = score) 
        
        # return
        out
    }
}

error_function <- function(cond, pars, error_obj = NULL, fun_name){
    text <- capture.output(cond)
    out_msg <- paste0("Error in ", fun_name , "iteration: ", i, "\n\n", "Parameters are:\n\n",
                      paste0(paste0(names(pars), ": ", pars[1, ]), collapse = "\n"),
                      "\n\nThe error message is:\n", text, "\n\n",
                      ifelse(is.null(error_obj), 
                             "\nSee pars.rds for the saved parameters.\n\n", 
                             "\nSee error.rds for the last available object and pars.rds for the saved parameters.\n\n"),
                      "--------------------------------------------------------------------------------------\n")
    cat(out_msg, file = "error.txt", append = TRUE)
    saveRDS(pars, file = "pars.rds")
    if(!is.null(error_obj)) saveRDS(error_obj, file = "error.rds")
    return(NA)
}


#' Simulate N times, compute CIs, assess CIs
#'
#' Takes a data.frame of parameter configurations, the number replicates \code{N} and
#' the number of cores to simulate accordingly.
#'
#' @param grid data.frame with columns
#' \item{sampleSize}{Sample size of the trial}
#' \item{effect}{mean effect of simulated trials}
#' \item{I2}{Higgin's I^2 heterogeneity measure}
#' \item{k}{number of studies}
#' \item{dist}{distribution used to simulate the study effects from}
#' @param N number of replications to be done for each scenario.
#' @param cores number of CPUs to use for the simulation.
#' The default value is obtained from \code{detectCores()}. 
#' @return tibble in long format with columns
#' \item{sampleSize}{sample size}
#' \item{I2}{Higgin's I^2 heterogeneity measure}
#' \item{k}{number of studies}
#' \item{method}{CI method}
#' \item{measure}{measure to assess the CI}
#' \item{value}{value of the measure}
sim <- function(grid, N = 1e4, cores = detectCores(), seed = as.numeric(Sys.time())){
    
    types <- vapply(grid, typeof, character(1L), USE.NAMES = TRUE)
    
    # check grid columns
    stopifnot(is.data.frame(grid),
              c("sampleSize", "I2", "k", "dist", 
                "effect", "large", "heterogeneity", "bias") %in% names(grid),
              types["sampleSize"] %in% c("double", "numeric", "integer"),
              types["effect"] %in% c("double", "numeric", "integer"),
              types["I2"] %in% c("double", "numeric", "integer"),
              types["k"] %in% c("double", "numeric", "integer"),
              types["heterogeneity"] %in% c("character"),
              types["dist"] %in% c("character"),
              types["bias"] %in% c("character"),
              types["large"] %in% c("double", "numeric", "integer"))
    
    # check grid
    stopifnot(is.numeric(grid$k),
              is.numeric(grid$sampleSize),
              all(is.finite(grid$sampleSize)),
              is.numeric(grid$effect),
              is.numeric(grid$effect),
              all(is.finite(grid$effect)),
              is.numeric(grid$I2),
              all(0 <= grid$I2), all(grid$I2 < 1),
              is.character(grid$heterogeneity),
              all(!is.na(grid$heterogeneity)),
              is.character(grid$dist),
              all(!is.na(grid$dist)),
              is.numeric(grid$large),
              all(is.finite(grid$large)),
              all(grid$large %in% c(0,1,2)),
              is.character(grid$bias),
              all(!is.na(grid$bias)),
              
              all(grid$k >= grid$large))
    
    # check other arguments
    stopifnot(is.numeric(N), length(N) == 1L, 1 <= N,
              is.numeric(seed), length(seed) == 1L)
    
    # register parallel backend
    registerDoParallel(cores)
    on.exit(stopImplicitCluster())
    
    # run simulation
    foreach(j = seq_len(nrow(grid)), .options.RNG=seed, .errorhandling = "pass") %dorng% {
        
        if(file.exists("error.txt")){
            # if error happened, skip the rest of the loop iterations
            out <- "skipped"
        } else {
            
            cat("start", j, "of", nrow(grid), fill=TRUE)
            grid %>% slice(j) -> pars
            
            # av is a list with elements that are either a tibble or NA (the latter in case of an error)
            av <- foreach(i = seq_len(N), .errorhandling = "pass") %do% {
                
                # If an error happened somewhere in the simulation return NA for all successive iterations
                if(file.exists("error.txt")) return(NA)
                
                # Repeat this N times. Simulate studies, calculate CIs, calculate measures
                res <- tryCatch({
                    #if(i == 17) stop("Error in iteration 17")
                    simREbias(k = pars$k, sampleSize = pars$sampleSize,
                              effect = pars$effect, I2 = pars$I2,
                              heterogeneity = pars$heterogeneity,
                              dist = pars$dist, bias = pars$bias,
                              large = pars$large, check_inputs = FALSE)
                }, error = function(cond) error_function(cond = cond, pars = pars, fun_name = "simREbias") )
                CIs <- tryCatch({
                    if(length(res) == 1L && is.na(res)) NA else sim2CIs(x = res)
                }, error = function(cond) error_function(cond = cond, pars = pars, error_obj = res, fun_name = "sim2CIs"))
                out <- tryCatch({
                    if(length(res) == 1L && is.na(CIs)) NA else CI2measures(x = CIs, pars = pars)
                    }, error = function(cond) error_function(cond = cond, pars = pars, error_obj = res, fun_name = "CI2measures"))
                out
            }
            
            # summarize the N tibbles. If any list element is NA, return "failed".
            if(any(is.na(av))){
                out <- "failed"
            } else {
                # rbind data frames in list "av"
                av <- bind_rows(av)
                # summarize simulations
                out <- tryCatch({
                    bind_rows(
                        ## mean for all measures
                        av %>% group_by(method) %>%
                            summarize(across(everything(), 
                                             .fns = list(mean = function(x) mean(x, na.rm = FALSE)),
                                             .names = "{.col}_{.fn}"),
                                      .groups = "drop") %>% 
                            pivot_longer(cols = !method,
                                         names_to = "measure",
                                         values_to = "value") %>%
                            cbind(pars, .),
                        
                        ## summary statistics for gamma_min
                        av %>% filter(grepl("Harmonic Mean.*CI|k-Trials.*CI", method)) %>%
                            group_by(method) %>%
                            summarize(across("gammaMin", 
                                             .fns = list(min = function(x) min(x, na.rm = FALSE), 
                                                         firstQuart = function(x) quantile(x, probs = 0.25, na.rm = FALSE, names = FALSE),
                                                         median = function(x) median(x, na.rm = FALSE),
                                                         mean = function(x) mean(x, na.rm = FALSE),
                                                         thirdQuart = function(x) quantile(x, probs = 0.75, na.rm = FALSE, names = FALSE),
                                                         max = function(x) max(x, na.rm = FALSE)),
                                                         .names = "{.col}_{.fn}"),
                                      .groups = "drop") %>%
                            pivot_longer(cols = !method,
                                         names_to = "measure",
                                         values_to = "value") %>%
                            cbind(pars, .),
                        
                        ## relative frequency for n
                        av %>% filter(grepl("Harmonic Mean.*CI|k-Trials.*CI", method)) %>%
                            group_by(method) %>%
                            summarize(across("n", 
                                             .fns = list("1" = function(x) sum(x == 1),
                                                         "2" = function(x) sum(x == 2),
                                                         "3" = function(x) sum(x == 3),
                                                         "4" = function(x) sum(x == 4),
                                                         "5" = function(x) sum(x == 5),
                                                         "6" = function(x) sum(x == 6),
                                                         "7" = function(x) sum(x == 7),
                                                         "8" = function(x) sum(x == 8),
                                                         "9" = function(x) sum(x == 9),
                                                         "gt9" = function(x) sum(x >= 10)),
                                             .names = "{.col}_{.fn}"),
                                      .groups = "drop") %>%
                            pivot_longer(cols = !method,
                                         names_to = "measure",
                                         values_to = "value") %>%
                            filter(!is.na(value)) %>% 
                            cbind(pars, .)
                    )
                }, error = function(cond) error_function(cond = cond, pars = pars, error_obj = av, fun_name = "summarizing measures"))
            }
        }
        out
        
    } -> o
    attr(o, "seed") <- seed
    attr(o, "N") <- N
    o
}




## set parameter grid to be evaluated 
grid <- expand.grid(sampleSize = 50,                                 # sample size of trial
                    effect = 0.2,                                    # average effect, impacts selection bias
                    I2 = c(0, 0.3, 0.6, 0.9),                        # Higgin's I^2 heterogeneity measure
                    k = c(3, 5, 10, 20, 50),                         # number of studies
                    heterogeneity = c("additive", "multiplicative"), # The heterogeneity model that the studies are simulated from
                    dist = c("Gaussian", "t"),                       # distribution 
                    bias = c("none", "moderate", "strong"),          # bias
                    large = c(0, 1, 2),                              # number of large studies
                    stringsAsFactors = FALSE) %>% 
    arrange(desc(heterogeneity))


## run simulation, e.g., on the Rambo server of I-MATH
tic()
out <- sim(grid = grid, N = 1e4, cores = 120)
toc()

## save results
dir.create("RData", showWarnings=FALSE)
sessionInfo <- sessionInfo()
print(sessionInfo)
save(out, sessionInfo, file = "RData/simulate_all.RData")


