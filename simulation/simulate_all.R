## This script can be used to simulate results from a random effects mode with
## selection publication bias, compute several CIs for the effect estimates and
## assess the CIs using several metrics.
##
## Here we investigate the effect of inbalanced study sizes. To this end, we
## we make one or two studies ten times the size of the others.
##
## It provides the functions:
## - simREbias():   simulates data from a random effects model with selection bias
##                  for meta analyses, where the effects (theta) describe mean differences.  
## - simRE():       simulates data from a random effects model for meta analyses,
##                  where the effects (theta) describe mean differences.  
## - pAccept():     computes the probability of publishing a study under the assumption
##                  of a 'moderate' and 'strong' publication bias as mentioned in
##                  Henmi & Copas, 2009.
## - sim2CIs():     takes data from simRE() as input and computes several
##                  confidence intervals for the combined effect estimate
## - CI2measures(): takes output from simRE() and sim2CI() and computes several
##                  quality measures of the CIs (width, coverage, score)
## - sim():         run entire simulation: generate data, compute CIs, assess CIs
##
## Florian Gerber, florian.gerber@uzh.ch, Oct. 14, 2021
rm(list = ls())
library(meta)
source("ReplicationSuccess_extension.R")
library(tidyverse); theme_set(theme_bw())
library(doParallel)
library(doRNG)
library(RhpcBLASctl); blas_set_num_threads(1) # multi threading of BLAS
library(tictoc)
library(metafor)
library(sn)

#' Simulate effect estimates and their standard errors using a random effects model
#'
#' Simulate effect estimates and their standard error using a random effects model.
#' @param k number of trials
#' @param sampleSize sample size of the trial
#' @param effect effect size
#' @param I2 Higgin's I^2 heterogeneity measure
#' @param dist distribution to simulate the study effect. Either "t" or "Gaussian".
#' the sample size as specified by \code{sampleSize}.
#' @return a matrix \code{k} x 2 matrix with columns
#' \code{theta} (effect estimates) and
#' \code{se} (standard errors).
simRE <- function(k, sampleSize, effect, I2, dist=c("t", "Gaussian"),
                  large = 0) {
    stopifnot(length(k) == 1,
              length(sampleSize) == 1,
              is.numeric(effect),
              length(effect) == 1,
              is.finite(effect),
              length(I2) == 1,
              0 <= I2, I2 < 1,
              !is.null(dist),
              !is.null(large),
              is.numeric(large),
              length(large)==1,
              large %in% c(0,1,2))
    dist <- match.arg(dist)
    n <- rep(sampleSize, k)
    if(large == 1)
        n[1] <- 10 * n[1]
    if(large == 2)
        n[1:2] <- 10 * n[1:2]
    phi <- 1/(1 - I2)
    if(dist == "t") {
        ## the sn::rst(xi=0, omega, nu) distribution has variance
        ## omega^2 nu/(nu-2) (if nu>2)
        ## where nu is the degrees of freedom (dof).
        ## So if we want the variance to be tau^2, then 
        ## omega^2 = tau^2 * (nu-2)/nu
        ## We use nu=4 dof then omega^2 = tau^2/2, so half as
        ## large as the heterogeneity variance under normality. 
        theta <- rst(n = k, xi = effect, omega = sqrt(phi/n), nu = 4)
    } else {
        theta <- rnorm(n = k, mean = effect, sd = sqrt(phi*2/n))
    }
    ## theta[1:ceiling(r*k)] <- theta[1:ceiling(r*k)] + bias
    se <- sqrt(rchisq(n = k, df = 2*n - 2) / (n*(n - 1)))
    o <- cbind("theta" = theta, "se" = se)
    rownames(o) <- NULL
    o
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
pAccept <- function(theta, se, bias = c("moderate", "strong")){
    ## Begg & Mazumdar, Biometrics, 1994
    ## moderate bias: beta = 4, gamma = 3
    ## strong bias:   beta = 4, gamma = 1.5
    stopifnot(!is.null(bias))
    bias <- match.arg(bias)
    if(bias == "moderate"){
        beta <- 4
        gamma <- 3
    } else {
        beta <- 4
        gamma <- 1.5
    }
    stopifnot(length(beta) == 1, is.finite(beta), 0 < beta,
              length(gamma) == 1, is.finite(gamma), 0 < gamma)
    
    exp(-beta * (dnorm(-theta / se))^gamma)
}

#' Simulate effect estimates and their standard errors using a random effects model
#' under none, moderate, or strong publication bias
#'
#' @param k number of trials
#' @param sampleSize sample size of the trial
#' @param effect effect size
#' @param I2 Higgin's I^2 heterogeneity measure
#' @param dist distribution to simulate the study effect. Either "t" or "Gaussian".
#' @param large A number in \code{c(0,1,2)} indicating the number of studies that have ten times
#' the sample size as specified by \code{sampleSize}. Publication bias is only applied to the smaller
#' studies with sample size specified by \code{sampleSize}. 
#' @param bias either 'none', 'moderate' or 'strong' as used in Henmi & Copas (2010).
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
#' simREbias(4, sampleSize= 50, effect=.2, I2=.3, dist="Gaussian", large=2, bias = "moderate")
#' simREbias(4, sampleSize= 50, effect=.2, I2=.3, dist="Gaussian", large=1, bias = "strong")
simREbias <- function(k, sampleSize, effect, I2, dist = c("t", "Gaussian"), large = 0,
                      bias = c("none", "moderate", "strong"),
                      verbose = TRUE){
    stopifnot(!is.null(bias))
    bias <- match.arg(bias)

    if(bias == "none")
        return(simRE(k=k, sampleSize=sampleSize, effect = effect, I2=I2, dist=dist, large=large))


    ## first ignore the 'large' 
    o <- simRE(k=k*3, sampleSize=sampleSize, effect=effect, I2=I2, dist=dist)
    pa <- pAccept(theta = o[,"theta"], se = o[,"se"], bias = bias)
    keep <- rbernoulli(n = k*3, p = pa)
    while(k > sum(keep)) {
        if(verbose)
            cat(".")
        o2 <- simRE(k=k*3, sampleSize=sampleSize, effect = effect, I2=I2, dist=dist,  large=0)
        pa2 <- pAccept(theta = o2[,"theta"], se = o2[,"se"], bias = bias)
        keep2 <- rbernoulli(n = k*3, p = pa2)
        o <- rbind(o, o2)
        keep <- c(keep, keep2)
    }
    o <- o[keep,][1:k,]

    ## add large studies
    if(large == 1){
        oLarge <- simRE(k=k, sampleSize=sampleSize, effect = effect, I2=I2, dist=dist, large=large)
        o <- rbind(oLarge[1,], o[-1,]) 
    }
    if(large == 2){
        oLarge <- simRE(k=k, sampleSize=sampleSize, effect = effect, I2=I2, dist=dist, large=large)
        o <- rbind(oLarge[1:2,], o[-c(1,2),]) 
    }
    o    
}


#' Confidence intervals from effect estimates and their standard errors
#'
#' Takes the output of \code{simRE} and returns CIs for the combined effect using the
#' indicated methods.
#' @param x matrix output from \code{simRE}.
#' @return a tibble with columns \code{lower}, \code{upper}, and \code{method}.
sim2CIs <- function(x){
    ## Henmy & Copas confidence Interval
    HC <- metafor::hc(object = metafor::rma(yi = x[, "theta"], sei = x[, "se"]))

    ## standard metagen with REML estimation of tau
    REML <- metagen(TE = x[, "theta"], seTE = x[, "se"], sm = "MD", 
                    method.tau = "REML")

    ## Hartung & Knapp
    HK <- metagen(TE = x[, "theta"], seTE = x[, "se"], sm = "MD", 
                  method.tau = "REML", hakn = TRUE)

    ## HMeam2sided
    if(nrow(x) <= 5) {
        HM2 <- hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                            alternative = "two.sided")
    } else {
        HM2 <- list(CI = cbind(lower = NA, upper = NA))
    }

    ## HMeanNone
    HM <- hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                       alternative = "none")

    ## HMeanNone_tau2
    HM_tau2 <- hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                           tau2 = REML$tau2, alternative = "none")

    ## HMeanNone_phi
    HM_phi <- hMeanChiSqCIphi(thetahat = x[, "theta"], se = x[, "se"], 
                               alternative = "none")

    tibble(lower = c(HC$ci.lb,
                     REML$lower.random,
                     HK$lower.random,
                     HM2$CI[,"lower"], 
                     HM$CI[,"lower"],
                     HM_tau2$CI[,"lower"],
                     HM_phi$CI[,"lower"]),
           upper = c(HC$ci.ub,
                     REML$upper.random,
                     HK$upper.random,
                     HM2$CI[,"upper"], 
                     HM$CI[,"upper"],
                     HM_tau2$CI[,"upper"],
                     HM_phi$CI[,"upper"]),
           method = c("Henmy & Copas",
                      "REML",
                      "Hartung & Knapp",
                      "Harmonic Mean two sided",
                      rep("Harmonic Mean", nrow(HM$CI)),
                      rep("Harmonic Mean Additive", nrow(HM_tau2$CI)),
                      rep("Harmonic Mean Multiplicative", nrow(HM_phi$CI))))
}



#' Computes quality measures for CIs
#'
#' @param x a tibble with columns \code{lower}, \code{upper}, and \code{method}
#' as obtained from \code{sim2CIs}.
#' @param effect effect size.
#' @return a tibble with columns 
#' \item{\code{method}}{method}
#' \item{\code{width}}{with of the intervals}
#' \item{\code{coverage}}{covarage of the true value 0}
#' \item{\code{score}}{interval score as defined in Gneiting and Raftery (2007)}
#' \item{\code{n}}{Number of intervals}
CI2measures <- function(x, effect) {
    methods <- unique(x$method)
    foreach(i = seq_along(methods), .combine = rbind) %do% {
        x %>% filter(method == methods[i]) %>%
            select(lower, upper) %>%
            as.matrix() ->
            x_sub

        { x_sub[,"upper"] - x_sub[,"lower"] } %>%
            sum(.) ->
            width
        
        as.numeric(any(x_sub[,"lower"] <= effect & effect <= x_sub[,"upper"])) ->
            coverage
        
        width + (2/0.05) * min(abs(x_sub[,"lower"]), abs(x_sub[,"upper"])) * (1 - coverage) ->
            score
        
        tibble(method = methods[i], width = width, coverage = coverage, score = score,
               n = nrow(x_sub))
    }
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
    stopifnot(is.data.frame(grid),
              c("sampleSize", "I2", "k", "dist", "effect") %in% names(grid),
              is.numeric(N), length(N) == 1, 1 <= N)
    registerDoParallel(cores)
    foreach(j = seq_len(nrow(grid)), .combine = rbind, .options.RNG=seed) %dorng% {
        cat("start", j, "of", nrow(grid), fill=TRUE)
        grid %>% slice(j) -> pars

        foreach(i = seq_len(N), .combine = rbind) %do% {
            res <- simREbias(k = pars$k, sampleSize = pars$sampleSize,
                             effect = pars$effect, I2 = pars$I2,
                             dist = pars$dist, bias = pars$bias)
            CIs <- sim2CIs(x = res)
            CI2measures(x = CIs, effect = pars$effect)
        } -> av

        ## compute the mean values of each measure
        av %>% group_by(method) %>%
            summarize(width_mean = mean(width),
                      coverage_mean = mean(coverage),
                      score_mean = mean(score),
                      n = mean(n),
                      ## width_sd = sd(width),
                      ## coverage_sd = sd(coverage),
                      ## score_sd = sd(score)
                      ) %>%
            gather(key = "measure", value = "value", width_mean, coverage_mean,
                   score_mean, n
                   ## width_sd, coverage_sd, score_sd
                   ) %>%
            cbind(grid %>% slice(j), .) -> out
        out
    } -> o
    attr(o, "seed") <- seed
    attr(o, "N") <- N
    o
}




## set parameter grid to be evaluated 
grid <- expand.grid(sampleSize = 50,                  # sample size of trial
                    effect = 0.2,                     # average effect, impacts selection bias
                    I2 = c(0, 0.3, 0.6, 0.9),         # Higgin's I^2 heterogeneity measure
                    k = c(2, 3, 5, 10, 20),           # number of studies
                    dist = c("Gaussian", "t"),        # distribution 
                    bias = c("none", "moderate", "strong"),
                    large = c(0, 1, 2),
                    stringsAsFactors = FALSE)        


## run simulation, e.g., on the Rambo server of I-MATH
tic()
out <- sim(grid = grid, N=10000, cores = 120)
toc()


## save results
dir.create("RData", showWarnings=FALSE)
sessionInfo <- sessionInfo()
print(sessionInfo)
save(out, sessionInfo, file = "RData/simulate_all.RData")
