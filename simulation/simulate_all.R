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
source("ReplicationSuccess_extension_LH.R")
library(tidyverse); theme_set(theme_bw())
library(rlang)
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
#' @param heterogeneity The heterogeneity model, the studies are simulated from.
#' Either "additive" or "multiplicative".
#' @param dist distribution to simulate the study effect. Either "t" or "Gaussian".
#' the sample size as specified by \code{sampleSize}.
#' @return a matrix \code{k} x 2 matrix with columns
#' \code{theta} (effect estimates) and
#' \code{se} (standard errors).
simRE <- function(k, sampleSize, effect, I2, 
                  heterogeneity = c("additive", "multiplicative"),
                  dist=c("t", "Gaussian"),
                  large) {
    
    # get args
    dist <- match.arg(dist)
    n <- rep(sampleSize, k)
    heterogeneity <- match.arg(heterogeneity)
    # include large studies
    if(large == 1)
        n[1] <- 10 * n[1]
    if(large == 2)
        n[1:2] <- 10 * n[1:2]
    # stuff for additive model
    if(heterogeneity == "additive"){
        eps2 <- 1/k * sum(2/n)
        tau2 <- eps2 * I2/(1 - I2)
        if(dist == "t") {
            ## the sn::rst(xi=0, omega, nu) distribution has variance
            ## omega^2 nu/(nu-2) (if nu>2)
            ## where nu is the degrefees of freedom (dof).
            ## So if we want the variance to be tau^2, then 
            ## omega^2 = tau^2 * (nu-2)/nu
            ## We use nu=4 dof then omega^2 = tau^2/2, so half as
            ## large as the heterogeneity variance under normality. 
            delta <- rst(n = k, xi = effect, omega = sqrt(tau2/2), nu = 4)
        } else {
            delta <- rnorm(n = k, mean = effect, sd = sqrt(tau2))
        }
        theta <- rnorm(n = k, mean = delta, sd = sqrt(2/n))
        ## theta[1:ceiling(r*k)] <- theta[1:ceiling(r*k)] + bias
        se <- sqrt(rchisq(n = k, df = 2*n - 2) / (n*(n - 1)))
        o <- cbind("theta" = theta, "se" = se, "delta" = delta)
    } else { ## multiplicative model
        phi <- 1/(1 - I2)
        if(dist == "t") {
            ## the sn::rst(xi=0, omega, nu) distribution has variance
            ## omega^2 nu/(nu-2) (if nu>2)
            ## where nu is the degrees of freedom (dof).
            ## So if we want the variance to be tau^2, then 
            ## omega^2 = tau^2 * (nu-2)/nu
            ## We use nu=4 dof then omega^2 = tau^2/2, so half as
            ## large as the heterogeneity variance under normality.
            ## sample sequentially with marginal variance equal to
            ## (phi-1)*2/n + 2/n = phi*2/n 
            delta <- rst(n = k, xi = effect, omega = sqrt((phi-1)/n), nu = 4)
            theta <- rnorm(n = k, mean = delta, sd = sqrt(2/n))
        } else {  ## Gaussian, sample directly from marginal
            theta <- rnorm(n = k, mean = effect, sd = sqrt(phi*2/n))
        }
        ## theta[1:ceiling(r*k)] <- theta[1:ceiling(r*k)] + bias
        se <- sqrt(rchisq(n = k, df = 2*n - 2) / (n*(n - 1)))
        o <- cbind("theta" = theta, "se" = se, "delta" = NA_real_)
    }
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
#' #' @param heterogeneity The heterogeneity model, the studies are simulated from.
#' Either "additive" or "multiplicative".
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
simREbias <- function(k, sampleSize, effect, I2,
                      heterogeneity = c("additive", "multiplicative"),
                      dist = c("t", "Gaussian"), 
                      large,
                      bias = c("none", "moderate", "strong"),
                      verbose = TRUE){
    # input checks
    stopifnot(length(k) == 1,
              length(sampleSize) == 1,
              is.numeric(effect),
              length(effect) == 1,
              is.finite(effect),
              length(I2) == 1,
              0 <= I2, I2 < 1,
              length(heterogeneity) == 1,
              !is.null(dist),
              !is.null(large),
              is.numeric(large),
              length(large)==1,
              large %in% c(0,1,2),
              length(bias) == 1,
              !is.null(bias))
    
    bias <- match.arg(bias)

    if(bias == "none"){
        o <- simRE(k=k, sampleSize=sampleSize, effect = effect, I2=I2, heterogeneity=heterogeneity, dist=dist, large=large)
        attr(o, "heterogeneity") <- heterogeneity
        attr(o, which = "effect") <- effect
        return(o)
    }
        
    ## first ignore the 'large' 
    o <- simRE(k=k*3, sampleSize=sampleSize, effect=effect, I2=I2, heterogeneity=heterogeneity, dist=dist, large=0)
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
        oLarge <- simRE(k=k, sampleSize=sampleSize, effect = effect, I2=I2, heterogeneity=heterogeneity, dist=dist, large=large)
        o <- rbind(oLarge[1,], o[-1,]) 
    }
    if(large == 2){
        oLarge <- simRE(k=k, sampleSize=sampleSize, effect = effect, I2=I2, heterogeneity=heterogeneity, dist=dist, large=large)
        o <- rbind(oLarge[1:2,], o[-c(1,2),]) 
    }
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
    ## Henmy & Copas confidence Interval
    HC <- metafor::hc(object = metafor::rma(yi = x[, "theta"], sei = x[, "se"], 
                                            control = list(maxiter = 1000, stepadj = 0.5)))

    ## standard metagen with REML estimation of tau
    REML <- metagen(TE = x[, "theta"], seTE = x[, "se"], sm = "MD", 
                    method.tau = "REML", 
                    control = list(maxiter = 1000, stepadj = 0.5))

    ## Hartung & Knapp
    HK <- metagen(TE = x[, "theta"], seTE = x[, "se"], sm = "MD", 
                  method.tau = "REML", hakn = TRUE,
                  control = list(maxiter = 1000, stepadj = 0.5))

    ## HMean2sided
    if(nrow(x) <= 5) {
        HM2 <- hMeanChiSqCI(thetahat = x[, "theta"], se = x[, "se"], 
                            alternative = "two.sided")
    } else {
        HM2 <- list(CI = cbind(lower = NA_real_, upper = NA_real_))
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
    
    tib <- tibble(lower = c(HC$ci.lb,
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
    out <- list(CIs = tib,
                model = attributes(x)$heterogeneity,
                gamma = tibble("method" = c("Harmonic Mean", "Harmonic Mean Additive", "Harmonic Mean Multiplicative"),
                               "gamma_min" = c(min(HM$gamma[,2]), min(HM_tau2$gamma[,2]), min(HM_phi$gamma[,2])),
                               "x_gamma_min" = c(HM$gamma[which.min(HM$gamma[,2]),1], HM_tau2$gamma[which.min(HM_tau2$gamma[,2]),1], 
                                                 HM_phi$gamma[which.min(HM_phi$gamma[,2]),1])),
                theta = x[, "theta"],
                delta = x[, "delta"],
                effect = attributes(x)$effect)
    out
}



#' Computes quality measures for CIs
#'
#' @param x a list with elements \code{CIs}, \code{model}, \code{gamma}, \code{theta}, 
#' \code{delta} and \code{effect} as obtained from \code{sim2CIs}.
#' @return a tibble with columns 
#' \item{\code{method}}{method}
#' \item{\code{width}}{with of the intervals}
#' \item{\code{coverage}}{covarage of the true value 0}
#' \item{\code{score}}{interval score as defined in Gneiting and Raftery (2007)}
#' \item{\code{coverage_effects}}{Proportion of study effects covered by the interval(s).}
#' \item{\code{n}}{Number of intervals}
CI2measures <- function(x) {
    methods <- unique(x$CIs$method)
    foreach(i = seq_along(methods), .combine = rbind) %do% {
        x$CIs %>% filter(method == methods[i]) %>%
            select(lower, upper) %>%
            as.matrix() ->
            x_sub
        
        if(x$model == "additive"){
            vapply(x$theta, function(theta){
                any(x_sub[,"lower"] <= theta & theta <= x_sub[,"upper"])
            }, logical(1L)) %>%
                mean() ->
                coverage_effects
        } else {
            coverage_effects <- NA_real_
        }
        
        if(methods[i] %in% paste0("Harmonic Mean", c("", " Additive", " Multiplicative"))){
            gamma_min <- x$gamma %>% filter(method == methods[i]) %>% pull(gamma_min)
        } else {
            gamma_min <- NA_real_
        }

        { x_sub[,"upper"] - x_sub[,"lower"] } %>%
            sum(.) ->
            width
        
        as.numeric(any(x_sub[,"lower"] <= x$effect & x$effect <= x_sub[,"upper"])) ->
            coverage
        
        width + (2/0.05) * min(abs(x_sub[,"lower"]), abs(x_sub[,"upper"])) * (1 - coverage) ->
            score
        
        if(all(is.na(x_sub))){
            n <- NA_real_
        } else {
            n <- nrow(x_sub)
        }
        
       out <- tibble(method = methods[i], 
                     coverage = coverage, coverage_effects = coverage_effects,
                     gammaMin = gamma_min, n = n,
                     score = score, width = width)
        
        # return
        out
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
              c("sampleSize", "I2", "k", "dist", 
                "effect", "large", "heterogeneity", "bias") %in% names(grid),
              is.numeric(N), length(N) == 1, 1 <= N)
    
    registerDoParallel(cores)
    on.exit(stopImplicitCluster())
    
    foreach(j = seq_len(nrow(grid)), .options.RNG=seed, .errorhandling = "pass") %dorng% {
        
        if(file.exists("error.txt")){
            # if error happened, skip the rest of the loop iterations
            out <- "skipped"
        } else {
            
            cat("start", j, "of", nrow(grid), fill=TRUE)
            grid %>% slice(j) -> pars
            
            # av is either a tibble or NULL
            av <- tryCatch({
                foreach(i = seq_len(N), .combine = rbind) %do% {
                    res <- simREbias(k = pars$k, sampleSize = pars$sampleSize,
                                     effect = pars$effect, I2 = pars$I2,
                                     heterogeneity = pars$heterogeneity,
                                     dist = pars$dist, bias = pars$bias,
                                     large = pars$large)
                    CIs <- sim2CIs(x = res)
                    CI2measures(x = CIs)
                } -> av
                av
            }, 
            error = function(cond){
                text = capture.output(cond)
                cat("Parameters are:\n",
                    paste0(paste0(names(pars), ":", pars[1, ]), collapse = "\n"),
                    "\n\n", "The error message is:\n",
                    text, "\n\n",
                    file = "error.txt", append = TRUE)
                return(NULL)
            })
            
            # out is either a tibble or "failed"
            if(is.null(av)){
                out <- "failed"
            } else {
                out <- tryCatch({
                    ## compute mean for everything
                    bind_rows(
                        # mean for all measures
                        av %>% group_by(method) %>%
                            summarize(across(everything(), 
                                             .fns = list(mean = function(x){if(any(is.na(x))){NA_real_} else {mean(x, na.rm = TRUE)}}),
                                             .names = "{.col}_{.fn}"),
                                      .groups = "drop"
                            ) %>% 
                            pivot_longer(cols = !method,
                                         names_to = "measure",
                                         values_to = "value"
                                         ## width_sd, coverage_sd, score_sd
                            ) %>%
                            cbind(pars, .),
                        
                        # summary statistics for gamma_min
                        av %>% filter(method %in% c("Harmonic Mean", "Harmonic Mean Additive", "Harmonic Mean Multiplicative")) %>%
                            group_by(method) %>%
                            summarize(across("gammaMin", 
                                             .fns = list(min = function(x){if(any(is.na(x))){NA_real_} else {min(x, na.rm = TRUE)}}, 
                                                         firstQuart = function(x){if(any(is.na(x))){NA_real_} else {quantile(x, probs = 0.25, na.rm = TRUE, names = FALSE)}},
                                                         median = function(x){if(any(is.na(x))){NA_real_} else {median(x, na.rm = TRUE)}},
                                                         mean = function(x){if(any(is.na(x))){NA_real_} else {mean(x, na.rm = TRUE)}},
                                                         thirdQuart = function(x){if(any(is.na(x))){NA_real_} else {quantile(x, probs = 0.75, na.rm = TRUE, names = FALSE)}},
                                                         max = function(x){if(any(is.na(x))){NA_real_} else {max(x, na.rm = TRUE)}}),
                                             .names = "{.col}_{.fn}"),
                                      .groups = "drop"
                            ) %>%
                            pivot_longer(cols = !method,
                                         names_to = "measure",
                                         values_to = "value"
                            ) %>%
                            cbind(pars, .),
                        
                        # relative frequency for n
                        av %>% filter(method %in% c("Harmonic Mean", "Harmonic Mean Additive", "Harmonic Mean Multiplicative")) %>%
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
                                      .groups = "drop"
                            ) %>%
                            pivot_longer(cols = !method,
                                         names_to = "measure",
                                         values_to = "value"
                            ) %>%
                            filter(!is.na(value)) %>% 
                            cbind(pars, .)
                    )
                },
                error = function(cond){
                    text <- capture.output(cond)
                    cat("Parameters are:\n",
                        paste0(paste0(names(pars), ":", pars[1, ]), collapse = "\n"),
                        "\n\n", "The error message is:\n",
                        text, "\n\n",
                        "The av object has been saved to 'error.rds'.",
                        file = "error.rds", append = TRUE)
                    saveRDS(av, file = "error.rds")
                    return("failed")
                })
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
                    k = c(2, 3, 5, 10, 20),                          # number of studies
                    heterogeneity = c("additive", "multiplicative"), # The heterogeneity model that the studies are simulated from
                    dist = c("Gaussian", "t"),                       # distribution 
                    bias = c("none", "moderate", "strong"),          # bias
                    large = c(0, 1, 2),                              # number of large studies
                    stringsAsFactors = FALSE) %>% 
    arrange(desc(heterogeneity))


## run simulation, e.g., on the Rambo server of I-MATH
tic()
out <- sim(grid = grid, N=10000, cores = 120)
toc()


## save results
dir.create("RData", showWarnings=FALSE)
sessionInfo <- sessionInfo()
print(sessionInfo)
save(out, sessionInfo, file = "RData/simulate_all.RData")


