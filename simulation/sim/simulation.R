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
## Felix Hofmann, felix.hofmann2@uzh.ch, Oct. 26, 2021

## Setup
# rm(list = ls())
remotes::install_github("felix-hof/confMeta")
# remotes::install_github("SamCH93/confMeta")
library(confMeta)
library(doParallel) 
library(doRNG)
library(RhpcBLASctl)
blas_set_num_threads(1) # multi threading of BLAS

## Load functions from other files
source("utils.R")
source("study_simulation.R")
source("studies2cis.R")
source("cis2measures.R")
source("measures2summary.R")

################################################################################
#                          Wrap everything in one function                     #
################################################################################

#' Simulate N times, compute CIs, assess CIs
#'
#' Takes a data.frame of parameter configurations, the number replicates
#' \code{N} and the number of cores to simulate accordingly.
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
sim <- function(
    grid,
    N = 1e4,
    cores = detectCores(),
    seed = as.numeric(Sys.time())
) {

    types <- vapply(grid, typeof, character(1L), USE.NAMES = TRUE)

    # check grid columns
    stopifnot(
        is.data.frame(grid),
        c(
            "sampleSize", "I2", "k", "dist",
            "effect", "large", "heterogeneity", "bias"
        ) %in% names(grid),
        types["sampleSize"] %in% c("double", "numeric", "integer"),
        types["effect"] %in% c("double", "numeric", "integer"),
        types["I2"] %in% c("double", "numeric", "integer"),
        types["k"] %in% c("double", "numeric", "integer"),
        types["heterogeneity"] %in% c("character"),
        types["dist"] %in% c("character"),
        types["bias"] %in% c("character"),
        types["large"] %in% c("double", "numeric", "integer")
    )

    # check grid
    stopifnot(
        is.numeric(grid$k),
        is.numeric(grid$sampleSize), all(is.finite(grid$sampleSize)),
        is.numeric(grid$effect), is.numeric(grid$effect),
        all(is.finite(grid$effect)),
        is.numeric(grid$I2), all(0 <= grid$I2), all(grid$I2 < 1),
        is.character(grid$heterogeneity), all(!is.na(grid$heterogeneity)),
        is.character(grid$dist), all(!is.na(grid$dist)),
        is.numeric(grid$large), all(is.finite(grid$large)),
        all(grid$large %in% c(0, 1, 2)),
        is.character(grid$bias), all(!is.na(grid$bias)),
        all(grid$k >= grid$large)
    )

    # check other arguments
    stopifnot(
        is.numeric(N), length(N) == 1L, 1 <= N,
        is.numeric(seed), length(seed) == 1L
    )

    # register parallel backend
    registerDoParallel(cores)
    on.exit(stopImplicitCluster())

    # run simulation
    o <- foreach::foreach(
        j = seq_len(nrow(grid)),
        .options.RNG = seed,
        .errorhandling = "pass"
    ) %dorng% {

        # Since we run this in parallel, we have no other means to check
        # whether an error occured than to check for the existence of the
        # file that is written by the error_function. If it exists, an
        # error happened and we can skip the rest of the iterations.
        if (file.exists("error.txt")) {
            # if error happened, skip the rest of the loop iterations
            out <- "skipped"
        } else {
            cat("start", j, "of", nrow(grid), fill = TRUE)
            pars <- grid[j, ]

            # av is a list with elements that are either a data.frame or NA
            # (the latter in case of an error)
            av <- vector("list", length = N)
            # Also keep track of the average probabilities when simulating
            # publication bias
            pb <- pars$bias != "none"
            p_accept <- if (pb) vector("numeric", length = N) else NULL
            for (i in seq_len(N)) {
                # system(paste0("printf 'j=", j, ", i=", i, "\n'"))
                # Repeat this N times. Simulate studies, calculate CIs,
                # calculate measures
                res <- sim_effects(pars = pars, i = i)
                if (pb) p_accept[i] <- attributes(res)$p_accept
                CIs <- calc_ci(x = res, pars = pars, i = i)
                av[[i]] <- calc_measures(x = CIs, pars = pars, i = i)
            }

            # summarize the N tibbles.
            # If any list element is NA, return "failed".
            if (any(is.na(av))) {
                out <- "failed"
            } else {
                # rbind data frames in list `av`
                df <- do.call("rbind", av)
                attr(df, "N") <- N
                attr(df, "effect") <- pars$effect
                # calculate the summary measures
                out <- calc_summary_measures(
                    df = df,
                    p_accept = p_accept,
                    pars = pars,
                    i = j
                )
                out
            }
        }
        # return output
        out
    }

    # # rbind ci_meas lists together
    o <- do.call("rbind", o)

    # set some attributes and return
    attr(o, "seed") <- seed
    attr(o, "N") <- N
    o
}

################################################################################
#            Calculating the output measures for each individual CI            #
################################################################################

## set parameter grid to be evaluated
grid <- expand.grid(
    # sample size of trial
    sampleSize = 50,
    # average effect, impacts selection bias
    effect = c(0.1, 0.2, 0.5),
    # Higgin's I^2 heterogeneity measure
    I2 = c(0, 0.3, 0.6, 0.9),
    # number of studies
    k = c(3, 5, 10, 20, 50),
    # The heterogeneity model that the studies are simulated from
    heterogeneity = c("additive"),
    # distribution
    dist = c("Gaussian", "sn"),
    # bias
    bias = c("none", "moderate", "strong"),
    # number of large studies
    large = c(0, 1, 2),
    stringsAsFactors = FALSE
)

# For testing
# grid <- grid[floor(seq(1, 1080, length.out = 16)), ]
# N <- 4
# i <- 1
# j <- 588

## run simulation, e.g., on the Rambo server of I-MATH
start <- Sys.time()
out <- sim(grid = grid, N = 2.5e3, cores = 200)
# out <- sim(grid = grid, N = 4, cores = 15)
end <- Sys.time()
print(end - start)

## save results
dir.create("RData", showWarnings = FALSE)
sessionInfo <- sessionInfo()
print(sessionInfo)
save(out, sessionInfo, file = "RData/simulate_all.RData")
