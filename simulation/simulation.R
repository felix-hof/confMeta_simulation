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
rm(list = ls())
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

## Function to generate file paths if intermediate data is saved
generate_path <- function(data_path, pars) {
    file.path(
        data_path,
        paste0(
            paste(
                vapply(
                    seq_along(pars),
                    function(x) paste(names(pars)[x], pars[[x]], sep = "_"),
                    character(1L)
                ),
                collapse = "-"
            ),
            ".rds"
        )
    )
}

# Extract the skewness of the studies from the data frame
extract_skewness <- function(df) {
    meth <- with(df, sub("^.+_(.+)$", "\\1", method))
    skewness <- df$skewness_data
    keep <- !duplicated(meth)
    stats::setNames(
        skewness[keep],
        meth[keep]
    )
}

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
#' @param seed The seed value for the random number generator.
#' @param save_data Indicates whether the individual CIs are stored. If TRUE,
#' all the CIs from each simulation is stored to disk. If FALSE, only the
#' summarized measures for each simulation scenario are saved to disk.
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
    seed = as.numeric(Sys.time()),
    save_data = FALSE) {
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
        is.numeric(seed), length(seed) == 1L,
        isTRUE(save_data) || isFALSE(save_data)
    )

    # register parallel backend
    registerDoParallel(cores)
    on.exit(stopImplicitCluster())

    ## Create the directory where we want to save the results
    dir.create("RData", showWarnings = FALSE)
    if (save_data) {
        data_path <- file.path("RData", "CIs")
        dir.create(data_path, showWarnings = FALSE)
    }

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
            cat("skipping", j, "of", nrow(grid), fill = TRUE)
            out <- "skipped"
        } else {
            cat("starting", j, "of", nrow(grid), fill = TRUE)
            pars <- grid[j, ]

            # av is a list with elements that are either a data.frame or NA
            # (the latter in case of an error)
            av <- vector("list", length = N)
            # Create a second list that we save the intermediate CIs in
            int_cis <- vector("list", length = N)
            # Also keep track of the average probabilities when simulating
            # publication bias
            pb <- pars$bias != "none"
            p_accept <- if (pb) vector("numeric", length = N) else NULL
            for (i in seq_len(N)) {
                # st <- Sys.time()
                # Repeat this N times. Simulate studies, calculate CIs,
                # calculate measures
                res <- sim_effects(pars = pars, i = i)
                if (pb) p_accept[i] <- attributes(res)$p_accept
                CIs <- calc_ci(x = res, pars = pars, i = i)
                if (save_data) {
                    int_cis[[i]] <- list(
                        CI = merge.data.frame(
                            CIs$estimates,
                            CIs$CIs[, c("lower", "upper", "method")],
                            by = "method",
                            all = TRUE,
                            sort = FALSE
                        ),
                        data_skewness = extract_skewness(CIs$CIs),
                        aucc = CIs$aucc,
                        aucc_ratio = CIs$aucc_ratio
                    )
                }
                av[[i]] <- calc_measures(x = CIs, pars = pars, i = i)
            }

            # save the intermediate CIs to disk
            if (save_data) {
                ci_data <- list(pars = pars, cis = int_cis)
                saveRDS(
                    object = ci_data,
                    file = generate_path(data_path = data_path, pars = pars)
                )
            }

            # summarize the N tibbles.
            # If any list element is NA, return "failed".
            if (any(is.na(av))) {
                out <- "failed"
            } else {
                # rbind data frames in list `av`
                df <- do.call("rbind", av)
                # anyNA(df$value)
                # df |> filter(is.na(value)) |> View()
                attr(df, "N") <- N
                attr(df, "effect") <- pars$effect
                attr(df, "median") <- attributes(res)$median
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

    # Keep the RNG attribute since rbind deletes attributes
    rng_att <- attributes(o)$rng

    # print size
    cat(
        paste0(
            "Results size is ",
            format(object.size(o), units = "MB", quote = FALSE),
            "."
        ),
        fill = TRUE
    )

    # Since we experienced some issues with the `save` function, we try to save
    # the output here as well -> but not on math servers as they only have 5GB
    # of disk space
    # saveRDS(o, file = "RData/o.rds")

    # rbind ci_meas lists together
    o <- tryCatch(
        {
            do.call("rbind", o)
        },
        error = function(cond) {
            # set some attributes and return
            attr(o, "seed") <- seed
            attr(o, "rng") <- rng_att
            attr(o, "N") <- N
            save(o, file = "RData/partial_sim.RData")
            NULL
        }
    )

    if (is.null(o)) {
        stop(
            paste0(
                "Error in creation of resulting dataframe. ",
                "Check list elements of the list saved in ",
                "'RData/partial_sim.RData'."
            )
        )
    }

    # set some attributes and return
    attr(o, "seed") <- seed
    attr(o, "rng") <- rng_att
    attr(o, "N") <- N
    o
}

################################################################################
#            Calculating the output measures for each individual CI            #
################################################################################

## set parameter grid to be evaluated
grid <- expand.grid(
    # sample size of trial
    sampleSize = 50L,
    # average effect, impacts selection bias
    effect = 0.2, # c(0.1, 0.2, 0.5),
    # Higgin's I^2 heterogeneity measure
    I2 = c(0, 0.3, 0.6, 0.9),
    # number of studies
    k = c(3L, 5L, 10L, 20L, 50L),
    # The heterogeneity model that the studies are simulated from
    heterogeneity = c("additive"),
    # distribution
    dist = c("Gaussian", "snr", "snl"),
    # bias
    bias = "none",
    # bias = c("none", "moderate", "strong"),
    # number of large studies
    large = c(0L, 1L, 2L),
    stringsAsFactors = FALSE
)

# j <- which(
#     with(
#         grid,
#         {
#             effect == 0.2 &
#                 I2 == 0.9 &
#                 k == 3L &
#                 dist == "snr" &
#                 bias == "none" &
#                 large == 1L
#         }
#     )
# )
# seed <- 42
# save_data <- TRUE

# Set some parameters based on available computing machines
machine <- Sys.info()["nodename"]
if (machine == "T14s") {
    # On my machine, run this with N = 5, on 14 cores, and on a smaller grid.
    # This is only used for development, debugging, and testing.
    N <- 50
    cores <- 15
    # grid <- grid[floor(seq(1, 1080, length.out = 160)), ]
} else if (machine == "david") {
    # The math institute has a server called "david" with 80 CPU cores
    N <- 1e4
    cores <- 60
} else if (machine == "rambo") {
    # The math institute has a server called "rambo" with 128 CPU cores
    N <- 1e4
    cores <- 124
} else if (machine == "box") {
    N <- 2e4
    cores <- 230
} else {
    stop(
        paste0(
            "Unknown host. Please add this machine to the if-statement",
            " and specify the amount of times each scenario should be run",
            " as well as the number of CPU cores."
        )
    )
}

# i <- 1
# j <- 5

## run simulation, e.g., on the Rambo server of I-MATH
cat(paste0("Running every scenario ", N, " times."), fill = TRUE)
cat(
    paste0("In total there are ", nrow(grid), " simulation scenarios."),
    fill = TRUE
)
cat(paste0("Simulation will be run on ", cores, " CPU cores."), fill = TRUE)
start <- Sys.time()
out <- sim(grid = grid, N = N, cores = cores, save_data = TRUE)
end <- Sys.time()
run_time <- end - start
cat(
    paste0(
        "Running simulation required ",
        round(run_time, 2),
        " ",
        attributes(run_time)$units,
        "."
    ),
    fill = TRUE
)
attr(out, which = "runtime") <- run_time

# Check some of the results as cause of the error:
# "Error in attr(o, "seed") <- seed : attempt to set attribute on NULL"
# This error happens when in the call to rbind when we try to create the data
# frame from a list. It turned out that the reason for this was that in the
# cases when one or more of these CIs did not exist, the confMeta function
# 'get_ci' which did return a misspecified list.

# files <- list.files(path = "RData/CIs", full.names = TRUE)
# has_na <- logical(length(files))
# na_cis <- list()
# for (i in seq_along(files)) {
#     dat <- readRDS(file = files[i])
#     cis <- lapply(dat$cis, "[[", i = "CI")
#     idx <- sapply(cis, \(x) anyNA(x$lower))
#     if (any(idx)) {
#         na_cis <- append(na_cis, cis[idx])
#         has_na[i] <- TRUE
#     }
# }


## save results
sessionInfo <- sessionInfo()
print(sessionInfo)
save(out, sessionInfo, file = "RData/simulate_all.RData")
