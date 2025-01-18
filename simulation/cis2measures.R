################################################################################
#            Calculating the output measures for each individual CI            #
################################################################################

# Calculate the skewness from the confidence intervals
calc_ci_skewness <- function(cis, estimates, ci_exists) {
    if (!ci_exists) {
        NA_real_
    } else {
        up <- cis[, "upper"]
        low <- cis[, "lower"]
        est <- estimates
        (up + low - 2 * est) / (up - low)
    }
}

# Dummy function for adding the skewness calculated from the data
# this was calculated before already
add_data_skewness <- function(data_skewness) {
    data_skewness
}

calc_coverage_effects <- function(cis, delta, ci_exists) {
    # if ci does not exist
    if (!ci_exists) {
        0
    } else {
        mean(
            vapply(
                delta,
                function(d) {
                    any(cis[, "lower"] <= d & d <= cis[, "upper"])
                },
                logical(1L)
            )
        )
    }
}

# calculate how many times at least one study-specific effect
# is covered
calc_coverage_effects_min1 <- function(cis, delta, ci_exists) {
    if (!ci_exists) {
        0
    } else {
        found <- FALSE
        for (d in delta) {
            if (any(cis[, "lower"] <= d & d <= cis[, "upper"])) {
                found <- TRUE
                break
            }
        }
        as.numeric(found)
    }
}

# calculate whether all deltas covered by CI
calc_coverage_effects_all <- function(cis, delta, ci_exists) {
    if (!ci_exists) {
        0
    } else {
        as.numeric(
            all(
                vapply(
                    delta,
                    function(d) {
                        any(
                            cis[, "lower"] <= d &
                            d <= cis[, "upper"]
                        )
                    },
                    logical(1L)
                )
            )
        )
    }
}

# calculate whether interval covers future study
# (only for hMean, hMean_additive, hMean_multiplicative)
calc_coverage_prediction <- function(cis, pars, ci_exists) {
    if (!ci_exists) {
        0
    } else {
        new_study <- simREbias(
            k = 1, sampleSize = pars$sampleSize,
            effect = pars$effect, I2 = pars$I2,
            heterogeneity = pars$heterogeneity, dist = pars$dist, large = 0,
            bias = "none"
        )[, "delta"]
        as.numeric(
            any(
                cis[, "lower"] <= new_study & new_study <= cis[, "upper"]
            )
        )
    }
}

# calculate the width
calc_width <- function(cis, ci_exists) {
    if (!ci_exists) {
        0
    } else {
        sum(cis[, "upper"] - cis[, "lower"])
    }
}

# calculate wether the ci covers the true effect
calc_coverage_true <- function(cis, effect, ci_exists) {
    if (!ci_exists) {
        0
    } else {
        as.numeric(any(cis[, "lower"] <= effect & effect <= cis[, "upper"]))
    }
}

# return gamma min
calc_gamma_min <- function(gamma) {
    if (all(is.na(gamma))) {
        NA_real_
    } else {
        gamma[, "gamma_min"]
    }
}

# return the score
calc_score <- function(cis, effect, ci_exists) {
    if (!ci_exists) {
        NA_real_
    } else {
        abs_ci <- abs(cis)
        calc_width(cis = cis, ci_exists = ci_exists) + (2 / 0.05) +
        min(abs_ci[, "lower"], abs_ci[, "upper"]) *
        1 - calc_coverage_true(cis, effect, ci_exists)
    }
}

# returns the number of CIs
calc_n <- function(cis, ci_exists) {
    if (!ci_exists) {
        0
    } else {
        nrow(cis)
    }
}

# returns the squared difference between the estimate and the effect
# in case of multiple effect estimates, NA is returned
calc_sq_diff <- function(estimates, effect) {
    if (length(estimates) != 1L) {
        NA_real_
    } else {
        (effect - estimates)^2
    }
}

# if there is only one estimate, returns it, otherwise returns NA
calc_estimate <- function(estimates) {
    if (length(estimates) != 1L) {
        NA_real_
    } else {
        estimates
    }
}

# based on the current method, what measures do we need to calculate
get_measures <- function(is_ci, is_pi, is_new) {

    # list all functions for measures
    all_meas <- list(
        coverage_mean = quote({
            calc_coverage_true(
                cis = cis,
                effect = mean,
                ci_exists = ci_exists
            )
        }),
        coverage_median = quote({
            calc_coverage_true(
                cis = cis,
                effect = median,
                ci_exists = ci_exists
            )
        }),
        coverage_effects = quote({
            calc_coverage_effects(
                cis = cis,
                delta = delta,
                ci_exists = ci_exists
            )
        }),
        coverage_effects_min1 = quote({
            calc_coverage_effects_min1(
                cis = cis,
                delta = delta,
                ci_exists = ci_exists
            )
        }),
        coverage_effects_all = quote({
            calc_coverage_effects_all(
                cis = cis,
                delta = delta,
                ci_exists = ci_exists
            )
        }),
        coverage_prediction = quote({
            calc_coverage_prediction(
                cis = cis,
                pars = pars,
                ci_exists = ci_exists
            )
        }),
        ci_skewness = quote({
            calc_ci_skewness(
                cis = cis,
                estimates = estimates,
                ci_exists = ci_exists
            )
        }),
        data_skewness = quote({
            add_data_skewness(data_skewness = data_skewness)
        }),
        width = quote({
            calc_width(
                cis = cis,
                ci_exists = ci_exists
            )
        }),
        gamma_min = quote({
            calc_gamma_min(gamma = gamma)
        }),
        score = quote({
            calc_score(
                cis = cis,
                effect = effect,
                ci_exists = ci_exists
            )
        }),
        n = quote({
            calc_n(
                cis = cis,
                ci_exists = ci_exists
            )
        }),
        p_max = quote({
            unique(p_max)
        }),
        aucc = quote({
            aucc
        }),
        aucc_ratio = quote({
            aucc_ratio
        }),
        mse = quote({
            calc_sq_diff(
                estimates = estimates,
                effect = effect
            )
        }),
        estimate = quote({
            calc_estimate(
                estimates = estimates
            )
        })
    )

    # Allocate memory for list
    calc_meas <- vector("list", length = 3L)
    counter <- 1L
    if (is_ci) {
        calc_meas[[counter]] <- c(
            "coverage_mean",
            "coverage_median",
            # "coverage_effects",
            # "coverage_effects_min1",
            # "coverage_effects_all",
            # "ci_skewness",
            # "data_skewness",
            "width",
            "score",
            # "mse",
            "estimate"
        )
        counter <- counter + 1L
    }
    if (is_pi) {
        calc_meas[[counter]] <- c(
            "coverage_prediction"
        )
        counter <- counter + 1L
    }
    if (is_new) {
        calc_meas[[counter]] <- c(
            "n",
            # "gamma",
            "p_max",
            "aucc",
            "aucc_ratio",
            "ci_skewness",
            "data_skewness"
        )
        counter <- counter + 1L
    }
    calc_meas <- do.call("c", calc_meas)

    # return a list with appropriate functions
    all_meas[calc_meas]
}

#' Computes quality measures for CIs
#'
#' @param x a list with elements \code{CIs}, \code{model}, \code{gamma},
#' \code{theta}, \code{delta} and \code{effect} as obtained from
#' \code{sim2CIs}.
#' @param pars a \code{data.frame} with one row. It should have column names
#' \code{k}, \code{sampleSize}, \code{effect}, \code{I2}, \code{heterogeneity},
#' \code{dist} and \code{large}. These are passed as parameters to
#' \link{simREbias()}. In order to simulate another study in order to assess
#' prediction interval coverage.
#' @return a tibble with columns
#' \item{\code{method}}{method}
#' \item{\code{width}}{with of the intervals}
#' \item{\code{coverage}}{covarage of the true value 0}
#' \item{\code{score}}{interval score as defined in Gneiting and Raftery (2007)}
#' \item{\code{coverage_effects}}{Proportion of study effects covered by
#' the interval(s).}
#' \item{\code{n}}{Number of intervals}
CI2measures <- function(x, pars) {

    # get the method corresponding to each row of the CIs
    ci_df <- x$CIs
    ci_method <- ci_df$method
    # get the indices of unique methods and their names
    unique_methods_idx <- !duplicated(ci_method)
    unique_methods <- ci_method[unique_methods_idx]
    # for each method, assess whether it is a CI, whether it is new,
    # whether it is a PI
    ci_exists <- ci_df$ci_exists[unique_methods_idx]
    is_ci <- ci_df$is_ci[unique_methods_idx]
    is_pi <- ci_df$is_pi[unique_methods_idx]
    is_new <- ci_df$is_new[unique_methods_idx]
    # convert the CIs to matrix
    cis <- as.matrix(ci_df[c("lower", "upper")])
    # # get gamma
    # gamma_methods <- x$gamma$method
    # gamma <- as.matrix(x$gamma[c("gamma_min", "x_gamma_min")])
    # get p_max
    p_max <- x$p_max$p_max
    p_max_method <- x$p_max$method
    # get estimates
    estimates <- x$estimates$estimate
    estimates_method <- x$estimates$method
    # get aucc
    aucc <- x$aucc$aucc
    aucc_method <- x$aucc$method
    # get aucc_ratio
    aucc_ratio <- x$aucc_ratio$aucc_ratio
    aucc_ratio_method <- x$aucc_ratio$method
    # get the deltas
    delta <- x$delta
    # get the effect (mean)
    effect <- x$effect
    # get the median
    median <- x$median
    # get the skewness from the data (gamma)
    data_skewness <- ci_df$skewness_data

    out <- vector("list", length(unique_methods))
    # loop over each of the methods
    for (k in seq_along(unique_methods)) {
        # Get the current method and characteristics
        curr_method <- unique_methods[k]
        curr_is_ci <- is_ci[k]
        curr_is_pi <- is_pi[k]
        curr_is_new <- is_new[k]
        # Get the unevaluated function calls for the current method
        curr_meas <- get_measures(
            is_ci = curr_is_ci,
            is_pi = curr_is_pi,
            is_new = curr_is_new
        )
        # create the argument list
        arg_list <- list(
            cis = cis[ci_method == curr_method, , drop = FALSE],
            mean = effect,
            median = median,
            ci_exists = ci_exists[k],
            delta = delta,
            pars = pars,
            data_skewness = data_skewness[k],
            p_max = if (curr_is_new) {
                p_max[p_max_method == curr_method]
            } else {
                NULL
            },
            aucc = if (curr_is_new) {
                aucc[aucc_method == curr_method]
            } else {
                NULL
            },
            aucc_ratio = if (curr_is_new) {
                aucc_ratio[aucc_ratio_method == curr_method]
            } else {
                NULL
            },
            # pass estimate always here
            # if more than one estimate, return NA in calc_diff function
            estimates = estimates[estimates_method == curr_method]
            # gamma = gamma[gamma_methods == curr_method, , drop = FALSE]
        )
        # Calculate the measures
        out[[k]] <- data.frame(
            value = vapply(
                curr_meas,
                function(expr, arg_list) {
                    eval(expr, envir = arg_list)
                },
                arg_list = arg_list,
                numeric(1L)
            ),
            measure = names(curr_meas),
            method = curr_method,
            is_ci = curr_is_ci,
            is_pi = curr_is_pi,
            is_new = curr_is_new,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    }

    do.call(
        "rbind",
        out
    )
}

## wrapper function across CI2measures with error handling
calc_measures <- function(x, pars, i) {
    # if error in function before, just return NA back.
    # Otherwise, run sim2CIs. If an error happens, log everything and
    # return NA
    if (length(x) == 1L && is.na(x)) {
        NA
    } else {
        tryCatch(
            {
                CI2measures(x = x, pars = pars)
            },
            error = function(cond) {
                error_function(
                    cond = cond,
                    pars = pars,
                    error_obj = x,
                    fun_name = "calc_measures",
                    i = i
                )
                NA
            }
        )
    }
}
