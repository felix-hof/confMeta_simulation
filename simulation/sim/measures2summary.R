################################################################################
#            Calculate summary measures from the individual measures           #
################################################################################

# Calculate the mean measure for each of the method and measure subgroups
get_mean_stats <- function(df) {
    # Summarise all the measures except gamma_min & p_max
    df_sub <- subset(
        df,
        !(measure %in% c("gamma_min", "p_max", "mse", "estimate", "ci_skewness"))
    )
    # This works since mean of a vector of all NAs is just NA
    fun <- function(x) {
        c("value" = mean(x, na.rm = TRUE), "prop" = sum(!is.na(x)) / length(x))
    }
    mean_stats <- stats::aggregate(
        value ~ measure + method + is_ci + is_pi + is_new,
        FUN = fun,
        data = df_sub,
        na.action = identity
    )
    mean_stats <- within(
        mean_stats,
        {
            prop <- value[, "prop"]
            value <- value[, "value"]
        }
    )
    mean_stats$usage <- "mean_plot"
    mean_stats$stat_fun <- "mean"
    mean_stats
}

# Calculate the summary measures for gamma_min
get_gamma_stats <- function(df, col_order) {
    # Run all the functions of f_list_gamma
    # To the correct subset of df and add
    # some information regarding stat and
    # usage
    df_sub <- subset(df, measure %in% c("gamma_min", "p_max"))
    # list of functions that calculate the desired summary measures
    f_list_gamma <- list(
        min = function(x) {
            c(
                "value" = min(x, na.rm = TRUE),
                "prop" = sum(!is.na(x)) / length(x)
            )
        },
        firstQuart = function(x) {
            c(
                "value" = quantile(x, probs = 0.25, names = FALSE),
                "prop" = sum(!is.na(x)) / length(x)
            )
        },
        median = function(x) {
            c(
                "value" = median(x, na.rm = TRUE),
                "prop" = sum(!is.na(x)) / length(x)
            )
        },
        mean = function(x) {
            c(
                "value" = mean(x, na.rm = TRUE),
                "prop" = sum(!is.na(x)) / length(x)
            )
        },
        thirdQuart = function(x) {
            c(
                "value" = quantile(x, probs = 0.75, names = FALSE),
                "prop" = sum(!is.na(x)) / length(x)
            )
        },
        max = function(x) {
            c(
                "value" = max(x, na.rm = TRUE),
                "prop" = sum(!is.na(x)) / length(x)
            )
        }
    )
    # loop over the functions and calculate the summary measures
    gamma_stats <- lapply(
        seq_along(f_list_gamma),
        function(i) {
            res <- stats::aggregate(
                value ~ measure + method + is_ci + is_pi + is_new,
                FUN = f_list_gamma[[i]],
                data = df_sub,
                na.action = identity
            )
            res <- within(
                res,
                {
                    prop <- value[, "prop"]
                    value <- value[, "value"]
                }
            )
            res$usage <- "summary_stat_plot"
            res$stat_fun <- names(f_list_gamma)[i]
            res
        }
    )
    gamma_stats <- do.call("rbind", gamma_stats)
    gamma_stats
}

# Calculate the summary measures for n
# We want to count the number of intervals:
# This function generates a list of functions
# that check how many times the number of
# intervals is equal to n. Here, n is vectorised
# and there is also a function that counts the number
# of times where the number of intervals exceeds
# max(n).
make_f_n <- function(n) {
    maxn <- max(n)
    f_n <- lapply(
        n,
        function(one_n) {
            # force(one_n)
            function(x) {
                c(
                    "value" = sum(x == one_n),
                    "prop" = sum(!is.na(x)) / length(x)
                )
            }
        }
    )
    f_n <- append(
        f_n,
        list(
            function(x) {
                c(
                    "value" = sum(x > maxn),
                    "prop" = sum(!is.na(x)) / length(x)
                )
            }
        )
    )
    names(f_n) <- c(as.character(n), paste0("gt", max(n)))
    f_n
}

# Count how many times the number of intervals is equal to n,
# where n is the same as the input to `make_f_n`.
get_n_stats <- function(df, col_order) {
    # get the counting functions
    f_list_n <- make_f_n(0:9)
    # get the number of intervals
    df_sub <- subset(df, measure == "n")
    # loop over the counting functions
    n_stats <- lapply(
        seq_along(f_list_n),
        function(i, f_list_n, df_sub) {
            res <- stats::aggregate(
                value ~ measure + method + is_ci + is_pi + is_new,
                FUN = f_list_n[[i]],
                data = df_sub,
                na.action = identity
            )
            res <- within(
                res,
                {
                    prop <- value[, "prop"]
                    value <- value[, "value"]
                }
            )
            res$usage <- "n_plot"
            res$stat_fun <- names(f_list_n)[i]
            res
        },
        f_list_n = f_list_n,
        df_sub = df_sub
    )
    n_stats <- do.call("rbind", n_stats)
    n_stats[col_order]
}

# This function calculates the MSE, variance and bias
# the estimates
get_bias_var_stats <- function(df, col_order) {
    # get the effect
    effect <- attributes(df)$effect
    # store the method specific info
    df_sub_rest <- unique(subset(df, select = -c(measure, value)))
    # get the point estimates for all methods and convert to factor
    df_sub_est <- subset(df, measure == "estimate", select = c(method, value))
    df_sub_est$method <- factor(df_sub_est$method)
    # Calculate the bias, i.e. mean_estimate - true_effect
    mean_est <- stats::aggregate(
        value ~ method,
        FUN = function(x) {
            c(
                "value" = mean(x, na.rm = TRUE),
                "prop" = sum(!is.na(x)) / length(x)
            )
        },
        data = df_sub_est,
        na.action = identity
    )
    mean_est <- within(
        mean_est,
        {
            prop <- value[, "prop"]
            value <- value[, "value"]
        }
    )
    bias <- within(mean_est, value <- value - effect)
    # Calculate variance, i.e. var(estimates)
    var_est <- stats::aggregate(
        value ~ method,
        FUN = function(x) {
            c(
                "value" = if (all(is.na(x))) {
                    NA_real_
                } else {
                    # similar to var(x, na.rm = TRUE) but
                    # uses 1/n instead of 1/(n-1)
                    y <- x[!is.na(x)]
                    n <- length(y)
                    n^(-1) * sum((y - mean(y))^2)
                },
                "prop" = sum(!is.na(x)) / length(x)
            )
        },
        data = df_sub_est,
        na.action = identity
    )
    var_est <- within(
        var_est,
        {
            prop <- value[, "prop"]
            value <- value[, "value"]
        }
    )
    # Calculate MSE, i.e. (1/n) * sum((estimate - true_effect)^2)
    mse_est <- stats::aggregate(
        value ~ method,
        FUN = function(x, effect) {
            c(
                "value" = mean((x - effect)^2, na.rm = TRUE),
                "prop" = sum(!is.na(x)) / length(x)
            )
        },
        data = df_sub_est,
        effect = effect,
        na.action = identity
    )
    mse_est <- within(
        mse_est,
        {
            prop <- value[, "prop"]
            value <- value[, "value"]
        }
    )
    # Check
    # test <- Reduce(
    #     f = function(x, y) merge(x, y, all = TRUE, by = "method"),
    #     x = list(bias, var_est, mse_est)
    # )
    # names(test)[2:4] <- c("bias", "var", "mse")
    # out <- cbind(sum = with(test, var + bias^2), mse =  test$mse)
    # out[, 1L] - out[, 2L]
    ll <- list("bias" = bias, "var" = var_est, "mse" = mse_est)
    do.call(
        "rbind",
        lapply(
            X = seq_along(ll),
            FUN = function(i, ll, df_sub_rest, col_order) {
                nm <- names(ll)[i]
                d <- cbind(data.frame(measure = nm), df_sub_rest)
                out <- merge(d, ll[[i]], by = "method", sort = FALSE)
                out$usage <- "mean_plot"
                out$stat_fun <- nm
                out[col_order]
            },
            ll = ll,
            df_sub_rest = df_sub_rest,
            col_order = col_order
        )
    )
}

# This function calculates the mean probability of acceptance.
# This is only important when we simulate with publication bias.
# If there is no publication bias, `p_accept` will be NULL
get_paccept_stats <- function(p_accept, col_order) {
    if (is.null(p_accept)) {
        NULL
    } else {
        data.frame(
            measure = "p_accept",
            method = "none",
            is_ci = FALSE,
            is_pi = FALSE,
            is_new = FALSE,
            value = mean(p_accept),
            prop = 1,
            usage = "paccept_plot",
            stat_fun = "mean",
            stringsAsFactors = FALSE,
            row.names = NULL
        )[col_order]
    }
}

get_skewness_stats <- function(df, pars, col_order) {
    # Filter skewness
    df_sub <- subset(df, measure %in% c("ci_skewness", "data_skewness"))
    # Add sign and remove skewness == 0 because for cohen's kappa we need the
    # by two table and sign must be positive or negative
    rs <- within(
        df_sub,
        {
            sign <- sign(value)
            sign <- replace(sign, sign == 0, NA_real_)
        }
    )
    # Get everything into wide format
    rs <- tidyr::pivot_wider(
        data = rs,
        names_from = "measure",
        values_from = c("value", "sign"),
        values_fn = list
    )
    # calculate kappa and the correlation
    out <- vapply(
        seq_len(nrow(rs)),
        function(i) {
            with(
                rs[i, ],
                {
                    kap_ci <- sign_ci_skewness[[1L]]
                    kap_data <- sign_data_skewness[[1L]]
                    cor_ci <- value_ci_skewness[[1L]]
                    cor_data <- value_data_skewness[[1L]]
                    kappa <- if (any(is.na(c(kap_ci, kap_data)))) {
                        NA_real_ 
                    } else {
                        psych::cohen.kappa(cbind(kap_ci, kap_data))$kappa
                    }
                    correlation <- cor(cor_ci, cor_data)
                    c(kappa, correlation)
                }
            )
        },
        numeric(2L)
    )

    cbind(
        data.frame(
            value = c(out[1, ], out[2, ]),
            measure = rep(c("kappa", "correlation"), each = ncol(out))
        ),
        subset(
            rs,
            select = -c(
                sign_ci_skewness, value_ci_skewness,
                sign_data_skewness, value_data_skewness
            )
        )
    )
}


# This is a wrapper function that calls all of the other
# functions above.
get_summary_measures <- function(df, p_accept, pars) {
    # Get all of the statistics and put them into a data frame
    # Also add parameters
    mean_stats <- get_mean_stats(df = df)
    col_order <- names(mean_stats)
    gamma_stats <- get_gamma_stats(df = df, col_order = col_order)
    n_stats <- get_n_stats(df = df, col_order = col_order)
    bias_var_stats <- get_bias_var_stats(df = df, col_order = col_order)
    # This will be NULL if no publication bias
    p_accept_stats <- get_paccept_stats(
        p_accept = p_accept,
        col_order = col_order
    )

    # This works because in R, the return of rbind(df, NULL),
    # with df being a data.frame, is just df instead of an error
    res <- rbind(
        mean_stats,
        gamma_stats,
        n_stats,
        bias_var_stats,
        p_accept_stats,
        make.row.names = FALSE
    )
    # add the parameters and return
    p <- as.data.frame(lapply(pars, rep, times = nrow(res)))
    cbind(p, res)
}


## Wrapper that handles possible errors
calc_summary_measures <- function(df, p_accept, pars, i) {
    tryCatch({
        get_summary_measures(df = df, p_accept = p_accept, pars = pars)
    },
    error = function(cond) {
        error_function(
            cond = cond,
            pars = pars,
            error_obj = df,
            fun_name = "calc_summary_measures",
            i = i
        )
        NA
    })
}
