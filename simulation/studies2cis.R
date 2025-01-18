################################################################################
#               Conversion of individual studies to CIs                        #
################################################################################

# ------------------------------------------------------------------------------
#                   Helper function to calculate skewness from data
# ------------------------------------------------------------------------------

# Calculate the skewness from the estimates, SEs and the estimated tau2.
# This is essentially the same function as Weighted.Desc.Stat::w.skewness
calc_data_skewness <- function(thetahat, se, tau2){
    x <- thetahat
    mu <- 1 / (se^2 + tau2)
    mux <- mu * x
    summu <- sum(mu)
    wmean <- sum(mux) / summu
    wsd <- sqrt((sum(mux * x) / summu) - wmean^2)
    (sum(mu * (x - wmean)^3) / summu) / wsd^3
}


# ------------------------------------------------------------------------------
#                        Reference methods
# ------------------------------------------------------------------------------

## The following functions are used to fit a model to some individual studies
## and then extract the estimate and CIs. The idea is to make a function of the
## method, the individual studies, and their SEs. Then we can simply use this
## and pass only what we really need. All of these functions are only
## used within the main function called `get_classic_intervals()`. Thus, these
## functions contain the logic that ultimately return a data.frame with columns
## `lower`, `upper`, `method`, `ci_exists`.
## This main function is used in:
## - sim2CIs

## Fit model with Henmi & Copas
get_classic_obj_hc <- function(thetahat, se, tau, control) {
    metafor::hc(
        object = metafor::rma(
            yi = thetahat,
            sei = se,
            tau2 = tau^2,
            control = control
        )
    )
}
## Fit model with REML
get_classic_obj_reml <- function(thetahat, se, tau, control) {
    meta::metagen(
        TE = thetahat,
        seTE = se,
        sm = "MD",
        tau.preset = tau,
        # method.tau = "REML",
        control = control
    )
}
## Fit model with Hartung-Knapp
get_classic_obj_hk <- function(thetahat, se, tau, control) {
    meta::metagen(
        TE = thetahat,
        seTE = se,
        sm = "MD",
        tau.preset = tau,
        # method.tau = "REML",
        # hakn = TRUE,
        method.random.ci = "HK",
        control = control
    )
}
## Fit model with bayesmeta
tau_prior_bm <- function(x) {
    bayesmeta::dhalfnormal(x, scale = 0.3)
}
get_classic_obj_bm <- function(thetahat, se) {
    bayesmeta::bayesmeta(
        y = thetahat,
        sigma = se,
        tau.prior = tau_prior_bm
    )
}

## wrapper function (pass it a method ('hk', 'hc', ...)), some
## study estimates and standard errors and get the corresponding
## model object
get_classic_obj <- function(method, thetahat, se, tau, control) {
    switch(
        method,
        "hk" = get_classic_obj_hk(
            thetahat = thetahat,
            se = se,
            tau = tau,
            control = control
        ),
        "hc" = get_classic_obj_hc(
            thetahat = thetahat,
            se = se,
            tau = tau,
            control = control
        ),
        "reml" = get_classic_obj_reml(
            thetahat = thetahat,
            se = se,
            tau = tau,
            control = control
        ),
        "bm" = get_classic_obj_bm(
            thetahat = thetahat,
            se = se
        )
    )
}


## extract lower and upper bound for Hartung-Knapp and REML CI
## from the model object
get_classic_ci_reml <- get_classic_ci_hk <- function(obj) {
    setNames(
        with(obj, c(lower.random, upper.random, TE.random)),
        c("lower", "upper", "estimate")
    )
}
## extract lower and upper bound for Hartung-Knapp and REML PI
## from the model object
get_classic_pi_reml <- get_classic_pi_hk <- function(obj) {
    with(obj, c(lower.predict, upper.predict))
}
## extract lower and upper bound for Henmi & Copas CI
## from the model object
get_classic_ci_hc <- function(obj) {
    setNames(with(obj, c(ci.lb, ci.ub, beta)), c("lower", "upper", "estimate"))
}
## extract lower and upper bound for Bayesmeta
## from the model object
get_classic_ci_bm <- function(obj) {
    num <- unname(obj$summary[c(1L, 5:6), 3L])
    setNames(num[c(2, 3, 1)], c("lower", "upper", "estimate"))
}
## which function to call depending on the method
get_classic_interval <- function(method, obj) {
    switch(
        method,
        "hk_pi" = get_classic_pi_hk(obj = obj),
        "reml_pi" = get_classic_pi_reml(obj = obj),
        "hk_ci" = get_classic_ci_hk(obj = obj),
        "reml_ci" = get_classic_ci_reml(obj = obj),
        "hc_ci" = get_classic_ci_hc(obj = obj),
        "bm_ci" = get_classic_ci_bm(obj = obj)
    )
}

## This is the final wrapper function of all the above functions
## related to the calculation of CIs for the classic methods, i.e.
## all methods that use the meta, metafor, or bayesmeta package.
## This function is used in:
## - sim2CIs
get_classic_intervals <- function(methods, thetahat, se, tau2, control) {
    # set control options
    # control <- list(maxiter = 1e4, stepadj = 0.25)
    # mappings for method codes to object codes
    # as some methods use the same objects
    # (prediction & confidence intervals are in the same object)
    meth2obj_code <- c(
        "hc_ci" = "hc",
        "reml_ci" = "reml",
        "reml_pi" = "reml",
        "hk_ci" = "hk",
        "hk_pi" = "hk",
        "bm_ci" = "bm"
    )
    # which objects to fit
    obj_codes <- meth2obj_code[methods]
    obj_fit <- unique(obj_codes)
    # fit objects
    objs <- lapply(
        obj_fit,
        get_classic_obj,
        thetahat = thetahat,
        se = se,
        tau = sqrt(tau2),
        control = control
    )
    names(objs) <- obj_fit
    # get the CIs/PIs and calculate the data skewness
    ci <- matrix(NA_real_, nrow = length(methods), ncol = 3L)
    for (i in seq_along(methods)) {
        ci[i, ] <- get_classic_interval(
            method = methods[i],
            obj = objs[[obj_codes[i]]]
        )
    }
    # construct names for methods
    tau_est_method <- names(tau2)
    est_method <- function(obj_codes, tau_est_method) {
        get1 <- function(obj_code, tau2) {
            if (!(obj_code %in% c("hc", "bm"))) {
            # if (!(obj_code %in% c("bm"))) {
                tau_est_method
            } else if (obj_code == "bm") {
                "HN(0, 0.3)"
            } else if (obj_code == "hc") {
                "DL"
            }
        }
        vapply(obj_codes, get1, character(1L))
    }
    mt <- est_method(obj_codes = obj_codes, tau_est_method = tau_est_method)
    meth <- paste0(
        methods,
        "_additive_", # since all of these are additive
        mt
    )

    # convert this to df
    out <- data.frame(
        lower = ci[, 1L],
        upper = ci[, 2L],
        estimate = ci[, 3L],
        method = meth,
        ci_exists = TRUE,
        is_ci = grepl("_ci$", methods),
        is_pi = grepl("_pi$", methods), ## not used anymore
        is_new = FALSE,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    out
}


# ------------------------------------------------------------------------------
#                               New methods
# ------------------------------------------------------------------------------

## The following functions do something quite similar to the ones above,
## i.e. they are helpers called within get_new_intervals(), which itself
## is called in sim2CIs(). Thus, these functions contain the logic that is
## needed to ultimately return a data.frame with columns `lower`, `upper`,
## `method`, and `ci_exists`.
## These functions are used in:
## - sim2CIs

## list all  the p-value functions here:
## This should be equivalent to all the p_* functions in the confMeta package
## Also, note that the names of the list elements must be consistent with the
## names in the `sim2CIs` function
get_p_value_functions <- function(methods) {
    mt <- unique(sub("_.*$", "", methods))
    p_val_fcts <- list(
        # hMeanF = hMeanChiSqMu_f,
        # hMeanChisq = hMeanChiSqMu_chisq,
        edgington = confMeta::p_edgington,
        wilkinson = confMeta::p_wilkinson,
        pearson = confMeta::p_pearson,
        tippett = confMeta::p_tippett,
        fisher = confMeta::p_fisher
    )
    p_val_fcts[mt]
}

## list all the argument configurations here
## This relates to the type of heterogeneity, used to adjust the
## standard errors of the individual studies. These arguments are
## subsequently passed to the p-value functions in `get_p_value_functions()`
get_p_value_args <- function(tau2, phi) {
    list(
        # none = list(
        #     heterogeneity = "none",
        #     check_inputs = FALSE
        # ),
        additive = list(
            heterogeneity = "additive",
            tau2 = tau2,
            check_inputs = FALSE
        )#,
        # multiplicative = list(
        #     heterogeneity = "multiplicative",
        #     phi = phi,
        #     check_inputs = FALSE
        # )
    )
}


# ## This returns the x and y coords of the p-value function
# ## between the smallest effect in the CI and the largest
# ## effect in the CI.
# get_gamma_x_y <- function(gamma) {
#     gamma_exists <- all(is.na(gamma))
#     if (gamma_exists) {
#         c(x = NA_real_, y = NA_real_)
#     } else {
#         y <- gamma[, 2L]
#         x <- gamma[, 1L]
#         min_idx <- which.min(y)
#         c(x = x[min_idx], y = y[min_idx])
#     }
# }

# loop over each of the function/argument combination
## The function creates a grid of p-value functions and argument configurations,
## i.e. all combinations of the outputs of `get_p_value_functions()` and
## `get_p_value_args()`. It then loops over the grid and calculates the CI
## for each p-value function/argument combination.
get_new_ci_gamma <- function(thetahat, se, p_funcs, arguments, methods) {

    # Construct a grid from p-value-functions and argument
    # configurations
    g <- expand.grid(
        p_funcs_idx = seq_along(p_funcs),
        arguments_idx = seq_along(arguments)
    )
    # Extract columns since looping over
    # atomic vectors is faster
    p_funcs_idx <- g$p_funcs_idx
    p_arguments_idx <- g$arguments_idx

    # Get the number of FUN/arg combinations
    l <- nrow(g)

    # make names from the p_value function and the heterogeneity
    nms <- paste0(
        names(p_funcs)[p_funcs_idx],
        "_",
        names(arguments)[p_arguments_idx]
    )

    # Pre-allocate memory for output
    ci <- vector("list", l)
    names(ci) <- nms
    # gamma <- matrix(NA_real_, ncol = 2L, nrow = l)

    # compute CIs and other stuff
    for (k in seq_len(l)) {
        arg_idx <- p_arguments_idx[k]
        f_idx <- p_funcs_idx[k]
        # make the function for current config
        arg_list <- arguments[[arg_idx]]
        f <- p_funcs[[f_idx]]
        formals(f) <- modifyList(formals(f), arg_list)
        # make names
        function_name <- names(p_funcs)[f_idx]
        meth <- grep(paste0(function_name, "_.*$"), methods, value = TRUE)
        times <- length(meth)
        # compute CI
        res <- confMeta:::get_ci(
            estimates = thetahat,
            SEs = se,
            conf_level = 0.95,
            p_fun = f
        )
        mt <- paste0(
            meth,
            "_",
            arg_list$heterogeneity,
            "_",
            if (is.null(arg_list$tau2)) {
                NA_character_
            } else {
                names(arg_list$tau2)
            }
        )
        # CI
        CI <- data.frame(
            lower = rep(res$CI[, 1L], times = times),
            upper = rep(res$CI[, 2L], times = times),
            method = rep(mt, each = nrow(res$CI)),
            ci_exists = rep(
                ifelse(is.na(res$CI[, 1L]), FALSE, TRUE),
                times = times
            ),
            is_ci = NA_character_,
            is_pi = NA_character_,
            is_new = TRUE,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        CI$is_pi <- grepl("_pi_", CI$method)
        CI$is_ci <- grepl("_ci_", CI$method)
        # estimates
        es <- res$p_max[, 1L]
        estimates <- data.frame(
            estimate = rep(es, times = times),
            method = rep(mt, each = length(es)),
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        # p_max
        pm <- unique(res$p_max[, 2])
        p_max <- data.frame(
            p_max = rep(pm, times = times),
            method = rep(mt, each = length(pm)),
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        # AUCC and AUCC ratio
        aucc <- data.frame(aucc = res$aucc, method = mt)
        aucc_ratio <- data.frame(aucc_ratio = res$aucc_ratio, method = mt)
        # return list
        ci[[k]] <- list(
            "CI" = CI,
            "estimates" = estimates,
            "p_max" = p_max,
            "aucc" = aucc,
            "aucc_ratio" = aucc_ratio
        )
        # # gamma related stuff
        # gamma[k, ] <- get_gamma_x_y(res$gamma)
    }

    # rbind all of the list elements
    out <- lapply(
        c(
            "CI" = "CI",
            "estimates" = "estimates",
            "p_max" = "p_max",
            "aucc" = "aucc",
            "aucc_ratio" = "aucc_ratio"
        ),
        function(ci, df) {
            do.call(
                "rbind",
                append(lapply(ci, "[[", i = df), list(make.row.names = FALSE))
            )
        },
        ci = ci
    )

    out
}

## This function is the main function we use in order to
## calculate the CIs as well as the minima and maxima of the
## p-value function. It is a wrapper of all of the above
get_new_intervals <- function(
    methods,
    thetahat,
    se,
    phi,
    tau2
) {

    # get a list of p-value functions
    p_funcs <- get_p_value_functions(methods = methods)
    # get a list of arguments
    arguments <- get_p_value_args(phi = phi, tau2 = tau2)
    # get the cis
    get_new_ci_gamma(
        thetahat = thetahat,
        se = se,
        p_funcs = p_funcs,
        arguments = arguments,
        methods = methods
    )
}

################################################################################
#                            Calculating the CIs                               #
################################################################################

# estimate tau2 with different methods
get_tau2 <- function(thetahat, se, methods, control) {
    get_one_tau2 <- function(thetahat, se, method, control) {
        if (method == "none") {
            0
        } else {
            meta::metagen(
                TE = thetahat,
                seTE = se,
                method.tau = method,
                sm = "MD",
                control = control
            )$tau2
        }
    }
    vapply(
        methods,
        get_one_tau2,
        se = se,
        thetahat = thetahat,
        control = control,
        FUN.VALUE = double(1L)
    )
}

#' Confidence intervals from effect estimates and their standard errors
#'
#' Takes the output of \code{simRE} and returns CIs for
#' combined effect using the
#' indicated methods.
#' @param x matrix output from \code{simRE}.
#' @return a tibble with columns \code{lower}, \code{upper}, and \code{method}.
sim2CIs <- function(x) {

    # Get the thetahats and the SEs
    thetahat <- x[, 1L]
    se <- x[, 2L]
    delta <- x[, 3L]

    # Store the attributes of x
    att <- attributes(x)

    control <- list(
        stepadj = 0.25,
        maxiter = 1e4,
        tau2.max = 1e4
    )

    # estimate tau2 with these methods
    tau2_all <- get_tau2(
        # methods = c("none", "DL", "PM", "REML"),
        methods = c("none", "DL", "REML"),
        thetahat = thetahat,
        se = se,
        control = control
    )

    # get the CIs of the classic methods using the
    classic_methods <- setNames(
        lapply(
            seq_along(tau2_all),
            function(thetahat, se, tau2_all, i, control) {
                tau2 <- tau2_all[i]
                # Henmi-Copas is only implemented for DL
                methods <- if (names(tau2) == "DL") {
                    c("hk_ci", "hc_ci", "reml_ci")
                } else {
                    c("hk_ci", "reml_ci")
                }
                # get those from classic methods
                out_ints <- get_classic_intervals(
                    methods = methods,
                    thetahat = thetahat,
                    se = se,
                    tau2 = tau2,
                    control = control
                )
                out_ints$skewness_data <- calc_data_skewness(
                    thetahat = thetahat,
                    se = se,
                    tau2 = tau2
                )
                out_ints
            },
            thetahat = thetahat,
            se = se,
            tau2_all = tau2_all,
            control = control
        ),
        names(tau2_all)
    )
    classic_methods <- unique(do.call("rbind", classic_methods))

    # For the new methods, we need phi and tau2 to estimate heterogeneity
    # we estimate phi ourselves and for tau2 we check whether it has already
    # been computed. If not, estimate it ourselves
    ## change Leo April 10 2024
    phi <- confMeta::estimate_phi(estimates = thetahat, SEs = se)

    # Loop over tau2 estimation methods ["none", "DL", "PM", "REML"] and
    # calculate [CI, gamma, estimates, p_max]
    new_methods <- setNames(
        lapply(
            seq_along(tau2_all),
            function(thetahat, se, tau2_all, i, control) {
                tau2 <- tau2_all[i]
                # get those from new methods
                out_ints <- get_new_intervals(
                    methods = c(
                        ## For the sake of naming things:
                        ## -> do NOT use underscores("_") in the
                        ##    function names here!
                        ##    The underscores are only there to
                        ##    separate the function names as defined in
                        ##    `get_p_value_functions()` from the desired interval
                        ##    which can be either `ci` or `pi`.

                        # "hMeanF_ci",
                        # "hMeanChisq_ci",
                        "tippett_ci",
                        "wilkinson_ci",
                        "pearson_ci",
                        "edgington_ci",
                        "fisher_ci"
                    ),
                    thetahat = thetahat,
                    se = se,
                    # phi = phi,
                    phi = NULL,
                    tau2 = tau2
                )
                out_ints$CI$skewness_data <- calc_data_skewness(
                    thetahat = thetahat,
                    se = se,
                    tau2 = tau2
                )
                out_ints
            },
            thetahat = thetahat,
            se = se,
            tau2_all = tau2_all,
            control = control
        ),
        names(tau2_all)
    )
    new_methods_ci <- do.call(
        what = "rbind",
        args = lapply(new_methods, "[[", i = "CI")
    )
    new_methods_est <- do.call(
        what = "rbind",
        args = lapply(new_methods, "[[", i = "estimates")
    )

    # Assemble output
    CIs <- rbind(
        subset(classic_methods, select = -c(estimate)),
        new_methods_ci,
        make.row.names = FALSE
    )
    estimates <- rbind(
        subset(classic_methods, select = c(estimate, method)),
        new_methods_est,
        # do.call("rbind", lapply(new_methods, "[[", i = "estimates")),
        make.row.names = FALSE
    )
    p_max <- do.call(
        "rbind",
        append(
            lapply(new_methods, "[[", i = "p_max"),
            list(make.row.names = FALSE)
        )
    )
    aucc <- do.call(
        "rbind",
        append(
            lapply(new_methods, "[[", i = "aucc"),
            list(make.row.names = FALSE)
        )
    )
    aucc_ratio <- do.call(
        "rbind",
        append(
            lapply(new_methods, "[[", i = "aucc_ratio"),
            list(make.row.names = FALSE)
        )
    )

    list(
        CIs = CIs,
        estimates = estimates,
        p_max = p_max,
        aucc = aucc,
        aucc_ratio = aucc_ratio,
        model = att$heterogeneity, # this is the simulation model
        theta = thetahat,
        delta = delta,
        effect = att$effect,
        median = att$median
    )
}


## Wrapper function around sim2CIs() that, in case of errors, runs
## the error_function such that the simulation stops immediately
## instead of continuing to calculate all the other scenarios
calc_ci <- function(x, pars, i) {
    # if error in function before, just return NA back.
    # Otherwise, run sim2CIs. If an error happens, log everything and
    # return NA
    if (length(x) == 1L && is.na(x)) {
        NA
    } else {
        tryCatch(
            {
                sim2CIs(x = x)
            },
            error = function(cond) {
                error_function(
                    cond = cond,
                    pars = pars,
                    error_obj = x,
                    fun_name = "calc_ci",
                    i = i
                )
                NA
            }
        )
    }
}
