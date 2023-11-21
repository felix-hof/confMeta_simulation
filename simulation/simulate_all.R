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
# rm(list = ls())
remotes::install_github("felix-hof/confMeta")
library(confMeta)
library(doParallel)
library(doRNG)
library(RhpcBLASctl)
blas_set_num_threads(1) # multi threading of BLAS

################################################################################
#                              Helper functions                                #
################################################################################


# # Helper function to convert the id_strings into proper names
# # This function is used in:
# # - get_new_intervals
# make_name <- function(id_strings) {
#     str <- strsplit(id_strings, split = "_")
#     nms <- vapply(
#         str,
#         function(x) {
#             method <- switch(
#                 x[1L],
#                 "hMean" = "Harmonic Mean",
#                 "ktrials" = "k-Trials",
#                 "pearson" = "Pearson",
#                 "edgington" = "Edgington",
#                 "fisher" = "Fisher"
#             )
#             het <- switch(
#                 x[length(x)],
#                 "none" = "",
#                 "additive" = "Additive",
#                 "multiplicative" = "Multiplicative"
#             )
#             distr <- if (length(x) == 3L) {
#                 switch(
#                     x[2L],
#                     "f" = "(f)",
#                     "chisq" = "(chisq)"
#                 )
#             } else {
#                 NA_character_
#             }
#             out <- paste(
#                 method, het, "CI"
#             )
#             if (!is.na(distr)) out <- paste(out, distr)
#             out
#         },
#         character(1L)
#     )
#     gsub("\\s+", " ", nms)
# }

## This function repeats each element of `x` exactly `each`
## times. However, in contrast to the regular `rep()` function
## `each` can also be vector of the same length as 'x'.
## repeat each element of `x` `each` times
## rep2(c("a", "b", "c"), c(1, 2, 3)) --> c("a", "b", "b", "c", "c", "c")
rep2 <- function(x, each) {
    do.call(
        "c",
        mapply(FUN = rep, x = x, each = each, SIMPLIFY = FALSE)
    )
}

## This function is only used in tryCatch() blocks in case of an error.
## The function writes the name of the function that errored, the error
## message, and the parameters to a file on disk. It also saves an object
## to disk containing the function inputs. It is mainly intended for
## collecting logs in case of errors
## This function is used in:
## - sim_effects
error_function <- function(cond, pars, error_obj = NULL, fun_name, i) {
    text <- capture.output(cond)
    out_msg <- paste0(
        "Error in ", fun_name, " iteration: ", i, "\n\n",
        "Parameters are:\n\n",
        paste0(paste0(names(pars), ": ", pars[1, ]), collapse = "\n"),
        "\n\nThe error message is:\n", text, "\n\n",
        if (is.null(error_obj)) {
            "See pars.rds for the saved parameters."
        } else {
            paste0(
                "See error.rds for the last available object",
                " and pars.rds for the saved parameters."
            )
        },
        "\n\n\n",
        "---------------------------------------------------------------------",
        "-----------\n\n\n"
    )
    cat(out_msg, file = "error.txt", append = TRUE)
    saveRDS(pars, file = "pars.rds")
    if (!is.null(error_obj)) saveRDS(error_obj, file = "error.rds")
    return(NA)
}

## These functions just call hMeanChiSqMu under the hood but fix the argument
## `distr` to be either `"f"` or `"chisq"`. This is useful because we do not
## actually need to worry about including `distr = "f"` in the `args` argument
## when constructing expressions of the form `do.call("hMeanChiSqMu", args)`.
## These functions are used in:
## - sim2CIs
hMeanChiSqMu_f <- function(
    estimates,
    SEs,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "none",
    check_inputs = TRUE,
    w = rep(1, length(estimates))
) {
    confMeta::p_hmean(
        estimates = estimates,
        SEs = SEs,
        mu = mu,
        phi = phi,
        tau2 = tau2,
        heterogeneity = heterogeneity,
        alternative = alternative,
        check_inputs = check_inputs,
        w = w,
        distr = "f"
    )
}

hMeanChiSqMu_chisq <- function(
    estimates,
    SEs,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "none",
    check_inputs = TRUE,
    w = rep(1, length(estimates))
) {
    confMeta::p_hmean(
        estimates = estimates,
        SEs = SEs,
        mu = mu,
        phi = phi,
        tau2 = tau2,
        heterogeneity = heterogeneity,
        alternative = alternative,
        check_inputs = check_inputs,
        w = w,
        distr = "chisq"
    )
}

## The following functions are used to fit a model to some individual studies
## and then extract the estimate and CIs. The idea is to make a function of the
## method, the individual studies, and their SEs. Then we can simply use this
## and pass only what we really need. Most of these function function are only
## used within the main function called get_classic_intervals(). Thus, these
## functions contain the logic that ultimately return a data.frame with columns
## `lower`, `upper`, `method`, `ci_exists`.
## This main function is used in:
## - sim2CIs

## Fit model with Henmi & Copas
get_classic_obj_hc <- function(thetahat, se, ...) {
    metafor::hc(
        object = metafor::rma(yi = thetahat, sei = se, ...)
    )
}
## Fit model with REML
get_classic_obj_reml <- function(thetahat, se, ...) {
    meta::metagen(
        TE = thetahat,
        seTE = se,
        sm = "MD",
        method.tau = "REML",
        ...
    )
}
## Fit model with Hartung-Knapp
get_classic_obj_hk <- function(thetahat, se, ...) {
    meta::metagen(
        TE = thetahat,
        seTE = se,
        sm = "MD",
        method.tau = "REML",
        hakn = TRUE,
        ...
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

## wrapper function
get_classic_obj <- function(method, thetahat, se, ...) {
    switch(
        method,
        "hk" = get_classic_obj_hk(
            thetahat = thetahat,
            se = se,
            ... = ...
        ),
        "hc" = get_classic_obj_hc(
            thetahat = thetahat,
            se = se,
            ... = ...
        ),
        "reml" = get_classic_obj_reml(
            thetahat = thetahat,
            se = se,
            ... = ...
        ),
        "bm" = get_classic_obj_bm(
            thetahat = thetahat,
            se = se
        )
    )
}
## extract lower and upper bound for Hartung-Knapp and REML CI
get_classic_ci_reml <- get_classic_ci_hk <- function(obj) {
    setNames(
        with(obj, c(lower.random, upper.random, TE.random)),
        c("lower", "upper", "estimate")
    )
}
## extract lower and upper bound for Hartung-Knapp and REML PI
get_classic_pi_reml <- get_classic_pi_hk <- function(obj) {
    with(obj, c(lower.predict, upper.predict))
}
## extract lower and upper bound for Henmi & Copas CI
get_classic_ci_hc <- function(obj) {
    setNames(with(obj, c(ci.lb, ci.ub, beta)), c("lower", "upper", "estimate"))
}
## extract lower and upper bound for Bayesmeta
get_classic_ci_bm <- function(obj) {
    num <- unname(obj$summary[c(1L, 5:6), 3L])
    setNames(num[c(2, 3, 1)], c("lower", "upper", "estimate"))
}
## which function to call depending on the method
get_classic_interval <- function(method, obj) {
    switch(
        method,
        "hk_pi" = get_classic_pi_hk(obj = obj),     # unusable since no estimate
        "reml_pi" = get_classic_pi_reml(obj = obj), # unusable since no estimate
        "hk_ci" = get_classic_ci_hk(obj = obj),
        "reml_ci" = get_classic_ci_reml(obj = obj),
        "hc_ci" = get_classic_ci_hc(obj = obj),
        "bm_ci" = get_classic_ci_bm(obj = obj)
    )
}

## This is the final wrapper function of all the above functions
## related to the calculation of CIs for the classic methods, i.e.
## all methods that use the meta package.
## It also includes a part for bayesmeta
## This function is used in:
## - sim2CIs
get_classic_intervals <- function(methods, thetahat, se) {
    # set control options
    control <- list(maxiter = 1e4, stepadj = 0.25)
    # hash tables for methods2obj and methods2int
    # as some methods use the same objects
    meth2obj_code <- c(
        "hc_ci" = "hc",
        "reml_ci" = "reml",
        "reml_pi" = "reml",
        "hk_ci" = "hk",
        "hk_pi" = "hk",
        "bm_ci" = "bm"
    )
    # meth2name <- c(
    #     "hc_ci" = "Henmi & Copas CI",
    #     "reml_ci" = "REML CI",
    #     "reml_pi" = "REML PI",
    #     "hk_ci" = "Hartung & Knapp CI",
    #     "hk_pi" = "Hartung & Knapp PI",
    #     "bm_ci" = "Bayesmeta CI"
    # )
    # which objects to fit
    obj_codes <- meth2obj_code[methods]
    obj_fit <- unique(obj_codes)
    # fit objects
    objs <- lapply(
        obj_fit,
        get_classic_obj,
        thetahat = thetahat,
        se = se,
        control = control
    )
    names(objs) <- obj_fit
    # get the CIs
    ci <- matrix(NA_real_, nrow = length(methods), ncol = 3L)
    for (i in seq_along(methods)) {
        ci[i, ] <- get_classic_interval(
            method = methods[i],
            obj = objs[[obj_codes[i]]]
        )
    }
    # construct names for methods
    est_method <- function(obj_codes) {
        get1 <- function(obj_code) {
            switch(
                obj_code,
                "hk" = "REML",
                "hc" = "DL",
                "bm" = "HN(0, 0.3)",
            )
        }
        vapply(obj_codes, get1, character(1L))
    }
    mt <- est_method(obj_codes = obj_codes)
    meth <- paste0(
        methods,
        "_additive_", # since all of these are additive
        mt
    )

    # convert this to df
    # which methods to call
    out <- data.frame(
        lower = ci[, 1L],
        upper = ci[, 2L],
        estimate = ci[, 3L],
        method = meth,
        ci_exists = TRUE,
        is_ci = grepl("_ci$", methods),
        is_pi = grepl("_pi$", methods), ## no need for PIs
        is_new = FALSE,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    # attach attribute tau2
    # if ("reml" %in% obj_codes) {
    #     attr(out, which = "tau2") <- objs$reml$tau2
    # }
    # return
    out
}

## The following functions do something quite similar to the ones above,
## i.e. they are helpers called within get_new_intervals(), which itself
## is called in sim2CIs(). Thus, these functions contain the logic that is
## needed to ultimately return a data.frame with columns `lower`, `upper`,
## `method`, and `ci_exists`.
## These functions are used in:
## - sim2CIs

## list all  the p-value functions here
get_p_value_functions <- function(methods) {
    mt <- unique(sub("_.*$", "", methods))
    p_val_fcts <- list(
        hMeanF = hMeanChiSqMu_f,
        hMeanChisq = hMeanChiSqMu_chisq,
        ktrials = confMeta::p_ktrials,
        pearson = confMeta::p_pearson,
        edgington = confMeta::p_edgington,
        fisher = confMeta::p_fisher
    )
    p_val_fcts[mt]
}

## list all the argument configurations here
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


## This returns the x and y coords of the p-value function
## between the smallest effect in the CI and the largest
## effect in the CI.
get_gamma_x_y <- function(gamma) {
    gamma_exists <- all(is.na(gamma))
    if (gamma_exists) {
        c(x = NA_real_, y = NA_real_)
    } else {
        y <- gamma[, 2L]
        x <- gamma[, 1L]
        min_idx <- which.min(y)
        c(x = x[min_idx], y = y[min_idx])
    }
}

# loop over each of the function/argument combination
## calculate the CIs, return a data.frame
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
        # return list
        ci[[k]] <- list(
            "CI" = CI,
            "estimates" = estimates,
            "p_max" = p_max
        )
        # # gamma related stuff
        # gamma[k, ] <- get_gamma_x_y(res$gamma)
    }

    # rbind all of the list elements
    out <- lapply(
        c("CI" = "CI", "estimates" = "estimates", "p_max" = "p_max"),
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
## calculate the CI as well as the minima and maxima of the
## p-value function.
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
#                    Simulating effects and standard errors                    #
################################################################################

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
simRE <- function(
    k,
    sampleSize,
    effect,
    I2,
    heterogeneity,
    dist,
    large
) {

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
            delta <- sn::rst(n = k, xi = effect, omega = sqrt(tau2 / 2), nu = 4)
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
            delta <- sn::rst(n = k, xi = effect, omega = sqrt(tau2 / 2), nu = 4)
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
simREbias <- function(
    k,
    sampleSize,
    effect,
    I2,
    heterogeneity = c("additive", "multiplicative"),
    dist = c("t", "Gaussian"),
    large,
    bias = c("none", "moderate", "strong"),
    verbose = TRUE,
    check_inputs = TRUE
) {
    # input checks
    if (check_inputs) {
        stopifnot(
            length(k) == 1L,
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
            k >= large
        )
    }

    bias <- match.arg(bias)
    dist <- match.arg(dist)
    heterogeneity <- match.arg(heterogeneity)

    if (bias == "none") {
        o <- simRE(
            k = k, sampleSize = sampleSize, effect = effect, I2 = I2,
            heterogeneity = heterogeneity, dist = dist, large = large
        )
        ## add attributes and return
        attr(o, "heterogeneity") <- heterogeneity
        attr(o, which = "effect") <- effect
        attr(o, which = "p_accept") <- 1
        return(o)
    }

    ## first ignore the 'large'
    o <- simRE(
        k = k * 3, sampleSize = sampleSize, effect = effect, I2 = I2,
        heterogeneity = heterogeneity, dist = dist, large = 0
    )
    pa <- pAccept(theta = o[, "theta"], se = o[, "se"], bias = bias)
    keep <- rbinom(n = k * 3, size = 1, prob = pa)
    while (k > sum(keep)) {
        if (verbose) cat(".")
        o2 <- simRE(
            k = k * 3, sampleSize = sampleSize, effect = effect,
            I2 = I2, heterogeneity = heterogeneity, dist = dist,  large = 0
        )
        pa2 <- pAccept(theta = o2[, "theta"], se = o2[, "se"], bias = bias)
        keep2 <- rbinom(n = k * 3, size = 1, p = pa2)
        o <- rbind(o, o2)
        keep <- c(keep, keep2)
        pa <- c(pa, pa2)
    }
    o <- o[as.logical(keep), ][1:k, ]

    ## add large studies
    if (large != 0) {
        oLarge <- simRE(
            k = large, sampleSize = sampleSize, effect = effect,
            I2 = I2, heterogeneity = heterogeneity, dist = dist, large = large
        )
        o <- rbind(oLarge, o[-seq_len(large), ])
    }

    ## add attributes and return
    attr(o, which = "heterogeneity") <- heterogeneity
    attr(o, which = "effect") <- effect
    attr(o, which = "p_accept") <- mean(pa)
    o
}

## Wrapper function around simREbias() that, in case of errors, runs
## the error_function such that the simulation stops immediately
## instead of continuing to calculate all the other scenarios
sim_effects <- function(pars, i) {
    # run simREbias on the elements of a list/dataframe and return
    # if there is an error, call error_function
    tryCatch(
        {
            with(
                pars,
                simREbias(
                    k = k,
                    sampleSize = sampleSize,
                    effect = effect,
                    I2 = I2,
                    heterogeneity = heterogeneity,
                    dist = dist,
                    large = large,
                    bias = bias,
                    verbose = FALSE,
                    check_inputs = FALSE
                )
            )
        },
        error = function(cond) {
            error_function(
                cond = cond,
                pars = pars,
                fun_name = "sim_effects",
                i = i
            )
            NA
        }
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

get_ints_reml <- function(thetahat, se, tau2, control) {
    obj <- meta::metagen(
        TE = thetahat,
        seTE = se,
        sm = "MD",
        tau.preset = sqrt(tau2),
        control = control
    )
    ci <- get_classic_interval(method = "reml_ci", obj = obj)
    mt <- paste0("reml_ci_additive_", names(tau2))
    ci_df <- data.frame(
        lower = ci[1L],
        upper = ci[2L],
        method = mt,
        ci_exists = TRUE,
        is_ci = TRUE,
        is_pi = FALSE,
        is_new = FALSE,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    est_df <- data.frame(
        estimate = ci[3L],
        method = mt,
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    list("CI" = ci_df, "estimates" = est_df)
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

    # estimate tau2
    tau2_all <- get_tau2(
        methods = c("none", "DL", "PM", "REML"),
        thetahat = thetahat,
        se = se,
        control = control
    )

    # get the CIs of the classic methods
    classic_methods <- get_classic_intervals(
        methods = c("hk_ci", "hc_ci"), # "bm_ci"),
        thetahat = thetahat,
        se = se
    )

    # For the new methods, we need phi and tau2 to estimate heterogeneity
    # we estimate phi ourselves and for tau2 we check whether it has already
    # been computed. If not, estimate it ourselves
    # phi <- confMeta::estimate_phi(estimates = thetahat, SEs = se)
    # tau2 <- attributes(classic_methods)$tau2
    # if (is.null(tau2)) {
    #     tau2 <- confMeta::estimate_tau2(estimates = thetahat, SEs = se)
    # }

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
                        ## -> do NOT use underscores("_") in these
                        ##    names here!

                        # "hMeanF_ci",
                        # "hMeanChisq_ci",
                        # "ktrials_ci",
                        # "pearson_ci",
                        "edgington_ci",
                        "fisher_ci"
                    ),
                    thetahat = thetahat,
                    se = se,
                    # phi = phi,
                    phi = NULL,
                    tau2 = tau2
                )
                ints_reml <- get_ints_reml(
                    thetahat = thetahat,
                    se = se,
                    tau2 = tau2,
                    control = control
                )
                out_ints$CI <- rbind(out_ints$CI, ints_reml$CI)
                out_ints$estimates <- rbind(
                    out_ints$estimates,
                    ints_reml$estimates
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

    # Assemble output
    CIs <- rbind(
        subset(classic_methods, select = -c(estimate)),
        do.call("rbind", lapply(new_methods, "[[", i = "CI")),
        make.row.names = FALSE
    )
    estimates <- rbind(
        subset(classic_methods, select = c(estimate, method)),
        do.call("rbind", lapply(new_methods, "[[", i = "estimates")),
        make.row.names = FALSE
    )
    p_max <- do.call(
        "rbind",
        append(
            lapply(new_methods, "[[", i = "p_max"),
            list(make.row.names = FALSE)
        )
    )

    list(
        CIs = CIs,
        estimates = estimates,
        p_max = p_max,
        model = att$heterogeneity, # this is the simulation model
        theta = thetahat,
        delta = delta,
        effect = att$effect
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



################################################################################
#            Calculating the output measures for each individual CI            #
################################################################################

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
        coverage_true = quote({
            calc_coverage_true(
                cis = cis,
                effect = effect,
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
            "coverage_true",
            # "coverage_effects",
            # "coverage_effects_min1",
            # "coverage_effects_all",
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
            "p_max"
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
    # get the deltas
    delta <- x$delta
    # get the effect
    effect <- x$effect

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
            effect = effect,
            ci_exists = ci_exists[k],
            delta = delta,
            pars = pars,
            p_max = if (curr_is_new) {
                p_max[p_max_method == curr_method]
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

################################################################################
#            Calculating the output measures for each individual CI            #
################################################################################

# Calculate the mean measure for each of the method and measure subgroups
get_mean_stats <- function(df) {
    # Summarise all the measures except gamma_min & p_max
    df_sub <- subset(
        df,
        !(measure %in% c("gamma_min", "p_max", "mse", "estimate"))
    )
    mean_stats <- stats::aggregate(
        value ~ measure + method + is_ci + is_pi + is_new,
        FUN = mean,
        data = df_sub
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
    f_list_gamma <- list(
        min = min,
        firstQuart = function(x) {
            quantile(x, probs = 0.25, names = FALSE)
        },
        median = median,
        mean = mean,
        thirdQuart = function(x) {
            quantile(x, probs = 0.75, names = FALSE)
        },
        max = max
    )
    gamma_stats <- lapply(
        seq_along(f_list_gamma),
        function(i) {
            res <- stats::aggregate(
                value ~ measure + method + is_ci + is_pi + is_new,
                FUN = f_list_gamma[[i]],
                data = df_sub
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
make_f_n <- function(n) {
    # This function generates a list of functions
    # that check how many times the number of
    # intervals is equal to n. Here, n is vectorised
    # and there is also a function that counts the number
    # of times where the number of intervals exceeds
    # max(n).
    maxn <- max(n)
    f_n <- lapply(
        n,
        function(one_n) {
            # force(one_n)
            function(x) sum(x == one_n)
        }
    )
    f_n <- append(f_n, list(function(x) sum(x > maxn)))
    names(f_n) <- c(as.character(n), paste0("gt", max(n)))
    f_n
}

get_n_stats <- function(df, col_order) {
    f_list_n <- make_f_n(0:9)
    df_sub <- subset(df, measure == "n")
    n_stats <- lapply(
        seq_along(f_list_n),
        function(i, f_list_n, df_sub) {
            res <- stats::aggregate(
                value ~ measure + method + is_ci + is_pi + is_new,
                FUN = f_list_n[[i]],
                data = df_sub
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
        FUN = function(x) mean(x, na.rm = TRUE),
        data = df_sub_est
    )
    bias <- within(mean_est, value <- abs(value - effect))
    # Calculate variance, i.e. var(estimates)
    var_est <- stats::aggregate(
        value ~ method,
        FUN = function(x) var(x, na.rm = TRUE),
        data = df_sub_est
    )
    # Calculate MSE, i.e. (1/n) * sum((estimate - true_effect)^2)
    mse_est <- stats::aggregate(
        value ~ method,
        FUN = function(x, effect) {
            mean((x - effect)^2, na.rm = TRUE)
        },
        data = df_sub_est,
        effect = effect
    )
    # Check
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


# This is a wrapper function that calls all of the other
# functions above.
get_summary_measures <- function(df, pars) {
    # Get all of the statistics
    mean_stats <- get_mean_stats(df = df)
    col_order <- names(mean_stats)
    gamma_stats <- get_gamma_stats(df = df, col_order = col_order)
    n_stats <- get_n_stats(df = df, col_order = col_order)
    bias_var_stats <- get_bias_var_stats(df = df, col_order = col_order)

    res <- rbind(
        mean_stats,
        gamma_stats,
        n_stats,
        bias_var_stats
    )
    p <- as.data.frame(lapply(pars, rep, times = nrow(res)))
    cbind(p, res)
}


## Add another wrapper that handles possible errors
calc_summary_measures <- function(df, pars, i) {
    tryCatch({
        get_summary_measures(df = df, pars = pars)
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

            # av is a list with elements that are either a tibble or NA
            # (the latter in case of an error)
            av <- vector("list", length = N)
            p_accept <- vector("numeric", length = N)
            for (i in seq_len(N)) {
                # Repeat this N times. Simulate studies, calculate CIs,
                # calculate measures
                res <- sim_effects(pars = pars, i = i)
                p_accept[i] <- attributes(res)$p_accept
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
                ci_meas <- calc_summary_measures(
                    df = df,
                    pars = pars,
                    i = i
                )
                out <- rbind(
                    ci_meas,
                    cbind(
                        pars,
                        data.frame(
                            measure = "p_accept",
                            method = "none",
                            is_ci = FALSE,
                            is_pi = FALSE,
                            is_new = FALSE,
                            value = mean(p_accept),
                            usage = "paccept_plot",
                            stat_fun = "mean",
                            stringsAsFactors = FALSE
                        )
                    )
                )
            }
        }
        # return output
        out
    }

    # # rbind ci_meas lists together
    # ci_meas <- do.call("rbind", ci_meas)

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
    dist = c("Gaussian", "t"),
    # bias
    bias = c("none", "moderate", "strong"),
    # number of large studies
    large = c(0, 1, 2),
    stringsAsFactors = FALSE
)

## run simulation, e.g., on the Rambo server of I-MATH
start <- Sys.time()
out <- sim(grid = grid, N = 1e3, cores = 120)
# out <- sim(grid = grid[688:703, ], N = 2, cores = 15)
end <- Sys.time()
print(end - start)

## save results
dir.create("RData", showWarnings = FALSE)
sessionInfo <- sessionInfo()
print(sessionInfo)
save(out, sessionInfo, file = "RData/simulate_all.RData")
