################################################################################
#                              Helper functions                                #
################################################################################

## This function repeats each element of `x` exactly `each`
## times. However, in contrast to the regular `rep()` function
## `each` can also be vector of the same length as 'x'.
## repeat each element of `x` `each` times
## rep2(c("a", "b", "c"), c(1, 2, 3)) --> c("a", "b", "b", "c", "c", "c")
rep2 <- function(x, ...) {
    unname(
        do.call(
            "c",
            mapply(FUN = rep, x = x, ... = ..., SIMPLIFY = FALSE)
        )
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

# ## These functions just call p_hmean under the hood but fix the argument
# ## `distr` to be either `"f"` or `"chisq"`. This is useful because we do not
# ## actually need to worry about including `distr = "f"` in the `args` argument
# ## when constructing expressions of the form `do.call("p_hmean", args)`.
# ## These functions are used in:
# ## - sim2CIs
# hMeanChiSqMu_f <- function(
#     estimates,
#     SEs,
#     mu = 0,
#     phi = NULL,
#     tau2 = NULL,
#     heterogeneity = "none",
#     alternative = "none",
#     check_inputs = TRUE,
#     w = rep(1, length(estimates))
# ) {
#     confMeta::p_hmean(
#         estimates = estimates,
#         SEs = SEs,
#         mu = mu,
#         phi = phi,
#         tau2 = tau2,
#         heterogeneity = heterogeneity,
#         alternative = alternative,
#         check_inputs = check_inputs,
#         w = w,
#         distr = "f"
#     )
# }
#
# hMeanChiSqMu_chisq <- function(
#     estimates,
#     SEs,
#     mu = 0,
#     phi = NULL,
#     tau2 = NULL,
#     heterogeneity = "none",
#     alternative = "none",
#     check_inputs = TRUE,
#     w = rep(1, length(estimates))
# ) {
#     confMeta::p_hmean(
#         estimates = estimates,
#         SEs = SEs,
#         mu = mu,
#         phi = phi,
#         tau2 = tau2,
#         heterogeneity = heterogeneity,
#         alternative = alternative,
#         check_inputs = check_inputs,
#         w = w,
#         distr = "chisq"
#     )
# }
