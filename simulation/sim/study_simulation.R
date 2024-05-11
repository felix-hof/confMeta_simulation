################################################################################
#                    Simulating effects and standard errors                    #
################################################################################
# 20240409: LeoUpdate: skewNormal replaces t distribution
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
#' @param dist distribution to simulate the study effect. Either "sn" or
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

    myalpha <- 4
    ## specify mean effect, variance tau2 and shape parameter alpha
    ## compute parameters xi and omega of skew-normal distribution
    ## with shape parameter alpha (default = myalpha)
    paramSN <- function(
        effect,
        tau2,
        alpha = myalpha
    ) {
        delta <- alpha / sqrt(1 + alpha^2)
        omega <- sqrt(tau2) / sqrt(1 - 2 * delta^2 / pi)
        xi <- effect - omega * delta * sqrt(2/pi)
        res <- c(xi, omega)
        names(res) <- c("xi", "omega")
        return(res)
    }


    # get args
    n <- rep(sampleSize, k)
    # include large studies
    if (large != 0) n[seq_len(large)] <- n[seq_len(large)] * 10
    # stuff for additive model
    if (heterogeneity == "additive") {
        eps2 <- 1 / k * sum(2 / n)
        tau2 <- eps2 * I2 / (1 - I2)

        if (dist == "sn") {
            ## xi: location parameter
            ## omega: scale parameter
            ## alpha: slant/shape parameter
            ## notation from Wikipedia
            p <- paramSN(effect=effect, tau2=tau2)
            delta <- sn::rsn(n = k, xi = p["xi"], omega = p["omega"], alpha = myalpha)
        } else {
            delta <- rnorm(n = k, mean = effect, sd = sqrt(tau2))
        }
        theta <- rnorm(n = k, mean = delta, sd = sqrt(2 / n))

    } else { ## multiplicative model
        phi <- 1 / (1 - I2)
        eps2 <- 1 / k * sum(2 / n)
        tau2 <- eps2 * (phi - 1)
        if (dist == "sn") {
            ## xi: location parameter
            ## omega: scale parameter
            ## alpha: slant/shape parameter
            ## notation from Wikipedia
            p <- paramSN(effect=effect, tau2=tau2)
            delta <- sn::rsn(n = k, xi = p["xi"], omega = p["omega"], alpha = myalpha)
        } else {
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
#' @param dist distribution to simulate the study effect. Either "sn" or
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
    dist = c("sn", "Gaussian"),
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
        sl <- seq_len(large)
        o <- rbind(oLarge, o[-sl, ])
        pa <- c(rep(1, large), pa[-sl])
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
