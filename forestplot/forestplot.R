library(tidyverse); theme_set(theme_bw())

#' Forest plot including harmomic mean intervals
#'
#' @param studies a data.frame with column names "y", "lower", "upper" describing the
#' study effect, the lower, and the upper bound of the CIs.
#' @param randomEffect a data.frame with column names "lower" and "upper" describing the
#' study the lower and upper bound of the CI from the random effects model.
#' @param hMean a data.frame with column names "lower" and "upper" describing the
#' study the lower and upper bound(s) of the CI(s) from the harmonic mean model.
#' @param barHeight height of the CI ticks.
#' @param arrowHeight height of the CI arrows.
#' @param lwd line width.
#' @param cex point size.
#' @return a ggplot
#' @note Adjust the vertical spacing by choosing the plotting region accordingly, e.g.,
#' \code{pdf(height=3)}.
#' @examples
#' ## tests 1 -------------------------------------
#' studies <- tibble(y = c(1, -2.3),
#'                   lower = y - .5,
#'                   upper = y + .5,
#'                   names = paste("study", seq_along(y)))
#' 
#' hMean <- tibble(lower = studies$y - .2,
#'                 upper = studies$y + .2)
#' 
#' randomEffect <- tibble(lower = mean(studies$y) - .2,
#'                        upper = mean(studies$y) + .2)
#' 
#' forest(studies, hMean=hMean, randomEffect = randomEffect)
#' 
#' 
#' 
#' ## tests 2 -------------------------------------
#' studies <- tibble(y = rnorm(10),
#'                   lower = y - .5,
#'                   upper = y + .5,
#'                   names = paste("study", seq_along(y)))
#' 
#' hMean <- tibble(lower = quantile(studies$y, probs=c(0,.5))-.1,
#'                 upper = quantile(studies$y, probs=c(.6,1))+.1)
#' 
#' randomEffect <- tibble(lower = mean(studies$y) - c(0, .2),
#'                        upper = mean(studies$y) + c(0,.2))
#' 
#' forest(studies, hMean=hMean, randomEffect = randomEffect)
#' 
#' 
#' 
#' ## tests 3 --------------------------------------
#' studies <- tibble(y = rnorm(20),
#'                   lower = y - .5,
#'                   upper = y + .5,
#'                   names = paste("study", seq_along(y)))
#' 
#' hMean <- tibble(lower = c(-2.1, -.31),
#'                 upper = lower + c(1.5, 3))
#' 
#' randomEffect <- tibble(lower = mean(studies$y) - .2,
#'                        upper = mean(studies$y) + .2)
#' 
#' forest(studies, hMean=hMean, randomEffect = randomEffect,)
#' 
forest <- function(studies,
                   randomEffect = NULL,
                   hMean = NULL,
                   barHeight = 0.5,
                   arrowHeight = 0.2,
                   lwd=1.1,
                   cex = 2){

    stopifnot(is.data.frame(studies),
              c("y", "lower", "upper") %in% names(studies),
              is.data.frame(randomEffect),
              c("lower", "upper") %in% names(randomEffect),
              is.data.frame(hMean),
              c("lower", "upper") %in% names(hMean),
              is.numeric(barHeight), length(barHeight)==1,
              is.finite(barHeight), 0 <= barHeight,
              is.numeric(arrowHeight), length(arrowHeight)==1,
              is.finite(arrowHeight), 0 <= arrowHeight,
              is.numeric(lwd), length(lwd)==1,
              is.finite(lwd), 0 <= lwd,
              is.numeric(cex), length(cex)==1,
              is.finite(cex), 0 <= cex)
    
    diamond <- function(center, height, width) {
        base <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 2) 
        trans <- rbind(base[1,]*width, base[2,]*height) + center
        geom_polygon(data = as.data.frame(t(trans)), mapping = aes(x = V1, y = V2)) 
    }

    
    mmin <- min(studies$lower, randomEffect$lower, hMean$lower)
    mmax <- max(studies$upper, randomEffect$upper, hMean$upper)

    
    p <- ggplot(data = studies,
                mapping = aes(x=y, y=rev(seq_along(y)))) + 
        geom_errorbar(width = barHeight, xmin = studies$lower, xmax = studies$upper, size=lwd) +
        geom_point(shape = 15) +
        geom_vline(xintercept = 0) +
        geom_text(aes(label = names), x = mmin - 1 ) +
        xlim(c(mmin - 1, mmax)) +
        scale_y_continuous(expand = c(.1, .1)) +
        theme_classic() +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.line.y=element_blank()) 
    
    if(!is.null(hMean)){
        p <- p + geom_segment(data = hMean,
                              mapping = aes(x = lower, xend = upper),
                              y= 0, yend=0,
                              arrow = arrow(ends="both", length = unit(arrowHeight, "in")),
                              color="blue", size=lwd) 
        p <- p + geom_text(label = "H mean", x = mmin - 1, y = 0) +
            geom_point(data = studies, mapping = aes(x=y), y=0, color="red", size=cex)
    }

    if(!is.null(randomEffect)){
        p <- p + diamond(center = c(mean(c(randomEffect$lower, randomEffect$upper)), -1),
                         height = .2, width = barHeight) 
        p <- p + geom_text(label = "RE", x = mmin - 1, y = -1 )
    }
    p
}



# ## tests 1 -------------------------------------
# studies <- tibble(y = c(1, -2.3),
#                   lower = y - .5,
#                   upper = y + .5,
#                   names = paste("study", seq_along(y)))
# 
# hMean <- tibble(lower = studies$y - .2,
#                 upper = studies$y + .2)
# 
# randomEffect <- tibble(lower = mean(studies$y) - .2,
#                        upper = mean(studies$y) + .2)
# 
# forest(studies, hMean=hMean, randomEffect = randomEffect)
# 
# 
# 
# # ## tests 2 -------------------------------------
# studies <- tibble(y = rnorm(10),
#                   lower = y - .5,
#                   upper = y + .5,
#                   names = paste("study", seq_along(y)))
# 
# hMean <- tibble(lower = quantile(studies$y, probs=c(0,.5))-.1,
#                 upper = quantile(studies$y, probs=c(.6,1))+.1)
# 
# randomEffect <- tibble(lower = mean(studies$y) - c(0, .2),
#                        upper = mean(studies$y) + c(0,.2))
# 
# forest(studies, hMean=hMean, randomEffect = randomEffect)



# ## tests 3 --------------------------------------
# studies <- tibble(y = rnorm(20),
#                   lower = y - .5,
#                   upper = y + .5,
#                   names = paste("study", seq_along(y)))
# 
# hMean <- tibble(lower = c(-2.1, -.31),
#                 upper = lower + c(1.5, 3))
# 
# randomEffect <- tibble(lower = mean(studies$y) - .2,
#                        upper = mean(studies$y) + .2)
# 
# forest(studies, hMean=hMean, randomEffect = randomEffect)
# 
