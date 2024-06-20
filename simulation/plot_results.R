# create figures from simulation output

# load libraries
rm(list = ls())
library(rlang)
library(dplyr)
library(ragg)
library(ggplot2)
library(patchwork)
theme_set(theme_bw())

# helper function to construct plot titles ----
make_title <- function(df) {
    nms <- names(df)
    nms <- vapply(
        strsplit(nms, ""),
        function(x) {
            x[1] <- toupper(x[1])
            paste0(x, collapse = "")
        },
        character(1L)
    )
    if (any(grepl("Effect", nms))) {
        nms[grepl("Effect", nms)] <- ifelse(
            grep("Effect", nms) == 1,
            "theta~\"",
            "\"~theta~\""
        )
    }
    if (any(grepl("I2", nms))) {
        nms[grepl("I2", nms)] <- ifelse(
            grep("I2", nms) == 1,
            "I^2~\"",
            "\"~I^2~\""
        )
    }
    if (any(nms == "Large")) nms[nms == "Large"] <- "No. of large studies"
    if (any(nms == "K")) nms[nms == "K"] <- "No. of studies"
    if (any(nms == "Heterogeneity")) {
        nms[nms == "Heterogeneity"] <- "Simulation model"
    }
    vals <- unname(unlist(df[1, ]))
    args <- ifelse(
        grepl("theta", nms[1]),
        paste0(
            "bquote(",
            paste0(paste0(nms, ": ", vals), collapse = ", "),
            "\")"
        ),
        paste0(
            "bquote(\"",
            paste0(paste0(nms, ": ", vals), collapse = ", "),
            "\")"
        )
    )
    return(eval(parse(text = args)))
}

# improv for now (TODO: Handle this in sim already)
convert_names <- function(id_strings) {
    st <- strsplit(id_strings, split = "_")
    get_names <- function(str) {
        meth <- c(
            "hMeanF" = "Harmonic Mean (f)",
            "hMeanChisq" = "Harmonic Mean (chisq)",
            "wilkinson" = "Wilkinson",
            "pearson" = "Pearson",
            "fisher" = "Fisher",
            "tippett" = "Tippett",
            "edgington" = "Edgington",
            "hc" =  "Henmi & Copas",
            "hk" = "Hartung & Knapp",
            "reml" = "Random effects"
        )
        type <- c(
            "ci" = "CI",
            "pi" = "PI"
        )
        # het <- c(
        #     "none" = "",
        #     "additive" = "Additive",
        #     "multiplicative" = "Multiplicative"
        # )
        tau_meth <- c(
            "none" = "none",
            "DL" = "DL",
            "PM" = "PM",
            "REML" = "REML"
        )
        paste0(
            meth[str[1L]],
            # het[str[3L]],
            " (",
            tau_meth[str[4L]],
            ") ",
            type[str[2L]]
        )
        # m <- meth[str[1L]]
        # h <- het[str[l]]
        # out <- if (h == "") paste(m, "CI") else paste(m, h, "CI")
        # if (l == 3L) out <- paste0(out, " (", str[2L], ")")
        # out
    }

    # p_accept_plots have method == "none"
    out <- rep("none", length(id_strings))
    is_none <- id_strings == "none"
    out[!is_none] <- vapply(st[!is_none], get_names, character(1L))
    out
}

to_char <- function(df) {
    out <- lapply(
        df,
        function(col) {
            as.character(col)
        }
    )
    as.data.frame(out)
}

# Load data ----
load(paste0("RData/simulate_all.RData"))

# head(out)
# tail(out)

# improv for now (TODO: once this is handled in sim, delete it here)
out <- out |>
    mutate(
        method = convert_names(method),
        I2 = factor(I2),
        k = factor(k)
    )
# anyNA(out$value)
# out |> filter(is.na(value), is_new) |> View()
# out |>
#     filter(
#         # dist == "Gaussian",
#         measure == "coverage_true",
#         I2 == 0.9,
#         effect == 0.1,
#         bias == "none"
#     ) |>
#     View()

# Prepare data (meanplots) ----
out_meanplots <- out |>
    filter(usage == "mean_plot")

# Prepare data (gammamax) ----
# out_gammamax <- out |>
#     filter(usage == "summary_stat_plot") |>
#     mutate(
#         measure = case_when(
#             stat_fun == "min" ~ "Minimum",
#             stat_fun == "firstQuart" ~ "1. Quartile",
#             stat_fun == "mean" ~ "Mean",
#             stat_fun == "median" ~ "Median",
#             stat_fun == "thirdQuart" ~ "3. Quartile",
#             stat_fun == "max" ~ "Maximum"
#         ),
#         measure = ordered(
#             measure,
#             levels = c(
#                 "Maximum", "3. Quartile",
#                 "Median", "Mean",
#                 "1. Quartile", "Minimum"
#             )
#         )
#     ) |>
#     arrange(k, dist, bias, large, heterogeneity, I2, method, measure)

# Prepare data (frequency) ----
# out_n <- out |>
#     filter(usage == "n_plot") |>
#     mutate(
#         measure = gsub("gt", "> ", stat_fun),
#         measure = ordered(measure, levels = rev(c("> 9", as.character(9:0)))),
#         value = value / attributes(out)$N
#     )

# Prepare data (summary) ----
out_sum <- out |>
    filter(measure %in% c(
        "coverage_mean", "coverage_median", "width",
        "correlation", "bias_mean", "bias_median"

    )
    )

# Prepare data (acceptance probability) ----
out_pubbias <- out |>
    filter(usage == "paccept_plot")

# Prepare data (skewness) ----
out_skew <- out |>
    filter(usage == "skewness_plot")
# Set output directory -----
# out_dir <- "~/test/Institution/harmonic_mean/figs_new"
# out_dir <- "~/ownCloud/Institution/harmonic_mean/test"
out_dir <- "newFigures"

# styles for the plots
th <- theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10),
    panel.grid = element_line(
        color = "lightgrey",
        linewidth = 0.5,
        linetype = "solid"
    )
)


#' #' Helper function to plot means
#' #' @param data obtained by simulate_all.R and saved in RData/simulate_all.RData
#' #' @param measure CI measure to plot
#' #' @param by make facets based on "method" or "I2".
plotPanels <- function(
    data,
    measure = c(
        "coverage_mean", "coverage_median", "coverage_effects",
        "coverage_effects_all",
        "coverage_effects_min1", "coverage_prediction", "n",
        "width", "score", "mse_mean", "mse_median", "var_mean",
        "var_median", "bias_mean", "bias_median", "p_accept",
        "data_skewness", "kappa", "correlation"
    ),
    by = c("method", "I2"),
    plotVar
) {

    measure <- match.arg(measure)
    by <- match.arg(by)
    plotVar <- match.arg(plotVar, choices = c("value", "prop"))
    # plotVar <- match.arg(plotVar, choices = c("value"))

    if (by == "method") {
        p <- data %>%
            # Set order of plots
            {
                if (measure == "coverage_prediction") {
                    order <- c(
                        "Harmonic Mean CI (chisq)",
                        "Harmonic Mean Additive CI (chisq)",
                        "Harmonic Mean Multiplicative CI (chisq)",
                        "Hartung & Knapp PI",
                        "Edgington CI",
                        "Edgington Additive CI",
                        "Edgington Multiplicative CI",
                        "REML PI",
                        "Fisher CI",
                        "Fisher Additive CI",
                        "Fisher Multiplicative CI"
                    )
                    .[] |>
                        filter(method %in% order) |>
                        mutate(method = factor(method, levels = order))
                } else {
                    order <- c(
                        "Edgington (none) CI",
                        "Edgington (REML) CI",
                        "Edgington (DL) CI",
                        "Fisher (none) CI",
                        "Fisher (REML) CI",
                        "Fisher (DL) CI",
                        "Pearson (none) CI",
                        "Pearson (REML) CI",
                        "Pearson (DL) CI",
                        "Tippett (none) CI",
                        "Tippett (REML) CI",
                        "Tippett (DL) CI",
                        "Wilkinson (none) CI",
                        "Wilkinson (REML) CI",
                        "Wilkinson (DL) CI",
                        "Hartung & Knapp (none) CI",
                        "Hartung & Knapp (REML) CI",
                        "Hartung & Knapp (DL) CI",
                        "Random effects (none) CI",
                        "Random effects (REML) CI",
                        "Random effects (DL) CI",
                        "Henmi & Copas (DL) CI"
                    )
                    .[] |>
                        filter(method %in% order) |>
                        mutate(method = factor(method, levels = order))
                }
            } |>
            ggplot(
                mapping = aes(x = k, y = eval(as.name(plotVar)), color = I2)
            ) +
            facet_wrap(~method) +
            geom_point() +
            stat_summary(fun = "mean", geom = "line", aes(group = I2)) +
            scale_color_discrete(name = expression(I^2)) +
            xlab("# studies")

        if (plotVar == "value") {
            p <- p + ylab(me)
        } else {
            p <- p + ylab("P(success)")
        }
    }
    if (by == "I2") {
        p <- data |>
            mutate(I2 = as.character(I2)) |>
            ggplot(mapping = aes(x = k, y = value, color = method)) +
            facet_wrap(~ I2, labeller = label_bquote(I^2 == .(I2))) +
            geom_point() +
            stat_summary(fun = "mean", geom = "line", aes(group = method)) +
            xlab("# studies")
    }

    if (stringr::str_detect(measure, "coverage")) {
        p <- p +
            geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5)
    }

    p
}

# for debugging purposes: Turn warnings into errors
options(warn = 2)



## Mean plots ------------------------------------------------------------------

out_path <- file.path(out_dir, "meanplots")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
for (i in unique(out_meanplots$measure)) {
    dir.create(file.path(out_path, i), showWarnings = FALSE)
}

opts <- list(
    dist = out_meanplots |> pull(dist) |> unique(),
    bias = out_meanplots |> pull(bias) |> unique(),
    large = out_meanplots |> pull(large) |> unique(),
    effect = out_meanplots |> pull(effect) |> unique() #,
    # heterogeneity = out_meanplots |> pull(heterogeneity) |> unique()
)
measure_opts <- out_meanplots |> pull(measure) |> unique()

list_seq <- seq_along(opts)

for (x in list_seq) { # loop over summary (eg. dist)
    current <- names(opts)[x]
    current_levels <- opts[[current]]
    cat(
        "\n",
        "Currently constructing plots for:",
        current,
        paste0("(", paste0(current_levels, collapse = ", "), ")"),
        "\n"
    )
    grid_others <- expand.grid(
        opts[list_seq[list_seq != x]],
        stringsAsFactors = FALSE
    )
    # loop over all combinations of other parameters (e.g. bias, large, effect)
    for (y in seq_len(nrow(grid_others))) {
        # filter the data by current parameter combination
        filters <- lapply(
            seq_along(grid_others),
            function(z) {
                lhs <- names(grid_others)[z]
                op <- quote(`==`)
                rhs <- grid_others[y, z]
                expr(`!!`(op)(!!sym(lhs), !!rhs))
            }
        )
        img_data <- out_meanplots |> filter(!!!filters)
        # View(img_data)
        for (me in measure_opts) { # loop over different measures
            img_data_me_ss <- img_data |> filter(measure == me)
            img_data_me_val <- img_data_me_ss |> filter(!is.na(value))
            if (!identical(img_data_me_ss, img_data_me_val)) {
                stop("Unexpected NAs.")
            }
            ll <- list(
                "value" = list(
                    data = img_data_me_val,
                    plotVar = "value"
                )#,
                # "ss" = list(
                #     data = img_data_me_ss,
                #     plotVar = "prop"
                # )
            )
            out <- lapply(
                ll,
                function(l, current_levels, current, me, th) {
                    pvar <- as.name(l$plotVar)
                    e <- bquote(c(min(.(pvar)), max(.(pvar))))
                    ylim <- eval(e, envir = l$data)
                    plots <- lapply(
                        current_levels,
                        function(z, l, current, me, th) {
                            # title
                            title <- l$data |>
                                select(
                                    -k, -I2, -method, -measure,
                                    -value, -all_of(current),
                                    -sampleSize, -usage, -stat_fun,
                                    -is_ci, -is_pi, -is_new, -prop
                                ) |>
                                distinct() |>
                                bind_cols(!!current := z) %>%
                                select(order(colnames(.))) |>
                                make_title()
                            p <- l$data |>
                                filter(!!sym(current) == z) |>
                            # data <- l$data |>
                            #     filter(!!sym(current) == z)
                            # View(data)
                            # data |> filter(grepl("Henmi", method)) |> View()
                            # data
                                plotPanels(
                                    measure = me,
                                    by = "method",
                                    plotVar = l$plotVar
                                ) +
                                ggtitle(eval(title)) +
                                th
                            if (ylim[1L] != ylim[2L]) {
                                p <- p + ylim(ylim)
                            }
                            cond_hline_0 <- me == "bias" &&
                                l$plotVar == "value" &&
                                0 > ylim[1] &&
                                0 < ylim[2]
                            if (cond_hline_0) {
                                p <- p +
                                    geom_hline(
                                        yintercept = 0,
                                        lty = 2,
                                        alpha = 0.5
                                    )
                            }
                            p
                        },
                        l = l,
                        current = current,
                        me = me,
                        th = th
                    )
                    wrap_plots(plots)
                },
                current_levels = current_levels,
                current = current,
                me = me,
                th = th
            )
            res <- wrap_plots(out, guides = "collect", ncol = 1) &
                theme(legend.position = "bottom")
            # generate the file names
            fname <- paste0(
                me, "/", toupper(current), "_",
                paste0(
                    paste0(
                        names(grid_others[y, ]),
                        "_", to_char(grid_others[y, ])
                    ),
                    collapse = "_"
                ),
                ".png"
            )
            filename <- file.path(out_path, fname)
            cat("\33[2K\rMaking file: ", filename)
            ggsave(
                filename = filename,
                plot = res,
                width = length(current_levels) * 6.5,
                height = 24,
                units = "in",
                device = ragg::agg_png
            )
        }
    }
}

## gamma_max summary statistics ------------------------------------------------

# out_path <- file.path(out_dir, "max_pH")
# dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
#
# opts <- list(
#     dist = unique(out_gammamax$dist),
#     bias = unique(out_gammamax$bias),
#     large = unique(out_gammamax$large),
#     # heterogeneity =  unique(out_gammamax$heterogeneity),
#     effect = unique(out_gammamax$effect),
#     I2 =  unique(out_gammamax$I2)
# )
# measure_opts <- unique(out_gammamax$measure)
#
# list_seq <- seq_along(opts)
#
# # Set order of the plots
# ord <- c(
#     "Edgington (none) CI",
#     "Edgington (REML) CI",
#     "Edgington (DL) CI",
#     "Edgington (PM) CI",
#     "Fisher (none) CI",
#     "Fisher (REML) CI",
#     "Fisher (DL) CI",
#     "Fisher (PM) CI",
#     "Pearson (none) CI",
#     "Pearson (REML) CI",
#     "Pearson (DL) CI",
#     "Pearson (PM) CI",
#     "Random effects (none) CI",
#     "Random effects (REML) CI",
#     "Random effects (DL) CI",
#     "Random effects (PM) CI"
# )
#
# ## For now, subset to only those methods we actually want plots of
# # out_gammamax <- subset(out_gammamax, method %in% ord)
#
# for (x in list_seq) { # loop over summary (eg. dist)
#     current <- names(opts)[x]
#     current_levels <- opts[[current]]
#     cat(
#         "Currently constructing plots for:",
#         current,
#         paste0("(", paste0(current_levels, collapse = ", "), ")"),
#         "\n"
#     )
#     grid_others <- expand.grid(
#         opts[list_seq[list_seq != x]],
#         stringsAsFactors = FALSE
#     )
#     # loop over all combinations of other parameters (e.g. bias, large, he)
#     for (y in seq_len(nrow(grid_others))) {
#         # filter the data by current parameter combination
#         filters <- lapply(seq_along(grid_others), function(z) {
#             lhs <- names(grid_others)[z]
#             op <- quote(`==`)
#             rhs <- grid_others[y, z]
#             expr(`!!`(op)(!!sym(lhs), !!rhs))
#         })
#         img_data_ss <- out_gammamax |>
#             filter(!!!filters) |>
#             mutate(method = factor(method, levels = ord))
#         img_data_val <- img_data_ss |> filter(!is.na(value))
#         ll <- list(
#             "value" = list(
#                 data = img_data_val,
#                 plotVar = "value"
#             ) #,
#             # "ss" = list(
#             #     data = img_data_ss,
#             #     plotVar = "prop"
#             # )
#         )
#         out <- lapply(
#             ll,
#             function(l, current_levels, current, th) {
#                 pvar <- as.name(l$plotVar)
#                 e <- bquote(c(min(.(pvar)), max(.(pvar))))
#                 ylim <- eval(e, envir = l$data)
#                 plots <- lapply(
#                     current_levels,
#                     function(z, l, current, th) {
#                         # title
#                         title <- l$data |>
#                             select(
#                                 -k, -I2, -method, -measure,
#                                 -value, -all_of(current),
#                                 -sampleSize, -usage, -stat_fun,
#                                 -is_ci, -is_pi, -is_new, -prop
#                             ) |>
#                             distinct() |>
#                             bind_cols(!!current := z) %>%
#                             select(order(colnames(.))) |>
#                             make_title()
#                         l$data |>
#                             filter(!!sym(current) == z) |>
#                             ggplot(
#                                 aes(
#                                     x = k,
#                                     y = eval(as.name(l$plotVar)),
#                                     color = measure
#                                 )
#                             ) +
#                             facet_wrap(~method, ncol = 4) +
#                             geom_point() +
#                             stat_summary(
#                                 fun = "mean",
#                                 geom = "line",
#                                 aes(group = measure)
#                             ) +
#                             ylim(c(0, 1)) +
#                             labs(
#                                 x = "# studies",
#                                 y = if (l$plotVar == "value") {
#                                     "highest p-value"
#                                 } else {
#                                     "P(success)"
#                                 },
#                                 color = "Summary statistic",
#                                 title = eval(title)
#                             ) +
#                             th
#                     },
#                     l = l,
#                     current = current,
#                     th = th
#                 )
#                 wrap_plots(plots)
#             },
#             current_levels = current_levels,
#             current = current,
#             th = th
#         )
#         res <- wrap_plots(out, guides = "collect", ncol = 1) &
#             theme(legend.position = "bottom")
#         # generate the file names
#         fname <- paste0(
#             out_path, "/", toupper(current), "_",
#             paste0(
#                 paste0(
#                     names(grid_others[y, ]),
#                     "_", to_char(grid_others[y, ])
#                 ),
#                 collapse = "_"
#             ),
#             ".png"
#         )
#         # filename <- file.path(out_path, fname)
#         filename <- fname
#         cat("\33[2K\rMaking file: ", filename)
#         ggsave(
#             filename = filename,
#             plot = res,
#             width = length(current_levels) * 6.5,
#             height = 24,
#             units = "in",
#             device = ragg::agg_png
#         )
#     }
# }

## relative frequencies --------------------------------------------------------

# out_path <- file.path(out_dir, "rel_freq")
# dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
#
# # Adapt make_title function for now
# make_title <- function(df) {
#     nms <- names(df)
#     nms <- vapply(
#         strsplit(nms, ""),
#         function(x) {
#             x[1] <- toupper(x[1])
#             paste0(x, collapse = "")
#         },
#         character(1L)
#     )
#     if (any(grepl("Effect", nms))) {
#         nms[grepl("Effect", nms)] <- ifelse(
#             grep("Effect", nms) == 1,
#             "theta~\"",
#             "\"~theta~\""
#         )
#     }
#     if (any(grepl("I2", nms))) {
#         nms[grepl("I2", nms)] <- ifelse(
#             grep("I2", nms) == 1,
#             "I^2~\"",
#             "\"~I^2~\""
#         )
#     }
#     if (any(nms == "Large")) nms[nms == "Large"] <- "# large studies"
#     if (any(nms == "K")) nms[nms == "K"] <- "# studies"
#     if (any(nms == "Heterogeneity")) {
#         nms[nms == "Heterogeneity"] <- "Simulation model"
#     }
#     vals <- unname(unlist(df[1, ]))
#     args <- ifelse(
#         grepl("theta", nms[1]),
#         paste0(
#             "bquote(", paste0(paste0(nms, ": ", vals), collapse = ", "),
#             "\")"
#         ),
#         paste0(
#             "bquote(\"", paste0(paste0(nms, ": ", vals), collapse = ", "),
#             "\")"
#         )
#     )
#     return(eval(parse(text = args)))
# }
#
# # These are the options (i.e. what variables to compare and their levels)
# opts <- list(
#     dist = unique(out_n$dist),
#     bias = unique(out_n$bias),
#     large = unique(out_n$large),
#     effect = unique(out_n$effect),
#     # heterogeneity = unique(out_n$heterogeneity),
#     I2 = unique(out_n$I2),
#     k = unique(out_n$k)
# )
#
# # What measures do we have here (this ix the x axis)
# measure_opts <- unique(out_n$measure)
#
# list_seq <- seq_along(opts)
#
# # Set order of the plots
# ord <- c(
#     "Edgington (none) CI",
#     "Edgington (REML) CI",
#     "Edgington (DL) CI",
#     "Edgington (PM) CI",
#     "Fisher (none) CI",
#     "Fisher (REML) CI",
#     "Fisher (DL) CI",
#     "Fisher (PM) CI",
#     "Pearson (none) CI",
#     "Pearson (REML) CI",
#     "Pearson (DL) CI",
#     "Pearson (PM) CI",
#     "Random effects (none) CI",
#     "Random effects (REML) CI",
#     "Random effects (DL) CI",
#     "Random effects (PM) CI"
# )
#
# # ord <- c(
# #     "Edgington (none) CI",
# #     "Edgington (REML) CI",
# #     "Edgington (DL) CI",
# #     "Edgington (PM) CI",
# #     "Fisher (none) CI",
# #     "Fisher (REML) CI",
# #     "Fisher (DL) CI",
# #     "Fisher (PM) CI"
# # )
#
# ## For now, subset to only those methods we actually want plots of
# # out_n <- subset(out_n, method %in% ord)
#
#
# # loop over summary (eg. dist)
# for (x in list_seq) {
#     current <- names(opts)[x]
#     current_levels <- opts[[current]]
#     cat(
#         "Currently constructing plots for:",
#         current,
#         paste0("(", paste0(current_levels, collapse = ", "), ")"),
#         "\n"
#     )
#     grid_others <- expand.grid(
#         opts[list_seq[list_seq != x]],
#         stringsAsFactors = FALSE
#     )
#     # loop over all combinations of other parameters (e.g. bias, large, he)
#     for (y in seq_len(nrow(grid_others))) {
#         # filter the data by current parameter combination
#         filters <- lapply(seq_along(grid_others), function(z) {
#             lhs <- names(grid_others)[z]
#             op <- quote(`==`)
#             rhs <- grid_others[y, z]
#             expr(`!!`(op)(!!sym(lhs), !!rhs))
#         })
#         # out_n |> select(method) |> unique()
#         # ord
#         img_data <- out_n |>
#             filter(!!!filters) |>
#             mutate(
#                 method = factor(method, levels = ord)
#             )
#         lapply(current_levels, function(z) {
#             title <- img_data |>
#                 filter(!!sym(current) == z) |>
#                 select(
#                     -method, -measure, -value,
#                     -all_of(current), -sampleSize,
#                     -usage, -stat_fun, -is_ci, -is_pi,
#                     -is_new
#                 ) |>
#                 distinct() |>
#                 bind_cols(!!current := z) %>%
#                 select(order(colnames(.))) |>
#                 make_title()
#             img_data |>
#                 filter(!!sym(current) == z) |>
#                 ggplot(aes(x = measure, y = value)) +
#                 geom_col(fill = viridisLite::viridis(1)) +
#                 ylim(c(0, 1)) +
#                 facet_wrap(~method, ncol = 4) +
#                 labs(
#                     x = "# intervals",
#                     y = "relative frequency",
#                     title = eval(title)
#                 )
#         }) |>
#             wrap_plots(guides = "collect") &
#             theme(
#                 legend.position = "bottom",
#                 text = element_text(size = 9),
#                 plot.title = element_text(size = 10)
#             )
#         filename <- paste0(
#             out_path, "/", toupper(current), "_",
#             paste0(
#                 paste0(
#                     names(grid_others[y, ]), "_",
#                     to_char(grid_others[y, ])
#                 ),
#                 collapse = "_"
#             ),
#             ".png"
#         )
#         cat("\33[2K\rMaking file: ", filename)
#         ggsave(
#             filename = filename,
#             width = length(current_levels) * 7.5,
#             height = 12,
#             units = "in",
#             device = ragg::agg_png
#         )
#     }
# }

# p_accept plots ---------------------------------------------------------------

# delete the following once the simulation finished
out_pubbias <- out_pubbias |>
    filter(value != 1)

out_path <- file.path(out_dir, "p_accept")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
for (i in unique(out_pubbias$measure))
    dir.create(file.path(out_path, i), showWarnings = FALSE)

opts <- list(
    dist = out_pubbias |> pull(dist) |> unique(),
    bias = out_pubbias |> pull(bias) |> unique(),
    large = out_pubbias |> pull(large) |> unique(),
    effect = out_pubbias |> pull(effect) |> unique() #,
    # heterogeneity = out_pubbias |> pull(heterogeneity) |> unique()
)
measure_opts <- out_pubbias |> pull(measure) |> unique()

list_seq <- seq_along(opts)

for (x in list_seq) { # loop over summary (eg. dist)
    current <- names(opts)[x]
    current_levels <- opts[[current]]
    cat(
        "\n",
        "Currently constructing plots for:",
        current,
        paste0("(", paste0(current_levels, collapse = ", "), ")"),
        "\n"
    )
    grid_others <- expand.grid(
        opts[list_seq[list_seq != x]],
        stringsAsFactors = FALSE
    )
    # loop over all combinations of other parameters (e.g. bias, large, effect)
    for (y in seq_len(nrow(grid_others))) {
        # filter the data by current parameter combination
        filters <- lapply(seq_along(grid_others), function(z) {
            lhs <- names(grid_others)[z]
            op <- quote(`==`)
            rhs <- grid_others[y, z]
            expr(`!!`(op)(!!sym(lhs), !!rhs))
        })
        img_data <- out_pubbias %>% filter(!!!filters)
        me <- "p_accept"
        img_data_me <- img_data |> filter(measure == me)
        ylim <- with(img_data_me, c(min(value), max(value)))
        lapply(
            current_levels,
            function(z) {
                title <- img_data_me |>
                    select(
                        -k, -I2, -method, -measure,
                        -value, -all_of(current),
                        -sampleSize, -usage, -stat_fun,
                        -is_ci, -is_pi, -is_new
                    ) |>
                    distinct() |>
                    bind_cols(!!current := z) %>%
                    select(order(colnames(.))) |>
                    make_title()
                img_data_me |>
                    filter(!!sym(current) == z) |>
                    ggplot(aes(x = k, y = value, color = I2)) +
                    geom_point() +
                    stat_summary(
                        fun = "mean",
                        geom = "line",
                        aes(group = I2)
                    ) +
                    scale_color_discrete(name = expression(I^2)) +
                    ylim(ylim) +
                    xlab("# studies") +
                    ylab(me) +
                    ggtitle(eval(title)) +
                    theme(
                        text = element_text(size = 9),
                        plot.title = element_text(size = 10)
                    )
            }
        ) |>
            wrap_plots(guides = "collect") &
            theme(legend.position = "bottom")
        filename <- paste0(
            out_path, "/", toupper(current), "_",
            paste0(
                paste0(
                    names(grid_others[y, ]),
                    "_", to_char(grid_others[y, ])
                ),
                collapse = "_"
            ),
            ".png"
        )
        cat("\33[2K\rMaking file: ", filename)
        ggsave(
            filename = filename,
            width = length(current_levels) * 6.5,
            height = 12,
            units = "in",
            device = ragg::agg_png
        )
    }
}


# skewness plots ---------------------------------------------------------------

out_path <- file.path(out_dir, "skewness")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
for (i in unique(out_skew$measure)) {
    dir.create(file.path(out_path, i), showWarnings = FALSE)
}

opts <- list(
    dist = out_skew |> pull(dist) |> unique(),
    bias = out_skew |> pull(bias) |> unique(),
    large = out_skew |> pull(large) |> unique(),
    effect = out_skew |> pull(effect) |> unique() #,
    # heterogeneity = out_skew |> pull(heterogeneity) |> unique()
)
measure_opts <- out_skew |> pull(measure) |> unique()

list_seq <- seq_along(opts)

for (x in list_seq) { # loop over summary (eg. dist)
    current <- names(opts)[x]
    current_levels <- opts[[current]]
    cat(
        "\n",
        "Currently constructing plots for:",
        current,
        paste0("(", paste0(current_levels, collapse = ", "), ")"),
        "\n"
    )
    grid_others <- expand.grid(
        opts[list_seq[list_seq != x]],
        stringsAsFactors = FALSE
    )
    # loop over all combinations of other parameters (e.g. bias, large, effect)
    for (y in seq_len(nrow(grid_others))) {
        # filter the data by current parameter combination
        filters <- lapply(
            seq_along(grid_others),
            function(z) {
                lhs <- names(grid_others)[z]
                op <- quote(`==`)
                rhs <- grid_others[y, z]
                expr(`!!`(op)(!!sym(lhs), !!rhs))
            }
        )
        img_data <- out_skew |> filter(!!!filters)
        # View(img_data)
        for (me in measure_opts) { # loop over different measures
            img_data_me_ss <- img_data |> filter(measure == me)
            img_data_me_val <- img_data_me_ss |> filter(!is.na(value))
            # if (!identical(img_data_me_ss, img_data_me_val)) {
            #     stop("Unexpected NAs.")
            # }
            ll <- list(
                "value" = list(
                    data = img_data_me_val,
                    plotVar = "value"
                )#,
                # "ss" = list(
                #     data = img_data_me_ss,
                #     plotVar = "prop"
                # )
            )
            out <- lapply(
                ll,
                function(l, current_levels, current, me, th) {
                    pvar <- as.name(l$plotVar)
                    e <- bquote(c(min(.(pvar)), max(.(pvar))))
                    ylim <- eval(e, envir = l$data)
                    plots <- lapply(
                        current_levels,
                        function(z, l, current, me, th) {
                            # title
                            title <- l$data |>
                                select(
                                    -k, -I2, -method, -measure,
                                    -value, -all_of(current),
                                    -sampleSize, -usage, -stat_fun,
                                    -is_ci, -is_pi, -is_new, -prop
                                ) |>
                                distinct() |>
                                bind_cols(!!current := z) %>%
                                select(order(colnames(.))) |>
                                make_title()
                            p <- l$data |>
                                filter(!!sym(current) == z) |>
                            # data <- l$data |>
                            #     filter(!!sym(current) == z)
                            # View(data)
                            # data |> filter(grepl("Henmi", method)) |> View()
                            # data
                                plotPanels(
                                    measure = me,
                                    by = "method",
                                    plotVar = l$plotVar
                                ) +
                                ggtitle(eval(title)) +
                                th
                            if (ylim[1L] != ylim[2L]) {
                                p <- p + ylim(ylim)
                            }
                            cond_hline_0 <- me == "bias" &&
                                l$plotVar == "value" &&
                                0 > ylim[1] &&
                                0 < ylim[2]
                            if (cond_hline_0) {
                                p <- p +
                                    geom_hline(
                                        yintercept = 0,
                                        lty = 2,
                                        alpha = 0.5
                                    )
                            }
                            p
                        },
                        l = l,
                        current = current,
                        me = me,
                        th = th
                    )
                    wrap_plots(plots)
                },
                current_levels = current_levels,
                current = current,
                me = me,
                th = th
            )
            res <- wrap_plots(out, guides = "collect", ncol = 1) &
                theme(legend.position = "bottom")
            # generate the file names
            fname <- paste0(
                me, "/", toupper(current), "_",
                paste0(
                    paste0(
                        names(grid_others[y, ]),
                        "_", to_char(grid_others[y, ])
                    ),
                    collapse = "_"
                ),
                ".png"
            )
            filename <- file.path(out_path, fname)
            cat("\33[2K\rMaking file: ", filename)
            ggsave(
                filename = filename,
                plot = res,
                width = length(current_levels) * 6.5,
                height = 24,
                units = "in",
                device = ragg::agg_png
            )
        }
    }
}


# Summary of meanplots ---------------------------------------------------------

# Debug
# out_sum_og <- out_sum

# Get the paramSN function to calculate true median
e <- new.env()
source(file = "study_simulation.R", local = e)
paramSN <- e$paramSN

# effect = 0.2
# k = 3
# I2 = 0.6
# dist = "snr"
# sampleSize = 50
# large = 2
# heterogeneity = "additive"
# alpha = 8

calc_median_sn <- function(effect, I2, k, dist, sampleSize, large, heterogeneity, alpha){

    if (is.factor(k)) k <- as.numeric(as.character(k))
    if (is.factor(I2)) I2 <- as.numeric(as.character(I2))
    args <- list(effect, I2, k, dist, sampleSize, large, heterogeneity)
    print(sapply(args, class))
    lengths <- sapply(args, length)
    if (!all(lengths == lengths[1])) stop("Length mismatch.")
    l <- lengths[1L]

    calc_one_median <- function(effect, dist, I2, k, sampleSize, large, heterogeneity, alpha) {
        if (dist == "Gaussian") return(effect)
        # alpha <- 8
        n <- rep(sampleSize, k)
        if (large != 0) n[seq_len(large)] <- n[seq_len(large)] * 10
        if (heterogeneity == "additive") {
            eps2 <- 1/ k * sum(2 / n)
            tau2 <- eps2 * I2 / (1 - I2)
            if (dist == "snl") alpha <- -alpha
            med <- paramSN(
                effect = effect,
                tau2 = tau2,
                alpha = alpha,
                median = TRUE
            )["median"]
        } else {
            stop("heterogeneity == 'multiplicative' is not implemented.")
        }
        med
    }

    res <- numeric(l)
    for (i in seq_len(l)) {
        res[i] <- calc_one_median(
            effect = effect[i],
            dist = dist[i],
            I2 = I2[i],
            k = k[i],
            sampleSize = sampleSize[i],
            large = large[i],
            heterogeneity = heterogeneity[i],
            alpha = 8
        )
    }
    res
}

# Try to create relative bias
within(
    out_sum |> filter(measure %in% c("bias_mean", "bias_median")),
    {
        # means are simple
        means <- effect
        # medians, not so much
        medians <- calc_median_sn(
                effect = effect,
                k = k,
                I2 = I2,
                dist = dist,
                sampleSize = sampleSize,
                large = large,
                heterogeneity = heterogeneity,
                alpha = 8
        )
    }
)

p <- paramSN(
    effect = out_sum$effect,
    tau2 = tau2,
    alpha = alpha,
    median = TRUE
)
e$simRE

# Calculate relative bias


# Set the order of plots
method_order <- c(
    # "Edgington (none) CI",
    # "Edgington (REML) CI",
    "Edgington (DL) CI",
    # "Fisher (none) CI",
    # "Fisher (REML) CI",
    # "Fisher (DL) CI",
    # "Pearson (none) CI",
    # "Pearson (REML) CI",
    # "Pearson (DL) CI",
    # "Tippett (none) CI",
    # "Tippett (REML) CI",
    # "Tippett (DL) CI",
    # "Wilkinson (none) CI",
    # "Wilkinson (REML) CI",
    # "Wilkinson (DL) CI",
    # "Random effects (none) CI",
    # "Random effects (REML) CI",
    "Random effects (DL) CI",
    # "Hartung & Knapp (none) CI",
    # "Hartung & Knapp (REML) CI",
    "Hartung & Knapp (DL) CI"#,
    # "Henmi & Copas (DL) CI"
)

# At this point out_sum only has "coverage_true".
# Here we use all methods, to kick some out just commment them above
out_sum <- subset(out_sum, method %in% method_order)

# Rename bias to PBias
out_sum <- within(out_sum, {PBias <- bias; rm(bias)})


# Remove certain levels from the plots
out_sum <- out_sum |> filter(dist != "snr")

# Convert publication bias and method to ordered factor
out_sum <- within(
    out_sum,
    {
        PBias <- factor(PBias, levels = c("none", "moderate", "strong"), ordered = TRUE)
        method <- factor(method, levels = method_order, ordered = TRUE)
    }
)



# Make directories
out_path <- file.path(out_dir, "summary")
overall_out_path <- file.path(out_path, "overall")
strat_out_path <- file.path(out_path, "stratified")
dir.create(overall_out_path, showWarnings = FALSE, recursive = TRUE)
dir.create(strat_out_path, showWarnings = FALSE, recursive = TRUE)

# Set types of plots
opts <- list(
    dist = unique(out_sum$dist),
    effect = unique(out_sum$effect),
    large = unique(out_sum$large),
    PBias = unique(out_sum$PBias)
)

msrs <- unique(out_sum$measure)


# Convert this to factor



# totitle <- function(string) {
#     stopifnot(
#         is.character(string),
#         length(string) == 1L
#     )
#     string <- gsub("_", " ", string)
#     string <- unlist(strsplit(string, split = " "))
#     string <- strsplit(string, split = "")
#     string <- vapply(string, function(y) {
#         y[1] <- toupper(y[1])
#         paste0(y, collapse = "")
#     }, character(1L))
#     paste0(string, collapse = " ")
# }

get_title_param <- function(param) {
    switch(
        param,
        "dist" = "distribution",
        "eff" = "effect",
        "large" = "number of large studies",
        "PBias" = "publication bias"
    )
}

get_title_me <- function(me) {
    switch(
        me,
        "coverage_mean" = "coverage of the mean",
        "coverage_median" = "coverage of the median",
        "bias_mean" = "bias of the mean",
        "bias_median" = "bias of the median",
        "width" = "width",
        "correlation" = "correlation",
    )
}

transparency <- 0.6

for (x in seq_along(msrs)) {
    # get the data we need on the y-axis
    me <- msrs[x]
    cur_data <- out_sum |> filter(measure == me)
    for (y in seq_along(opts)) {
        # Get facets
        param <- names(opts)[y]
        levels <- opts[[y]]
        # construct a list of names that we want to group the data by
        gb_o <- lapply(list("method", "k", "I2"), as.name)
        gb_o <- append(gb_o, list(as.name(param)))
        # Set some variables for the plot title, as these are used in all plots
        # title_me <- get_title_me(me = me)
        # title_param <- get_title_param(param = param)
        # plot_title <- paste0(
        #     "Minimum, mean, and maximum ",
        #     title_me,
        #     " by method and ",
        #     title_param
        # )
        # Do stuff for stratification
        ## if the facet variable is dist or PBias, also add the other one as a
        ## grouping variable for stratified plots
        is_dist <- param == "dist"
        is_PBias <- param == "PBias"
        do_strat <- is_dist || is_PBias
        if (do_strat) {
            # Generate the group by variables
            gb_s <- if (is_dist) {
                append(gb_o, list(as.name("PBias")))
            } else if (is_PBias) {
                append(gb_o, list(as.name("dist")))
            }
            # Determine the other parameter and levels
            other_param <- if (is_dist) "PBias" else "dist"
            other_levels <- opts[[other_param]]
            # Get the stratified data ready
            img_data_s <- cur_data |>
                group_by(!!!gb_s) |>
                summarise(
                    min = min(value),
                    mean = mean(value),
                    max = max(value),
                    .groups = "drop"
                )
            # Get the first parts of the plot title
            # plot_title_s <- paste0(
            #     plot_title,
            #     " with ",
            #     get_title_param(param = other_param),
            #     " = "
            # )
            # Loop over the other_param levels and create a plot for each
            for (l in other_levels) {
                filt <- bquote(expr = .(as.name(other_param)) == .(as.character(l)))
                dat_sub <- img_data_s |> filter(eval(filt))
                # title_s <- paste0(plot_title_s, l)
                ylimes_s <- with(dat_sub, c(min(min), max(max)))
                # Loop over the methods
                plots <- lapply(seq_along(method_order), function(m) {
                    meth <- method_order[m]
                    # filt <- bquote(expr = .(as.name(other_param)) == .(as.character(l)))
                    dat <- dat_sub |> filter(method == meth)
                    # if no data, don't produce a plot
                    if (nrow(dat) == 0L) return(NULL)
                    ref0.95 <- grepl("coverage", me, fixed = TRUE)
                    ref0 <- grepl("correlation|bias", me, fixed = FALSE)
                    p <- ggplot(dat, aes(x = k, y = mean, color = I2, group = I2))
                    if (ref0.95) {
                        yint <- 0.95
                    } else if (ref0) {
                        yint <- 0
                    }
                    if ((ref0.95 || ref0) && (yint > ylimes_s[1] && yint < ylimes_s[2])) {
                        p <- p +
                        geom_hline(yintercept = yint, color = "black", lty = 2, alpha = 0.5)
                    }
                    p <- p +
                    geom_point(position = position_dodge(width = 0.5)) +
                    geom_errorbar(
                        aes(ymin = min, ymax = max),
                        # position = "dodge",
                        position = position_dodge(width = 0.5),
                        width = 0.4
                        # alpha = transparency
                    ) +
                    # scale_color_discrete(name = expression(I^2)) +
                    facet_wrap(~eval(as.name(param)), ncol = 1L) +
                    ylim(ylimes_s) +
                    th +
                    theme(
                        legend.position = "bottom",
                        plot.title = element_text(hjust = 0.5)
                    ) #+
                    # labs(x = "k", y = me, color = bquote(I^2), title = meth)
                    if (m == 1L) {
                        p <- p + labs(x = "k", y = me, color = bquote(I^2), title = meth)
                    } else {
                        p <- p + labs(x = "k", y = NULL, color = bquote(I^2), title = meth)
                    }
                    p
                }) 
                plots <- plots[!sapply(plots, is.null)]
                plots <- plots |>
                    wrap_plots(guides = "collect", ncol = length(plots)) &
                    # plot_annotation(title = title_s) &
                    theme(legend.position = "bottom")
                filename <- file.path(
                    strat_out_path,
                    paste0(
                        toupper(param), "_", me, "_", other_param,
                        "_", as.character(l), ".png"
                    )
                )
                cat("\33[2K\rMaking file: ", filename, fill = TRUE)
                ggsave(
                    filename = filename,
                    device = ragg::agg_png,
                    plot = plots,
                    width = length(levels) * 5,
                    height = length(method_order) * 5,
                    units = "in"
                )
            }
        }
        img_data_o <- cur_data |>
            group_by(!!!gb_o) |>
            summarise(
                min = min(value),
                mean = mean(value),
                max = max(value),
                .groups = "drop"
            )
        # Some things for the overall plots
        ylimes <- c(min(img_data_o$min), max(img_data_o$max))
        plots <- lapply(seq_along(method_order), function(m) {
            meth <- method_order[m]
            dat <- img_data_o |>
                filter(method == meth)
            # if no data, don't produce a plot
            if (nrow(dat) == 0L) return(NULL)
            p <- ggplot(dat, aes(x = k, y = mean, color = I2, group = I2))
            ref0.95 <- grepl("coverage", me, fixed = TRUE)
            ref0 <- grepl("correlation|bias", me, fixed = FALSE)
            if (ref0.95) {
                yint <- 0.95
            } else if (ref0) {
                yint <- 0
            }
            if ((ref0.95 || ref0) && (yint > ylimes_s[1] && yint < ylimes_s[2])) {
                p <- p +
                geom_hline(yintercept = yint, color = "black", lty = 2, alpha = 0.5)
            }
            p <- p +
            geom_point(position = position_dodge(width = 0.5)) +
            geom_errorbar(
                aes(ymin = min, ymax = max),
                # position = "dodge",
                position = position_dodge(width = 0.5),
                width = 0.4
                # alpha = transparency
            ) +
            # scale_color_discrete(name = expression(I^2)) +
            facet_wrap(~eval(as.name(param)), ncol = 1L) +
            ylim(ylimes_s) +
            th +
            theme(
                legend.position = "bottom",
                plot.title = element_text(hjust = 0.5)
            ) #+
            # labs(x = "k", y = me, color = bquote(I^2), title = meth)
            if (m == 1L) {
                p <- p + labs(x = "k", y = me, color = bquote(I^2), title = meth)
            } else {
                p <- p + labs(x = "k", y = NULL, color = bquote(I^2), title = meth)
            }
            p
            # p <- p +
            # geom_point(position = position_dodge(width = 0.5)) +
            #     geom_errorbar(
            #         aes(ymin = min, ymax = max),
            #         # position = "dodge",
            #         position = position_dodge(width = 0.5),
            #         width = 0.4
            #         # alpha = transparency
            #     ) +
            #     # scale_color_discrete(name = expression(I^2)) +
            #     facet_wrap(~eval(as.name(param))) +
            #     ylim(ylimes) +
            #     th +
            #     theme(legend.position = "bottom") +
            #     labs(y = z, color = bquote(I^2))
        }) 
        plots <- plots[!sapply(plots, is.null)]
        plots <- plots |>
            wrap_plots(guides = "collect", ncol = 1L) &
            # plot_annotation(title = plot_title) &
            theme(legend.position = "bottom")
        filename <- file.path(
            overall_out_path,
            paste0(toupper(me), "_", param, ".png")
        )
        cat("\33[2K\rMaking file: ", filename, fill = TRUE)
        ggsave(
            filename = filename,
            device = ragg::agg_png,
            plot = plots,
            width = length(levels) * 5,
            height = length(method_order) * 5,
            units = "in"
        )
    }
}
