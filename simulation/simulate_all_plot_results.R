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
            "ktrials" = "k-Trials",
            "pearson" = "Pearson",
            "fisher" = "Fisher",
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
    vapply(st, get_names, character(1L))
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

# Prepare data (meanplots) ----
out_meanplots <- out %>%
    filter(usage == "mean_plot")

# Prepare data (gammamax) ----
out_gammamax <- out %>%
    filter(usage == "summary_stat_plot") %>%
    mutate(
        measure = case_when(
            stat_fun == "min" ~ "Minimum",
            stat_fun == "firstQuart" ~ "1. Quartile",
            stat_fun == "mean" ~ "Mean",
            stat_fun == "median" ~ "Median",
            stat_fun == "thirdQuart" ~ "3. Quartile",
            stat_fun == "max" ~ "Maximum"
        ),
        measure = ordered(
            measure,
            levels = c(
                "Maximum", "3. Quartile",
                "Median", "Mean",
                "1. Quartile", "Minimum"
            )
        )
    ) %>%
    arrange(k, dist, bias, large, heterogeneity, I2, method, measure)

# Prepare data (frequency) ----
out_n <- out %>%
    filter(usage == "n_plot") %>%
    mutate(
        measure = gsub("gt", "> ", stat_fun),
        measure = ordered(measure, levels = rev(c("> 9", as.character(9:0)))),
        value = value / attributes(out)$N
    )

# Prepare data (summary) ----
out_sum <- out_meanplots %>%
    filter(grepl("coverage", measure))

# Prepare data (summary) ----
out_pubbias <- out %>%
    filter(usage == "paccept_plot")

# Set output directory -----
# out_dir <- "~/test/Institution/harmonic_mean"
out_dir <- "~/ownCloud/Institution/harmonic_mean"

#' #' Helper function to plot means
#' #' @param data obtained by simulate_all.R and saved in RData/simulate_all.RData
#' #' @param measure CI measure to plot
#' #' @param by make facets based on "method" or "I2".
plotPanels <- function(
    data,
    measure = c(
        "coverage_true", "coverage_effects", "coverage_effects_all",
        "coverage_effects_min1", "coverage_prediction", "n",
        "width", "score", "mse", "var", "bias", "p_accept"
    ),
    by = c("method", "I2")
) {

    measure <- match.arg(measure)
    by <- match.arg(by)
    # data <- data %>%
    #     mutate(
    #         k = factor(k),
    #         I2 = factor(I2)
    #     )
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
                    .[] %>%
                        filter(method %in% order) %>%
                        mutate(method = factor(method, levels = order))
                } else {
                    order <- c(
                        "Edgington (none) CI",
                        "Edgington (REML) CI",
                        "Edgington (DL) CI",
                        "Edgington (PM) CI",
                        "Fisher (none) CI",
                        "Fisher (REML) CI",
                        "Fisher (DL) CI",
                        "Fisher (PM) CI",
                        "Random effects (none) CI",
                        "Random effects (REML) CI",
                        "Random effects (DL) CI",
                        "Random effects (PM) CI",
                        "Hartung & Knapp (REML) CI",
                        "Henmi & Copas (DL) CI"
                    )
                    .[] %>%
                        filter(method %in% order) %>%
                        mutate(method = factor(method, levels = order))
                }
            } %>%
            ggplot(mapping = aes(x = k, y = value, color = I2)) +
            facet_wrap(~method) +
            geom_point() +
            stat_summary(fun = "mean", geom = "line", aes(group = I2)) +
            scale_color_discrete(name = expression(I^2)) +
            xlab("# studies") +
            ylab(measure)
    }
    if (by == "I2") {
        p <- data %>%
            mutate(I2 = as.character(I2)) %>%
            ggplot(mapping = aes(x = k, y = value, color = method)) +
            facet_wrap(~ I2, labeller = label_bquote(I^2 == .(I2))) +
            geom_point() +
            stat_summary(fun = "mean", geom = "line", aes(group = method)) +
            xlab("# studies") +
            ylab(measure)
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

## 4.1
out_path <- file.path(out_dir, "figs_new/meanplots")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
for (i in unique(out_meanplots$measure))
    dir.create(file.path(out_path, i), showWarnings = FALSE)

opts <- list(
    dist = out_meanplots %>% pull(dist) %>% unique(),
    bias = out_meanplots %>% pull(bias) %>% unique(),
    large = out_meanplots %>% pull(large) %>% unique(),
    effect = out_meanplots %>% pull(effect) %>% unique() #,
    # heterogeneity = out_meanplots %>% pull(heterogeneity) %>% unique()
)
measure_opts <- out_meanplots %>% pull(measure) %>% unique()

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
        img_data <- out_meanplots %>% filter(!!!filters)
        for (me in measure_opts) { # loop over different measures
            img_data_me <- img_data %>% filter(measure == me)
            ylim <- with(img_data_me, c(min(value), max(value)))
            lapply(
                current_levels,
                function(z) {
                    title <- img_data_me %>%
                        select(
                            -k, -I2, -method, -measure,
                            -value, -all_of(current),
                            -sampleSize, -usage, -stat_fun,
                            -is_ci, -is_pi, -is_new
                        ) %>%
                        distinct() %>%
                        bind_cols(!!current := z) %>%
                        select(order(colnames(.))) %>%
                        make_title()
                    p <- img_data_me %>%
                        filter(!!sym(current) == z) %>%
                        plotPanels(measure = me, by = "method") +
                        ylim(ylim) +
                        ggtitle(eval(title)) +
                        theme(
                            text = element_text(size = 9),
                            plot.title = element_text(size = 10)
                        )
                    if (0 > ylim[1] & 0 < ylim[2]) {
                        p <- p +
                            geom_hline(yintercept = 0, lty = 2, alpha = 0.5)
                    }
                    p
                }
            ) %>%
                wrap_plots(guides = "collect") &
                theme(legend.position = "bottom")
            filename <- paste0(
                out_path, "/", me, "/", toupper(current), "_",
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
}

## gamma_max summary statistics ------------------------------------------------

out_path <- file.path(out_dir, "figs_new/max_pH")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

opts <- list(
    dist = unique(out_gammamax$dist),
    bias = unique(out_gammamax$bias),
    large = unique(out_gammamax$large),
    # heterogeneity =  unique(out_gammamax$heterogeneity),
    effect = unique(out_gammamax$effect),
    I2 =  unique(out_gammamax$I2)
)
measure_opts <- unique(out_gammamax$measure)

list_seq <- seq_along(opts)

# Set order of the plots
ord <- c(
    "Edgington (none) CI",
    "Edgington (REML) CI",
    "Edgington (DL) CI",
    "Edgington (PM) CI",
    "Fisher (none) CI",
    "Fisher (REML) CI",
    "Fisher (DL) CI",
    "Fisher (PM) CI",
    "Random effects (none) CI",
    "Random effects (REML) CI",
    "Random effects (DL) CI",
    "Random effects (PM) CI"
)

## For now, subset to only those methods we actually want plots of
# out_gammamax <- subset(out_gammamax, method %in% ord)



for (x in list_seq) { # loop over summary (eg. dist)
    current <- names(opts)[x]
    current_levels <- opts[[current]]
    cat(
        "Currently constructing plots for:",
        current,
        paste0("(", paste0(current_levels, collapse = ", "), ")"),
        "\n"
    )
    grid_others <- expand.grid(
        opts[list_seq[list_seq != x]],
        stringsAsFactors = FALSE
    )
    # loop over all combinations of other parameters (e.g. bias, large, he)
    for (y in seq_len(nrow(grid_others))) {
        # filter the data by current parameter combination
        filters <- lapply(seq_along(grid_others), function(z) {
            lhs <- names(grid_others)[z]
            op <- quote(`==`)
            rhs <- grid_others[y, z]
            expr(`!!`(op)(!!sym(lhs), !!rhs))
        })
        img_data <- out_gammamax %>%
            filter(!!!filters) %>%
            mutate(method = factor(method, levels = ord))
        lapply(
            current_levels,
            function(z) {
                title <- img_data %>%
                    filter(!!sym(current) == z) %>%
                    select(
                        -k, -method, -measure,
                        -value, -all_of(current),
                        -sampleSize, -usage, -stat_fun,
                        -is_pi, -is_ci, -is_new
                    ) %>%
                    distinct() %>%
                    bind_cols(!!current := z) %>%
                    select(order(colnames(.))) %>%
                    make_title()
                img_data %>%
                    filter(!!sym(current) == z) %>%
                    ggplot(aes(x = k, y = value, color = measure)) +
                    facet_wrap(~method, ncol = 4) +
                    geom_point() +
                    stat_summary(fun = "mean", geom = "line", aes(group = measure)) +
                    ylim(c(0, 1)) +
                    labs(
                        x = "# studies",
                        y = "highest p-value",
                        color = "Summary statistic",
                        title = eval(title)
                    )
            }
        ) %>%
            wrap_plots(guides = "collect") &
            theme(legend.position = "bottom",
                  text = element_text(size = 9),
                  plot.title = element_text(size = 10))
        filename <- paste0(
            out_path, "/", toupper(current), "_",
            paste0(
                paste0(names(grid_others[y, ]), "_", to_char(grid_others[y, ])),
                collapse = "_"
            ),
            ".png"
        )
        cat("\33[2K\rMaking file: ", filename)
        ggsave(
            filename = filename,
            width = length(current_levels) * 7,
            height = 12,
            units = "in",
            device = ragg::agg_png
        )
    }
}


## relative frequencies --------------------------------------------------------

out_path <- file.path(out_dir, "figs_new/rel_freq")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# Adapt make_title function for now
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
    if (any(nms == "Large")) nms[nms == "Large"] <- "# large studies"
    if (any(nms == "K")) nms[nms == "K"] <- "# studies"
    if (any(nms == "Heterogeneity")) {
        nms[nms == "Heterogeneity"] <- "Simulation model"
    }
    vals <- unname(unlist(df[1, ]))
    args <- ifelse(
        grepl("theta", nms[1]),
        paste0(
            "bquote(", paste0(paste0(nms, ": ", vals), collapse = ", "),
            "\")"
        ),
        paste0(
            "bquote(\"", paste0(paste0(nms, ": ", vals), collapse = ", "),
            "\")"
        )
    )
    return(eval(parse(text = args)))
}

# These are the options (i.e. what variables to compare and their levels)
opts <- list(
    dist = unique(out_n$dist),
    bias = unique(out_n$bias),
    large = unique(out_n$large),
    effect = unique(out_n$effect),
    # heterogeneity = unique(out_n$heterogeneity),
    I2 = unique(out_n$I2),
    k = unique(out_n$k)
)

# What measures do we have here (this ix the x axis)
measure_opts <- unique(out_n$measure)

list_seq <- seq_along(opts)

# Set order of the plots
ord <- c(
    "Edgington (none) CI",
    "Edgington (REML) CI",
    "Edgington (DL) CI",
    "Edgington (PM) CI",
    "Fisher (none) CI",
    "Fisher (REML) CI",
    "Fisher (DL) CI",
    "Fisher (PM) CI"
)

## For now, subset to only those methods we actually want plots of
# out_n <- subset(out_n, method %in% ord)


# loop over summary (eg. dist)
for (x in list_seq) {
    current <- names(opts)[x]
    current_levels <- opts[[current]]
    cat(
        "Currently constructing plots for:",
        current,
        paste0("(", paste0(current_levels, collapse = ", "), ")"),
        "\n"
    )
    grid_others <- expand.grid(
        opts[list_seq[list_seq != x]],
        stringsAsFactors = FALSE
    )
    # loop over all combinations of other parameters (e.g. bias, large, he)
    for (y in seq_len(nrow(grid_others))) {
        # filter the data by current parameter combination
        filters <- lapply(seq_along(grid_others), function(z) {
            lhs <- names(grid_others)[z]
            op <- quote(`==`)
            rhs <- grid_others[y, z]
            expr(`!!`(op)(!!sym(lhs), !!rhs))
        })
        img_data <- out_n %>%
            filter(!!!filters) %>%
            mutate(
                method = factor(method, levels = ord)
            )
        lapply(current_levels, function(z) {
            title <- img_data %>%
                filter(!!sym(current) == z) %>%
                select(
                    -method, -measure, -value,
                    -all_of(current), -sampleSize,
                    -usage, -stat_fun, -is_ci, -is_pi,
                    -is_new
                ) %>%
                distinct() %>%
                bind_cols(!!current := z) %>%
                select(order(colnames(.))) %>%
                make_title()
            img_data %>%
                filter(!!sym(current) == z) %>%
                ggplot(aes(x = measure, y = value)) +
                geom_col(fill = viridisLite::viridis(1)) +
                ylim(c(0, 1)) +
                facet_wrap(~method, ncol = 4) +
                labs(
                    x = "# intervals",
                    y = "relative frequency",
                    title = eval(title)
                )
        }) %>%
            wrap_plots(guides = "collect") &
            theme(
                legend.position = "bottom",
                text = element_text(size = 9),
                plot.title = element_text(size = 10)
            )
        filename <- paste0(
            out_path, "/", toupper(current), "_",
            paste0(
                paste0(
                    names(grid_others[y, ]), "_",
                    to_char(grid_others[y, ])
                ),
                collapse = "_"
            ),
            ".png"
        )
        cat("\33[2K\rMaking file: ", filename)
        ggsave(
            filename = filename,
            width = length(current_levels) * 7.5,
            height = 12,
            units = "in",
            device = ragg::agg_png
        )
    }
}

# p_accept plots ---------------------------------------------------------------

# delete the following once the simulation finished
out_pubbias <- out_pubbias |>
    filter(value != 1)

out_path <- file.path(out_dir, "figs_new/p_accept")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
for (i in unique(out_pubbias$measure))
    dir.create(file.path(out_path, i), showWarnings = FALSE)

opts <- list(
    dist = out_pubbias %>% pull(dist) %>% unique(),
    bias = out_pubbias %>% pull(bias) %>% unique(),
    large = out_pubbias %>% pull(large) %>% unique(),
    effect = out_pubbias %>% pull(effect) %>% unique() #,
    # heterogeneity = out_pubbias %>% pull(heterogeneity) %>% unique()
)
measure_opts <- out_pubbias %>% pull(measure) %>% unique()

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
        # for (me in measure_opts) { # loop over different measures
            img_data_me <- img_data %>% filter(measure == me)
            ylim <- with(img_data_me, c(min(value), max(value)))
            lapply(
                current_levels,
                function(z) {
                    title <- img_data_me %>%
                        select(
                            -k, -I2, -method, -measure,
                            -value, -all_of(current),
                            -sampleSize, -usage, -stat_fun,
                            -is_ci, -is_pi, -is_new
                        ) %>%
                        distinct() %>%
                        bind_cols(!!current := z) %>%
                        select(order(colnames(.))) %>%
                        make_title()
                    img_data_me %>%
                        filter(!!sym(current) == z) %>%
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
            ) %>%
                wrap_plots(guides = "collect") &
                theme(legend.position = "bottom")
            filename <- paste0(
                out_path, "/", me, "/", toupper(current), "_",
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
        # }
    }
}


# Summary of meanplots ---------------------------------------------------------

# TODO: edit filename (effect is not incorporated yet)

# # Make directories
# out_path <- file.path(out_dir, "figs_new/summary/")
# dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
#
# # Set types of plots
# opts <- list(
#     bias = unique(out_sum$bias),
#     meth = unique(out_sum$method),
#     meas = unique(out_sum$measure)
# )
#
# # Set the order of plots
# method_order <- c(
#     "Edgington (none) CI",
#     "Edgington (REML) CI",
#     "Edgington (DL) CI",
#     "Edgington (PM) CI",
#     "Fisher (none) CI",
#     "Fisher (REML) CI",
#     "Fisher (DL) CI",
#     "Fisher (PM) CI",
#     "Random effects (none) CI",
#     "Random effects (REML) CI",
#     "Random effects (DL) CI",
#     "Random effects (PM) CI",
#     "Hartung & Knapp (REML) CI",
#     "Henmi & Copas (DL) CI"
# )
#
# out_sum <- subset(out_sum, method %in% method_order)
#
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
#
# for (x in opts$meas) {
#     data <- out_sum %>%
#         filter(measure == x) %>%
#         group_by(k, bias, I2, method, effect, heterogeneity) %>%
#         summarise(
#             mean_value = mean(value),
#             max_value = max(value),
#             min_value = min(value),
#             .groups = "drop"
#         )
#     # split methods into additive, multiplicative, and rest
#     methods <- unique(data[c("method", "heterogeneity")])
#     add_idx <- data$heterogeneity == "additive"
#     mult_idx <- grepl("Multiplicative", methods, fixed = TRUE)
#     additive <- methods[add_idx]
#     multiplicative <- methods[mult_idx]
#     rest <- methods[!(add_idx | mult_idx)]
#     # list of methods for each plot
#     method_list <- list(
#         additive = c(additive, rest),
#         multiplicative = c(multiplicative, rest)
#     )
#     # subset data according to methods and fix order
#     lapply(seq_along(method_list), function(y) {
#         out <- subset(data, method %in% method_list[[y]])
#         # Get the ylim
#         ylimes <- c(min(out$min_value), max(out$max_value))
#         # set transparency
#         transparency <- 0.6
#         # get methods to loop over
#         methods <- unique(out$method)
#         # order them
#         methods <- method_order[method_order %in% methods]
#         # make plots
#         plots <- lapply(seq_along(methods), function(z) {
#             out %>%
#                 filter(method == methods[z]) %>%
#                 mutate(
#                     k = factor(k),
#                     I2 = factor(I2),
#                     bias = factor(
#                         bias,
#                         levels = c("none", "moderate", "strong")
#                     )
#                 ) %>%
#                 ggplot(aes(x = k, y = mean_value, color = I2, group = I2)) +
#                 geom_line() +
#                 geom_point() +
#                 geom_errorbar(
#                     aes(ymin = min_value, ymax = max_value),
#                     alpha = transparency
#                 ) +
#                 scale_color_discrete(name = expression(I^2)) +
#                 facet_wrap(~bias) +
#                 ylim(ylimes) +
#                 theme_bw() +
#                 theme(legend.position = "bottom") +
#                 labs(y = methods[z])
#         }) %>%
#             wrap_plots(., guides = "collect", nrow = length(.)) +
#             plot_annotation(title = totitle(x)) &
#             theme(legend.position = "bottom")
#         ggsave(
#             filename = paste0(
#                 out_path,
#                 x,
#                 "_",
#                 if (y == 1L) "additive" else "multiplicative",
#                 ".png"
#             ),
#             device = ragg::agg_png,
#             plot = plots,
#             width = length(unique(out$bias)) * 5,
#             height = length(unique(out$method)) * 5,
#             units = "in"
#         )
#
#     })
# }
