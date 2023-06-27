## create figures from simulation output

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

# Load data ----
load(paste0("RData/simulate_all.RData"))

# Prepare data (meanplots) ----
out_meanplots <- out %>%
    bind_rows() %>%
    filter(
        grepl("_mean$", measure),
        !grepl("^gammaMin", measure),
        !grepl("two sided", method),
        !grepl("\\(f\\)", method)
    ) %>%
    mutate(
        method = ifelse(
            grepl("REML", method),
            gsub("REML", "Random Effects, REML", method),
            method
        ),
        measure = gsub("_mean$", "", measure)
    )

# Prepare data (gammaMin) ----
out_gammamin <- out %>%
    bind_rows() %>%
    filter(
        grepl("^gammaMin", measure),
        grepl("Harmonic Mean|k-Trials|Pearson", method),
        !grepl("two sided", method),
        !grepl("\\(f\\)", method)
    ) %>%
    distinct() %>% # throw out the doubly included gammaMin_mean entries
    mutate(
        measure = gsub("^gammaMin_", "", measure),
        measure = case_when(
            measure == "min" ~ "Minimum",
            measure == "firstQuart" ~ "1. Quartile",
            measure == "mean" ~ "Mean",
            measure == "median" ~ "Median",
            measure == "thirdQuart" ~ "3. Quartile",
            measure == "max" ~ "Maximum"
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
    bind_rows() %>%
    filter(
        grepl("^n", measure),
        !grepl("_mean$", measure),
        !grepl("\\(f\\)", method)
    ) %>%
    mutate(
        measure = gsub("n_", "", measure),
        measure = gsub("gt", "> ", measure),
        measure = ordered(measure, levels = rev(c("> 9", as.character(9:1)))),
        value = value / 5000
    )

# Prepare data (summary) ----
out_sum <- out_meanplots %>%
    filter(grepl("coverage", measure),
           grepl(
               # "Harmonic Mean|k-Trials|Hartung|Henmy|REML",
               "Harmonic Mean|k-Trials|Pearson|Hartung|REML",
               method
           ),
           !grepl("\\(f\\)", method))

# Set output directory -----
out_dir <- "~/test/Institution/harmonic_mean"

#' #' Helper function to plot means
#' #' @param data obtained by simulate_all.R and saved in RData/simulate_all.RData
#' #' @param measure CI measure to plot
#' #' @param by make facets based on "method" or "I2".
plotPanels <- function(
        data,
        measure = c(
            "coverage_true", "coverage_effects", "coverage_effects_all",
            "coverage_effects_min1", "coverage_prediction", "n",
            "width", "score"
        ),
        by = c("method", "I2")
    ) {
    measure <- match.arg(measure)
    by <- match.arg(by)
    data2 <- data[data$measure == measure, ] %>%
        mutate(
            k = factor(k),
            I2 = factor(I2)
        )
    if (measure == "n") {
        data2 <- data2[grepl("Harmonic Mean|k-Trials|Pearson", data$method), ]
    }
    if (by == "method") {
        p <- data2 %>%
            # Set order of plots
            {
                if (measure == "coverage_prediction") {
                    .[] %>%
                        filter(grepl("Harmonic Mean|k-Trials|Pearson|PI", method)) %>%
                        mutate(
                            method = factor(
                                method,
                                levels = c(
                                    "Harmonic Mean CI (chisq)",
                                    "Harmonic Mean Additive CI (chisq)",
                                    "Harmonic Mean Multiplicative CI (chisq)",
                                    "Hartung & Knapp PI",
                                    # "Harmonic Mean CI (f)",
                                    # "Harmonic Mean Additive CI (f)",
                                    # "Harmonic Mean Multiplicative CI (f)",
                                    "k-Trials CI",
                                    "k-Trials Additive CI",
                                    "k-Trials Multiplicative CI",
                                    "Random Effects, REML PI",
                                    "Pearson CI",
                                    "Pearson Additive CI",
                                    "Pearson Multiplicative CI"
                                )
                            )
                        )
                } else {
                    .[] %>%
                        filter(!grepl("PI", method)) %>%
                        mutate(
                            method = factor(
                                method,
                                levels = c(
                                    "Harmonic Mean CI (chisq)",
                                    "Harmonic Mean Additive CI (chisq)",
                                    "Harmonic Mean Multiplicative CI (chisq)",
                                    "Hartung & Knapp CI",
                                    # "Harmonic Mean CI (f)",
                                    # "Harmonic Mean Additive CI (f)",
                                    # "Harmonic Mean Multiplicative CI (f)",
                                    "k-Trials CI",
                                    "k-Trials Additive CI",
                                    "k-Trials Multiplicative CI",
                                    "Henmy & Copas CI",
                                    "Pearson CI",
                                    "Pearson Additive CI",
                                    "Pearson Multiplicative CI",
                                    "Random Effects, REML CI"
                                )
                            )
                        )
                }
            } %>%
            ggplot(mapping = aes(x = k, y = value, color = I2)) +
            facet_wrap(~ method) +
            geom_point() +
            stat_summary(fun="mean", geom="line", aes(group=I2)) +
            scale_color_discrete(name = expression(I^2)) +
            xlab("# studies") +
            ylab(measure)
    }
    if (by == "I2") {
        p <- data2 %>% mutate(I2 = as.character(I2)) %>%
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

# filter_data <- function(measure, plot_type){
#     filter_data_meanplots <- function(measure) {
#         
#     }
#     filter_data_gamma <- function(measure) {
#         
#     }
#     filter_data_frequency <- function(measure) {
#         
#     }
#     filter_data_summary <- function(measure) {
#         
#     }
# }

## Mean plots ---------------------------------------------------------------------



## 4.1
out_path <- file.path(out_dir, "figs/meanplots")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
for (i in unique(out_meanplots$measure))
    dir.create(file.path(out_path, i), showWarnings = FALSE)

opts <- list(
    dist = out_meanplots %>% pull(dist) %>% unique(),
    bias = out_meanplots %>% pull(bias) %>% unique(),
    large = out_meanplots %>% pull(large) %>% unique(),
    heterogeneity = out_meanplots %>% pull(heterogeneity) %>% unique()
)
measure_opts <- out_meanplots %>% pull(measure) %>% unique()

list_seq <- seq_along(opts)

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
        img_data <- out_meanplots %>% filter(!!!filters)
        for (me in measure_opts) { # loop over different measures
            img_data %>%
                filter(measure == me) %>%
                {
                    if (me == "coverage_prediction") {
                        ylim <- .[] %>%
                            filter(
                                grepl("Harmonic Mean|PI|k-Trials", method)
                            ) %>%
                            pull(value) %>%
                            {
                                c(min(.), max(.))
                            }
                    } else {
                        ylim <- .[] %>%
                            filter(!grepl("PI", method)) %>%
                            pull(value) %>%
                            {
                                c(min(.), max(.))
                            }
                    }
                    lapply(current_levels, function(z) {
                        title <- .[] %>%
                            select(
                                -k, -I2, -method, -measure,
                                -value, -all_of(current),
                                -sampleSize
                            ) %>%
                            distinct() %>%
                            bind_cols(!!current := z) %>%
                            select(order(colnames(.))) %>%
                            make_title()
                        .[] %>%
                            filter(!!sym(current) == z) %>%
                            plotPanels(measure = me, by = "method") +
                            ylim(ylim) +
                            ggtitle(eval(title)) +
                            theme(
                                text = element_text(size = 9),
                                plot.title = element_text(size = 10)
                            )
                    })
                } %>%
                wrap_plots(guides = "collect") &
                theme(legend.position = "bottom")
            filename <- paste0(
                out_path, "/", me, "/", toupper(current), "_",
                paste0(
                    paste0(
                        names(grid_others[y, ]),
                        "_", grid_others[y, ]
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

## gamma_min summary statistics ------------------------------------------------

out_path <- file.path(out_dir, "figs/min_pH")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

opts <- list(
    dist = out_gammamin %>% pull(dist) %>% unique(),
    bias = out_gammamin %>% pull(bias) %>% unique(),
    large = out_gammamin %>% pull(large) %>% unique(),
    heterogeneity = out_gammamin %>% pull(heterogeneity) %>% unique(),
    I2 = out_gammamin %>% pull(I2) %>% unique()
)
measure_opts <- out_gammamin %>% pull(measure) %>% unique()

list_seq <- seq_along(opts)

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
        img_data <- out_gammamin %>%
            filter(!!!filters) %>%
            mutate(
                method = factor(
                    method,
                    levels = c(
                        "Harmonic Mean CI (chisq)",
                        "Harmonic Mean Additive CI (chisq)",
                        "Harmonic Mean Multiplicative CI (chisq)",
                        #"Harmonic Mean CI (f)",
                        #"Harmonic Mean Additive CI (f)",
                        #"Harmonic Mean Multiplicative CI (f)",
                        "k-Trials CI",
                        "k-Trials Additive CI",
                        "k-Trials Multiplicative CI",
                        "Pearson CI",
                        "Pearson Additive CI",
                        "Pearson Multiplicative CI"
                    )
                )
            )
        lapply(current_levels, function(z) {
            title <- img_data %>%
                filter(!!sym(current) == z) %>%
                select(
                    -k, -method, -measure,
                    -value, -all_of(current),
                    -sampleSize
                    ) %>%
                distinct() %>%
                bind_cols(!!current := z) %>%
                select(order(colnames(.))) %>%
                make_title()
            img_data %>%
                filter(!!sym(current) == z) %>%
                ggplot(aes(x = k, y = value, color = measure)) +
                geom_point() +
                geom_line() +
                ylim(c(0, 1)) +
                facet_wrap(~method) +
                labs(x = "# studies", y = "value",
                     color = "Summary statistic",
                     title = eval(title))
        }) %>%
            wrap_plots(guides = "collect") &
            theme(legend.position = "bottom",
                  text = element_text(size = 9),
                  plot.title = element_text(size = 10))
        ggsave(filename = paste0(
            out_path, "/", toupper(current), "_",
            paste0(
                paste0(names(grid_others[y, ]), "_", grid_others[y, ]),
                collapse = "_"
            ),
            ".png"),
               width = length(current_levels) * 7,
               height = 12,
               units = "in",
               device = ragg::agg_png)
    }
}


## relative frequencies --------------------------------------------------------

out_path <- file.path(out_dir, "figs/rel_freq/")
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

opts <- list(
    dist = out_n %>% pull(dist) %>% unique(),
    bias = out_n %>% pull(bias) %>% unique(),
    large = out_n %>% pull(large) %>% unique(),
    heterogeneity = out_n %>% pull(heterogeneity) %>% unique(),
    I2 = out_n %>% pull(I2) %>% unique(),
    k = out_n %>% pull(k) %>% unique()
)
measure_opts <- out_n %>% pull(measure) %>% unique()

list_seq <- seq_along(opts)

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
                method = factor(
                    method,
                    levels = c(
                        "Harmonic Mean CI (chisq)",
                        "Harmonic Mean Additive CI (chisq)",
                        "Harmonic Mean Multiplicative CI (chisq)",
                        #"Harmonic Mean CI (f)",
                        #"Harmonic Mean Additive CI (f)",
                        #"Harmonic Mean Multiplicative CI (f)",
                        "k-Trials CI",
                        "k-Trials Additive CI",
                        "k-Trials Multiplicative CI",
                        "Pearson CI",
                        "Pearson Additive CI",
                        "Pearson Multiplicative CI"
                    )
                )
            )
        lapply(current_levels, function(z) {
            title <- img_data %>%
                filter(!!sym(current) == z) %>%
                select(
                    -method, -measure, -value,
                    -all_of(current), -sampleSize
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
                facet_wrap(~method) +
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
        ggsave(
            filename = paste0(
                out_path, "/", toupper(current), "_",
                paste0(
                    paste0(
                        names(grid_others[y, ]), "_",
                        grid_others[y, ]), collapse = "_"
                ),
                ".png"
            ),
            width = length(current_levels) * 7.5,
            height = 12,
            units = "in",
            device = ragg::agg_png
        )
    }
}

# Summary of meanplots ---------------------------------------------------------

out_path <- file.path(out_dir, "figs/summary/")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

opts <- list(bias = out_sum %>% pull(bias) %>% unique(),
             meth = out_sum %>% pull(method) %>% unique(),
             meas = out_sum %>% pull(measure) %>% unique())

totitle <- function(string) {
    stopifnot(
        is.character(string),
        length(string) == 1L
    )
    string <- gsub("_", " ", string)
    string <- unlist(strsplit(string, split = " "))
    string <- strsplit(string, split = "")
    string <- vapply(string, function(y) {
        y[1] <- toupper(y[1])
        paste0(y, collapse = "")
    }, character(1L))
    paste0(string, collapse = " ")
}

for (x in opts$meas) {
    data <- out_sum %>%
        filter(measure == x) %>%
        {
            if (x == "coverage_prediction") {
                .[] %>%
                    filter(
                        grepl(
                            "(^Harmonic Mean .+ CI .+|^.+PI$|^k-Trials|^Pearson.*)",
                            method
                        )
                    )
            } else {
                .[] %>%
                    filter(!grepl("^.+PI$", method))
            }
        } %>%
        group_by(k, bias, I2, method) %>%
        summarise(
            mean_value = mean(value),
            max_value = max(value),
            min_value = min(value),
            .groups = "drop"
        )
    # split methods into additive, multiplicative, and rest
    methods <- data %>% pull(method) %>% unique()
    additive <- grep("Additive", methods, value = TRUE)
    multiplicative <- grep("Multiplicative", methods, value = TRUE)
    rest <- grep(
        "^((?!Harmonic|k-Trials|Pearson).)*$",
        methods,
        value = TRUE,
        perl = TRUE
    )
    # list of methods for each plot
    method_list <- list(
        additive = c(additive, rest),
        multiplicative = c(multiplicative, rest)
    )
    # subset data according to methods and fix order
    lapply(seq_along(method_list), function(y) {
        out <- subset(data, method %in% method_list[[y]])
        # Get the ylim
        ylimes <- c(min(out$min_value), max(out$max_value))
        # set transparency
        transparency <- 0.6
        # get methods to loop over
        methods <- unique(out$method)
        # order them
        method_order <- c(
            "Harmonic Mean CI (chisq)",
            "Harmonic Mean Additive CI (chisq)",
            "Harmonic Mean Multiplicative CI (chisq)",
            # "Harmonic Mean CI (f)",
            # "Harmonic Mean Additive CI (f)",
            # "Harmonic Mean Multiplicative CI (f)",
            "k-Trials CI",
            "k-Trials Additive CI",
            "k-Trials Multiplicative CI",
            "Pearson CI",
            "Pearson Additive CI",
            "Pearson Multiplicative CI",
            "Random Effects, REML CI",
            "Random Effects, REML PI",
            "Hartung & Knapp CI",
            "Hartung & Knapp PI",
            "Henmy & Copas CI"
        )
        methods <- method_order[method_order %in% methods]
        # make plots
        plots <- lapply(seq_along(methods), function(z) {
            out %>%
                filter(method == methods[z]) %>%
                mutate(
                    k = factor(k),
                    I2 = factor(I2),
                    bias = factor(
                        bias,
                        levels = c("none", "moderate", "strong")
                    )
                ) %>%
                ggplot(aes(x = k, y = mean_value, color = I2, group = I2)) +
                geom_line() +
                geom_point() +
                geom_errorbar(
                    aes(ymin = min_value, ymax = max_value),
                    alpha = transparency
                ) +
                scale_color_discrete(name = expression(I^2)) +
                facet_wrap(~bias) +
                ylim(ylimes) +
                theme_bw() +
                theme(legend.position = "bottom") +
                labs(y = methods[z])
        }) %>%
            wrap_plots(., guides = "collect", nrow = length(.)) +
            plot_annotation(title = totitle(x)) &
            theme(legend.position = "bottom")
        ggsave(
            filename = paste0(
                out_path, x, "_",
                if (y == 1L) "additive" else "multiplicative",
                ".png"
            ),
            device = ragg::agg_png,
            plot = plots,
            width = length(unique(out$bias)) * 5,
            height = length(unique(out$method)) * 5,
            units = "in"
        )

    })
}

