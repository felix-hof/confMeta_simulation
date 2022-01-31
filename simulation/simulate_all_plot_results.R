## create figures from simulation output

# load libraries
rm(list = ls())
library(rlang)
library(tidyverse); theme_set(theme_bw())
library(patchwork)

# helper function to construct plot titles ----
make_title <- function(df){
    nms <- names(df)
    nms <- vapply(strsplit(nms, ""), 
                  function(x){x[1] <- toupper(x[1]); paste0(x, collapse = "")},
                  character(1L))
    if(any(grepl("Effect", nms))){
        nms[grepl("Effect", nms)] <- ifelse(grep("Effect", nms) == 1, "theta~\"", "\"~theta~\"")
    }
    if(any(nms == "Large")){nms[nms == "Large"] <- "No. of large studies"}
    if(any(nms == "Heterogeneity")){nms[nms == "Heterogeneity"] <- "Simulation model"}
    vals <- unname(unlist(df[1, ]))
    args <- ifelse(grepl("theta", nms[1]), 
                   paste0("bquote(", paste0(paste0(nms, ": ", vals), collapse = ", "), "\")"),
                   paste0("bquote(\"", paste0(paste0(nms, ": ", vals), collapse = ", "), "\")"))
    return(eval(parse(text = args)))
}

# Load data ----
load(paste0("RData/simulate_all.RData"))

# Prepare data (meanplots) ----
out_meanplots <- out %>% 
    bind_rows() %>% 
    filter(grepl("_mean$", measure), 
           !grepl("^gammaMin", measure),
           !grepl("two sided", method)) %>% 
    mutate(method = ifelse(grepl("REML", method), gsub("REML", "Random Effects, default, REML", method), method),
           measure = gsub("_mean$", "", measure))

# Prepare data (gammaMin) ----
out_gammamin <- out %>% 
    bind_rows() %>% 
    filter(grepl("^gammaMin", measure), 
           grepl("Harmonic Mean", method),
           !grepl("two sided", method)) %>% 
    distinct() %>% # throw out the doubly included gammaMin_mean entries
    mutate(measure = gsub("^gammaMin_", "", measure),
           measure = case_when(measure == "min" ~ "Minimum",
                               measure == "firstQuart" ~ "1. Quartile",
                               measure == "mean" ~ "Mean",
                               measure == "median" ~ "Median",
                               measure == "thirdQuart" ~ "3. Quartile",
                               measure == "max" ~ "Maximum"),
           measure = ordered(measure, levels = c("Maximum", "3. Quartile", "Median", "Mean", "1. Quartile", "Minimum"))) %>%
    arrange(k, dist, bias, large, heterogeneity, I2, method, measure)

# Prepare data (frequency) ----
out_n <- out %>% 
    bind_rows() %>% 
    filter(grepl("^n", measure), 
           !grepl("_mean$", measure)) %>% 
    mutate(measure = gsub("n_", "", measure),
           measure = gsub("gt", "> ", measure),
           measure = ordered(measure, levels = rev(c("> 9", as.character(9:1)))),
           value = value/1e4)

#' Helper function to plot means
#' @param data obtained by simulate_all.R and saved in RData/simulate_all.RData
#' @param measure CI measure to plot
#' @param by make facets based on "method" or "I2".
plotPanels <- function(data,
                       measure = c("coverage_true", "coverage_effects",
                                   "coverage_prediction", "n", "width", "score"),
                       by = c("method", "I2")){
    measure <- match.arg(measure)
    by <- match.arg(by)
    data2 <- data[data$measure == measure, ] %>%
        mutate(k = factor(k),
               I2 = factor(I2))
    if(measure == "n"){
        data2 <- data2[grepl("Harmonic Mean", data$method), ]
    }
    if(by == "method"){
        data2 %>% 
            # Set order of plots
            {
                if(measure == "coverage_prediction"){
                    .[] %>% 
                        filter(grepl("Harmonic Mean|PI", method)) %>% 
                        mutate(method = factor(method, 
                                               levels = c("Harmonic Mean CI (chisq)", "Harmonic Mean Additive CI (chisq)",
                                                          "Harmonic Mean Multiplicative CI (chisq)", "Harmonic Mean CI (f)", 
                                                          "Harmonic Mean Additive CI (f)", "Harmonic Mean Multiplicative CI (f)",
                                                          "Hartung & Knapp PI", "Random Effects, default, REML PI")))
                } else {
                    .[] %>%
                        filter(!grepl("PI", method)) %>% 
                        mutate(method = factor(method, 
                                               levels = c("Harmonic Mean CI (chisq)", "Harmonic Mean Additive CI (chisq)",
                                                          "Harmonic Mean Multiplicative CI (chisq)", "Harmonic Mean CI (f)", 
                                                          "Harmonic Mean Additive CI (f)", "Harmonic Mean Multiplicative CI (f)",
                                                          "Hartung & Knapp CI", "Random Effects, default, REML CI", 
                                                          "Henmy & Copas CI")))
                }
            } %>%
            ggplot(mapping = aes(x = k, y = value, color = I2)) +
            facet_wrap(~ method) +
            geom_point() +
            stat_summary(fun="mean", geom="line", aes(group=I2)) +
            scale_color_discrete(name = expression(I^2)) +
            xlab("# studies") +
            ylab(measure) -> p
    }
    if(by == "I2"){
        data2 %>% mutate(I2 = as.character(I2)) %>%
        ggplot(mapping = aes(x = k, y = value, color = method)) +
            facet_wrap(~ I2, labeller = label_bquote(I^2 == .(I2))) +
            geom_point() +
            stat_summary(fun="mean", geom="line", aes(group=method)) +
            xlab("# studies") +
            ylab(measure) -> p
    }

    if(str_detect(measure, "coverage")){
        p <- p +
            geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5)
    }
    p
}


## Mean plots ---------------------------------------------------------------------



## 4.1
dir.create("figs/meanplots", showWarnings = FALSE, recursive = TRUE)

opts <- list(dist = out_meanplots %>% pull(dist) %>% unique(),
             bias = out_meanplots %>% pull(bias) %>% unique(),
             large = out_meanplots %>% pull(large) %>% unique(),
             heterogeneity = out_meanplots %>% pull(heterogeneity) %>% unique())
measure_opts <- out_meanplots %>% pull(measure) %>% unique()

list_seq <- seq_along(opts)

for(x in list_seq){ # loop over summary (eg. dist)
    current <- names(opts)[x]
    current_levels <- opts[[current]]
    cat("Currently constructing plots for:", current, paste0("(", paste0(current_levels, collapse = ", "), ")"), "\n")
    grid_others <- expand.grid(opts[list_seq[list_seq != x]], stringsAsFactors = FALSE)
    for(y in seq_len(nrow(grid_others))){ # loop over all combinations of other parameters (bias, large, he)
        # filter the data by current parameter combination
        filters <- lapply(seq_along(grid_others), function(z){
            lhs <- names(grid_others)[z]
            op <- quote(`==`)
            rhs <- grid_others[y, z]
            expr(`!!`(op)(!!sym(lhs), !!rhs))
        })
        img_data <- out_meanplots %>% filter(!!!filters)
        for(me in measure_opts){ # loop over different measures
            img_data %>% 
                filter(measure == me) %>% 
                {
                    if(me == "coverage_prediction"){
                        ylim <- .[] %>% 
                            filter(grepl("Harmonic Mean|PI", method)) %>% 
                            pull(value) %>% 
                            {c(min(.), max(.))}
                    } else {
                        ylim <- .[] %>%
                            filter(!grepl("PI", method)) %>% 
                            pull(value) %>% 
                            {c(min(.), max(.))}
                    }
                    lapply(current_levels, function(z){
                        title <- .[] %>% 
                            select(-k, -I2, -method, -measure, -value, -all_of(current), -sampleSize) %>% 
                            distinct() %>% 
                            bind_cols(!!current := z) %>% 
                            select(order(colnames(.))) %>% 
                            make_title()
                        .[] %>% 
                            filter(!!sym(current) == z) %>% 
                            plotPanels(measure = me, by="method") +
                            ylim(ylim) +
                            ggtitle(eval(title)) +
                            theme(text = element_text(size = 9),
                                  plot.title = element_text(size = 10))
                    })
                } %>% 
                wrap_plots(guides = "collect") &
                theme(legend.position = "bottom")
            ggsave(filename = paste0("figs/meanplots/", toupper(current), "_", 
                                     paste0(paste0(names(grid_others[y, ]), "_", grid_others[y, ]), collapse = "_"),
                                     "_", me, ".png"),
                   width = length(current_levels) * 6.5,
                   height = 12,
                   units = "in",
                   type = "cairo")
        } 
    }
}

## gamma_min summary statistics -----------------------------------------------------------

dir.create("figs/min_pH", showWarnings = FALSE, recursive = TRUE)

out_gammamin %>% select(dist) %>% unique() %>% pull() -> dist
out_gammamin %>% select(bias) %>% unique() %>% pull() -> bias
out_gammamin %>% select(large) %>% unique() %>% pull() -> large
out_gammamin %>% select(heterogeneity) %>% unique() %>% pull() -> heterogeneity
out_gammamin %>% select(I2) %>% unique() %>% pull() -> I2
#out_gammamin %>% select(measure) %>% unique() %>% pull() -> me

for(di in dist){
    for(la in large){
        for(bi in bias){
            for(he in heterogeneity){
                for(i in I2){
                    pars <- out_gammamin %>%
                        filter(dist == di, bias == bi, large == la, heterogeneity == he, I2 == i) %>%
                        ggplot(aes(x = k, y = value, colour = measure)) +
                        facet_wrap(~ method) +
                        geom_point() +
                        geom_line() +
                        xlab("# studies") +
                        labs(colour = "Statistic") +
                        ggtitle(bquote("dist:"~.(di)~", bias: "~.(bi)~", "~theta~": 0.2, "~I^2~": "~.(i)~", no. of large studies: "~.(la)~", simulation model: "~.(he)))
                    ggsave(filename = paste0("figs/min_pH/MIN-PH_", di, "_bias_", bi, "_large_", la,"_sim-mod_", he, "_I2_", i, ".png"),
                           width = 10,
                           height = 11,
                           units = "in")


                }
            }
        }
    }
}


## relative frequencies -----------------------------------------------------------

dir.create("figs/rel_freq", showWarnings = FALSE, recursive = TRUE)

out_n %>% select(dist) %>% unique() %>% pull() -> dist
out_n %>% select(bias) %>% unique() %>% pull() -> bias
out_n %>% select(large) %>% unique() %>% pull() -> large
out_n %>% select(heterogeneity) %>% unique() %>% pull() -> heterogeneity
out_n %>% select(I2) %>% unique() %>% pull() -> I2
out_n %>% select(k) %>% unique() %>% pull() -> k

for(di in dist){
    for(la in large){
        for(bi in bias){
            for(he in heterogeneity){
                for(i in I2){
                    for(K in k){
                        out_n %>%
                            filter(dist == di, bias == bi, large == la, heterogeneity == he, I2 == i, k == K) %>%
                            ggplot(aes(x = measure, y = value, fill = k)) +
                            facet_wrap(~ method) +
                            geom_col() +
                            xlab("# intervals") +
                            ylab("relative frequency") +
                            theme(legend.position = "none") +
                            ggtitle(bquote("# studies: "~.(K)~", dist:"~.(di)~", bias: "~.(bi)~", "~theta~": 0.2, "~I^2~": "~.(i)~", # large studies: "~.(la)~", simulation model: "~.(he)))
                        ggsave(filename = paste0("figs/rel_freq/k_", K, "_", di, "_bias_", bi, "_large_", la,"_sim-mod_", he, "_I2_", i, ".png"),
                               width = 10,
                               height = 11,
                               units = "in")


                    }
                }
            }
        }
    }
}


for(di in dist){
    for(la in large){
        for(bi in bias){
            for(he in heterogeneity){
                for(i in I2){
                    system(paste0("convert \\( ",
                                  "figs/rel_freq/k_2_",
                                  di,
                                  "_bias_", bi,
                                  "_large_", la,
                                  "_sim-mod_", he,
                                  "_I2_", i,
                                  ".png ",
                                  "figs/rel_freq/k_3_",
                                  di,
                                  "_bias_", bi,
                                  "_large_", la,
                                  "_sim-mod_", he,
                                  "_I2_", i,
                                  ".png ",
                                  "+append \\) ",
                                  "\\( ",
                                  "figs/rel_freq/k_5_",
                                  di,
                                  "_bias_", bi,
                                  "_large_", la,
                                  "_sim-mod_", he,
                                  "_I2_", i,
                                  ".png ",
                                  "figs/rel_freq/k_10_",
                                  di,
                                  "_bias_", bi,
                                  "_large_", la,
                                  "_sim-mod_", he,
                                  "_I2_", i,
                                  ".png ",
                                  "+append \\) ",
                                  "figs/rel_freq/k_20_",
                                  di,
                                  "_bias_", bi,
                                  "_large_", la,
                                  "_sim-mod_", he,
                                  "_I2_", i,
                                  ".png ",
                                  "-append figs/rel_freq/STUDIES_", di, "_bias_", bi,"_large_", la, "_sim-mod_", he, "_I2_", i, ".png "))
                }
            }
        }
    }
}


