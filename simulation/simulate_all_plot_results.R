## create figures from simulation output

rm(list = ls())
library(tidyverse); theme_set(theme_bw())

load(paste0("RData/simulate_all.RData"))

#' Helper function to plot the results
#' @param data obtained by simulate_all.R and saved in RData/simulate_all.RData
#' @param measure CI measure to plot
#' @param by make facets based on "method" or "I2".
plotPanels <- function(data,
                       measure = c("width", "coverage", "coverage_effects",
                                   "score", "n"),
                       by = c("method", "I2")){
    measure <- match.arg(measure)
    by <- match.arg(by)
    data[data$measure == measure, ] %>%
        mutate(k = factor(k),
               I2 = factor(I2)) -> data
    if(measure == "n")
        data <- data[data$method %in% c("Harmonic Mean", "Harmonic Mean Additive",
                                        "Harmonic Mean Multiplicative"), ]
    if(by == "method")
        data %>% mutate(method = paste(method)) %>%
            ggplot(mapping = aes(x = k, y = value, color = I2)) +
            facet_wrap(~ method) +
            geom_point() +
            stat_summary(fun="mean", geom="line", aes(group=I2)) +
            scale_color_discrete(name = expression(I^2)) +
            xlab("# studies") +
            ylab(measure) -> p        
    if(by == "I2")
        data %>% mutate(I2 = as.character(I2)) %>%
        ggplot(mapping = aes(x = k, y = value, color = method)) +
            facet_wrap(~ I2, labeller = label_bquote(I^2 == .(I2))) +
            geom_point() +
            stat_summary(fun="mean", geom="line", aes(group=method)) +
            xlab("# studies") +
            ylab(measure) -> p       

    if(str_detect(measure, "coverage")) {
        p <- p + 
            geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5)
    }
    p
}

## sim %>% 
##     filter(measure == "coverage_new_mean") %>%
##     arrange(method)




## 4.1
dir.create("figs", showWarnings = FALSE)
out %>%
    filter(method != "Harmonic Mean two sided") %>%
    mutate(method = ifelse(method == "REML", "Random Effects, default, REML", method),
           measure = gsub("_mean", "", measure))-> out2

out2 %>% select(measure) %>% unique() %>% pull() -> measure
out2 %>% select(dist) %>% unique() %>% pull() -> dist 
out2 %>% select(bias) %>% unique() %>% pull() -> bias 
out2 %>% select(large) %>% unique() %>% pull() -> large 
out2 %>% select(heterogeneity) %>% unique() %>% pull() -> heterogeneity 

grid <- expand.grid(measure = measure, dist = dist, bias = bias, large = large, 
                    heterogeneity = heterogeneity, stringsAsFactors = FALSE) %>% 
    filter(!(measure == "coverage_effects" & heterogeneity == "multiplicative"))




for(i in 1:nrow(grid)){
        print(grid[i,])
        out2[out2$dist == grid[i, "dist"] &
             out2$measure == grid[i, "measure"] &
             out2$bias == grid[i, "bias"] &
             out2$large == grid[i, "large"] &
             out2$heterogeneity == grid[i, "heterogeneity"], ] %>%
            plotPanels(data=., measure = grid[i, "measure"], by="method") +
            ggtitle(paste0("dist: ", grid[i, "dist"],
                           ", bias: ", grid[i, "bias"],
                           ", theta: 0.2, no. of large studies: ", grid[i, "large"],
                           ", simulation model: ", grid[i, "heterogeneity"], sep=""))
        ggsave(filename = paste0("figs/simulate_all",
                                 "_", grid[i, "dist"],
                                 "_large_", grid[i, "large"],
                                 "_bias_", grid[i, "bias"],
                                 "_sim-mod_", grid[i, "heterogeneity"], 
                                 "_", grid[i, "measure"], ".png"),
               width = 7,
               height = 6,
               units = "in")
}


## make summary pannels
## needs convert from ImageMagick-ims6.q16(1) available for Linux Ubuntu and other systems
## 1. bias none vs strong
for(di in dist)
    for(la in large)
        for(me in measure)
            for(he in heterogeneity)
                if(!(he == "multiplicative" & me == "coverage_effects")){
                    system(paste0("convert ",
                                  "figs/simulate_all",
                                  "_", di,
                                  "_large_", la,
                                  "_bias_none",
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/simulate_all",
                                  "_", di,
                                  "_large_", la,
                                  "_bias_moderate",
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/simulate_all",
                                  "_", di,
                                  "_large_", la,
                                  "_bias_strong",
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " +append ",
                                  "figs/simulate_all_summary_panel_BIAS",
                                  "_", di,
                                  "_large_", la,
                                  "_sim-mod_", he,
                                  "_", me, ".png"))
                }
                

## 2. no large trials 0, 1, 2
for(di in dist)
    for(bi in bias)
        for(me in measure)
            for(he in heterogeneity)
                if(!(he == "multiplicative" & me == "coverage_effects")){
                    system(paste0("convert ",
                                  "figs/simulate_all",
                                  "_", di,
                                  "_large_0",
                                  "_bias_", bi,
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/simulate_all",
                                  "_", di,
                                  "_large_1",
                                  "_bias_", bi, 
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/simulate_all",
                                  "_", di,
                                  "_large_2",
                                  "_bias_", bi,
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  " +append ",
                                  "figs/simulate_all_summary_panel_TRIAL_SIZE",
                                  "_", di,
                                  "_bias_", bi,
                                  "_sim-mod_", he,
                                  "_", me, ".png"))
                }

## 3. heterogeneity model used for simulation
for(di in dist)
    for(bi in bias)
        for(me in measure)
            for(la in large)
                if(me != "coverage_effects"){
                    system(paste0("convert ",
                                  "figs/simulate_all",
                                  "_", di,
                                  "_large_", la,
                                  "_bias_", bi,
                                  "_sim-mod_additive",
                                  "_", me, ".png",
                                  " ",
                                  "figs/simulate_all",
                                  "_", di,
                                  "_large_", la,
                                  "_bias_", bi, 
                                  "_sim-mod_multiplicative",
                                  "_", me, ".png",
                                  " ",
                                  " +append ",
                                  "figs/simulate_all_summary_panel_SIM-MOD",
                                  "_", di,
                                  "_bias_", bi,
                                  "_sim-mod_", he,
                                  "_", me, ".png"))
                }
