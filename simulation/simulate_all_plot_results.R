## create figures from simulation output

rm(list = ls())
library(tidyverse); theme_set(theme_bw())

load(paste0("RData/simulate_all.RData"))

out <- out %>% 
    bind_rows() #%>%
    # mutate(interval = ifelse(grepl("CI", method), "CI", "PI"),
    #        method = gsub("\\s*CI|\\s*PI", "", method))

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
    data[data$measure == measure, ] %>%
        mutate(k = factor(k),
               I2 = factor(I2)) -> data
    if(measure == "n")
        data <- data[grepl("Harmonic Mean", data$method), ]
        # data <- data[data$method %in% c("Harmonic Mean", "Harmonic Mean Additive",
        #                                 "Harmonic Mean Multiplicative"), ]
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


## Mean plots ---------------------------------------------------------------------



## 4.1
dir.create("figs", showWarnings = FALSE)
dir.create("figs/meanplots", showWarnings = FALSE)
out2 <- out %>%
    filter(!grepl("two sided", method),
           grepl("_mean$", measure), !grepl("gammaMin", measure)) %>%
    mutate(method = ifelse(grepl("REML", method), gsub("REML", "Random Effects, default, REML", method), method),
           measure = gsub("_mean$", "", measure))

out2 %>% select(measure) %>% unique() %>% pull() -> measure
out2 %>% select(dist) %>% unique() %>% pull() -> dist 
out2 %>% select(bias) %>% unique() %>% pull() -> bias 
out2 %>% select(large) %>% unique() %>% pull() -> large 
out2 %>% select(heterogeneity) %>% unique() %>% pull() -> heterogeneity
#out2 %>% select(interval) %>% unique() %>% pull() -> interval

grid <- expand.grid(measure = measure, dist = dist, bias = bias, large = large, 
                    heterogeneity = heterogeneity, 
                    #interval = interval, 
                    stringsAsFactors = FALSE)



single_plots <- character(nrow(grid))
for(i in 1:nrow(grid)){
    print(grid[i,])
    single_plots[i] <- paste0("figs/meanplots/",
                          grid[i, "dist"],
                          "_large_", grid[i, "large"],
                          "_bias_", grid[i, "bias"],
                          "_sim-mod_", grid[i, "heterogeneity"],
                          #"_", grid[i, "interval"],
                          "_", grid[i, "measure"], ".png")
    out2[out2$dist == grid[i, "dist"] &
             out2$measure == grid[i, "measure"] &
             out2$bias == grid[i, "bias"] &
             #out2$interval == grid[i, "interval"] &
             out2$large == grid[i, "large"] &
             out2$heterogeneity == grid[i, "heterogeneity"], ] %>%
        {
            if(grid[i, "measure"] == "coverage_prediction"){
                .[] %>% 
                    filter(grepl("Harmonic Mean|PI", method))
            } else {
                .[] %>% 
                    filter(!grepl("PI", method))
            }
        } %>% 
        plotPanels(data=., measure = grid[i, "measure"], by="method") +
        ggtitle(paste0("dist: ", grid[i, "dist"],
                       ", bias: ", grid[i, "bias"],
                       ", theta: 0.2, no. of large studies: ", grid[i, "large"],
                       #"interval: ", grid[i, "bias"],
                       ", simulation model: ", grid[i, "heterogeneity"], sep="")) +
        theme(plot.title = element_text(size = 10))
    ggsave(filename = single_plots[i],
           width = 12,
           height = 12,
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
                                  "figs/meanplots/",
                                  di,
                                  "_large_", la,
                                  "_bias_none",
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/meanplots/",
                                  di,
                                  "_large_", la,
                                  "_bias_moderate",
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/meanplots/",
                                  di,
                                  "_large_", la,
                                  "_bias_strong",
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " +append ",
                                  "figs/meanplots/BIAS",
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
                                  "figs/meanplots/",
                                  di,
                                  "_large_0",
                                  "_bias_", bi,
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/meanplots/",
                                  di,
                                  "_large_1",
                                  "_bias_", bi, 
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/meanplots/",
                                  di,
                                  "_large_2",
                                  "_bias_", bi,
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  " +append ",
                                  "figs/meanplots/TRIAL_SIZE",
                                  "_", di,
                                  "_bias_", bi,
                                  "_sim-mod_", he,
                                  "_", me, ".png"))
                }

## 3. heterogeneity model used for simulation: additive, multiplicative
for(di in dist)
    for(bi in bias)
        for(me in measure)
            for(la in large)
                if(me != "coverage_effects"){
                    system(paste0("convert ",
                                  "figs/meanplots/",
                                  di,
                                  "_large_", la,
                                  "_bias_", bi,
                                  "_sim-mod_additive",
                                  "_", me, ".png",
                                  " ",
                                  "figs/meanplots/",
                                  di,
                                  "_large_", la,
                                  "_bias_", bi, 
                                  "_sim-mod_multiplicative",
                                  "_", me, ".png",
                                  " ",
                                  " +append ",
                                  "figs/meanplots/SIM-MOD",
                                  "_", di,
                                  "_large_", la, 
                                  "_bias_", bi,
                                  "_", me, ".png"))
                }

## 4. Distribution: t, Gaussian
for(bi in bias)
    for(me in measure)
        for(la in large)
            for(he in heterogeneity)
                if(!(he == "multiplicative" & me == "coverage_effects")){
                    system(paste0("convert ",
                                  "figs/meanplots/Gaussian",
                                  "_large_", la,
                                  "_bias_", bi,
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  "figs/meanplots/t",
                                  "_large_", la,
                                  "_bias_", bi, 
                                  "_sim-mod_", he,
                                  "_", me, ".png",
                                  " ",
                                  " +append ",
                                  "figs/meanplots/DIST",
                                  "_large_", la, 
                                  "_bias_", bi,
                                  "_sim_mod_", he,
                                  "_", me, ".png"))
                }

unlink(single_plots)

## gamma_min summary statistics -----------------------------------------------------------

dir.create("figs/min_pH", showWarnings = FALSE)

out2 <- out %>% 
    bind_rows() %>%
    filter(grepl("gammaMin", measure), 
           grepl("Harmonic Mean|Harmonic Mean Additive|Harmonic Mean Multiplicative", method)) %>% 
    distinct() %>% 
    mutate(measure = gsub("^gammaMin_", "", measure),
           measure = case_when(measure == "min" ~ "Minimum",
                               measure == "firstQuart" ~ "1. Quartile",
                               measure == "mean" ~ "Mean",
                               measure == "median" ~ "Median",
                               measure == "thirdQuart" ~ "3. Quartile",
                               measure == "max" ~ "Maximum"),
           measure = ordered(measure, levels = c("Maximum", "3. Quartile", "Median", "Mean", "1. Quartile", "Minimum")))

out2 %>% select(dist) %>% unique() %>% pull() -> dist 
out2 %>% select(bias) %>% unique() %>% pull() -> bias 
out2 %>% select(large) %>% unique() %>% pull() -> large 
out2 %>% select(heterogeneity) %>% unique() %>% pull() -> heterogeneity
out2 %>% select(I2) %>% unique() %>% pull() -> I2 
#out2 %>% select(measure) %>% unique() %>% pull() -> me 

for(di in dist){
    for(la in large){
        for(bi in bias){
            for(he in heterogeneity){
                for(i in I2){
                    out2 %>% 
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

dir.create("figs/rel_freq", showWarnings = FALSE)

out2 <- out %>% 
    bind_rows() %>%
    filter(!grepl("gammaMin", measure),
           !grepl("mean", measure),
           method %in% c("Harmonic Mean", 
                         "Harmonic Mean Additive", 
                         "Harmonic Mean Multiplicative")) %>% 
    distinct() %>% 
    mutate(measure = gsub("n_", "", measure),
           measure = gsub("gt", "> ", measure),
           measure = ordered(measure, levels = rev(c("> 9", as.character(9:1)))),
           value = value/10000)

out2 %>% select(dist) %>% unique() %>% pull() -> dist 
out2 %>% select(bias) %>% unique() %>% pull() -> bias 
out2 %>% select(large) %>% unique() %>% pull() -> large 
out2 %>% select(heterogeneity) %>% unique() %>% pull() -> heterogeneity
out2 %>% select(I2) %>% unique() %>% pull() -> I2 
out2 %>% select(k) %>% unique() %>% pull() -> k

for(di in dist){
    for(la in large){
        for(bi in bias){
            for(he in heterogeneity){
                for(i in I2){
                    for(K in k){
                        out2 %>% 
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


