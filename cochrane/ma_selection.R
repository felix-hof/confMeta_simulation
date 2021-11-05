
rm(list = ls())
library(tidyverse)
library(ReplicationSuccess)


load("data/PubBias.RData") 
## 'data.ma': data of the meta analysis
## 'data':    results from studies


## subset to desired meta anaylses
data.ma %>% filter(estimable == "YES",
                   outcome.group == "efficacy",
           #        I2 >= 0,
                   nr.studies <= 10,
                   nr.results <= 10,
                   one.sign,
                   !hasDuplicatedEffects,
                   !duplicated(effect.key),
                   total1 > 0,
                   total2 > 0,
                   outcome.measure.merged == "MD",
                   outcome.flag == "CONT",
                   hasSubgroups == FALSE,
                   isSubgroup == FALSE,
                   unpublishedStudies == 0) -> dm


## get 'comparison.nr' and 'outcome.nr' from 'outsub.id'
dm %<>% 
    separate(outsub.id, c("to_remove", "comparison.nr", "outcome.nr"), remove = FALSE) %>%
    select(-to_remove) %>%
    mutate(comparison.nr = as.numeric(comparison.nr),
           outcome.nr = as.numeric(outcome.nr),
           id2 = paste(id, comparison.nr, outcome.nr))

## extract study estimates
data %>% filter(id %in% dm$id) %>%
    mutate(id2 = paste(id, comparison.nr, outcome.nr)) %>%
    right_join(dm, by = "id2") %>%
    select(id2, effect.se, effect.es, I2) %>%
    group_by(id2) %>%
    summarize(id2 = id2, I2 = I2, study = 1:n(),
              effect.es = effect.es, effect.se = effect.se) %>%
    ungroup() %>%
    filter(!is.na(effect.es)) ->
    study_effects
study_effects %>% summary()

## get number of disjunct hMeanChiSqCIs for each meta alaysis
study_effects %>% group_by(id2) %>%
    summarize(id2 = id2[1],
              nr_interval = nrow(hMeanChiSqCI(thetahat = effect.es, se = effect.se, alternative = "none")),
              I2 = I2[1]) %>%
    ungroup() ->
    nr_intervals
nr_intervals %>% summary()

nr_intervals %>% ggplot(mapping = aes(x = factor(nr_interval), y = I2)) +
    geom_boxplot()

nr_intervals %>% group_by(nr_interval) %>% summarize(n = n(), I2_median = median(I2))

## look at a study with 7 intervals and plot p-value function
study_effects %>% select(-I2) %>%
    right_join(nr_intervals, by = "id2") %>%
    filter(nr_interval == 7) %>%
    pull(id2) %>% unique()

study_effects %>% filter(id2 == "CD011308 1 3") %>% pull(effect.es) -> thetahat
study_effects %>% filter(id2 == "CD011308 1 3") %>% pull(effect.se) -> se

level <- 0.95; alpha <- 1 - level
muSeq <- seq(-6, 0, length.out = 1000)
pValueSeq <- hMeanChiSqMu(thetahat = thetahat, se = se,
                          alternative = "none", mu = muSeq)
(CIs <- hMeanChiSqCI(thetahat = thetahat, se = se, alternative = "none"))

plot(x = muSeq, y = pValueSeq, type = "l", panel.first = grid(lty = 1),
     xlab = expression(mu), ylab = "p-value")
abline(v = thetahat, h = alpha, lty = 2)
arrows(x0 = CIs[, 1], x1 = CIs[, 2], y0 = alpha,
       y1 = alpha, col = "darkgreen", lwd = 3, angle = 90, code = 3)

