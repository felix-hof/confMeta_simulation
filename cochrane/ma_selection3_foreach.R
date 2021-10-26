rm(list = ls())
library(tidyverse)
library(ReplicationSuccess, lib = "~/git/ReplicationSuccess/lib")
library(meta)
library(doParallel)
registerDoParallel(3)

load("data/PubData.RData") 
## 'data.ma': data of the meta analysis
## 'data':    results from studies


## subset to desired meta anaylses
data.ma %>% filter(estimable == "YES",
                   outcome.group == "efficacy",
           #        I2 >= 0,
                   nr.studies <= 10,
                   nr.results <= 10,
                   nr.studies >= 2,
                   nr.results >= 2,
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
study_effects %>%
    distinct(id2, effect.es, .keep_all=TRUE) %>%
    group_by(id2) %>%
    summarize(id2 = id2[1],
              ## nr_interval = nrow(hMeanChiSqCI(thetahat = effect.es,
              ##                                 se = effect.se, alternative = "none")$CI),
              ## gamma_mean = hMeanChiSqCI(thetahat = effect.es,
              ##                           se = effect.se, alternative = "none")$gammaMean,
              ## gamma_hmean = hMeanChiSqCI(thetahat = effect.es,
              ##                            se = effect.se, alternative = "none")$gammaHMean,
              ## I2 = metagen(TE = effect.es, seTE = effect.se)$I2,
              ## Q = metagen(TE = effect.es, seTE = effect.se)$Q,
              ## pval.Q = metagen(TE = effect.es, seTE = effect.se)$pval.Q,
              tau2 = metagen(TE = effect.es, seTE = effect.se)$tau2) %>%
    ungroup() -> res


study_effects %>%
    distinct(id2, effect.es, .keep_all=TRUE) -> study_effects2

out <- foreach(i = seq_along(unique(study_effects2$id2)), .combine = rbind) %dopar% {
    id2_i <- unique(study_effects2$id2)[i]
    study_effects2 %>% filter(id2 == id2_i) -> subs
    subs %>% pull(effect.es) -> effect.es 
    subs %>% pull(effect.se) -> effect.se 
    hm <- hMeanChiSqCI(thetahat = effect.es, se = effect.se, alternative = "none")
    mg <- metagen(TE = effect.es, seTE = effect.se)
    data.frame(id2 = id2_i,
               nr_interval = nrow(hm$CI),
               gamma_mean = hm$gammaMean,
               gamma_hmean = hm$gammaHMean,
               I2 = mg$I2,
               Q = mg$Q,
               pval.Q = mg$pval.Q)
}
 

study_effects %>%
    distinct(id2, effect.es, .keep_all=TRUE) %>%
    filter(id2 == "CD000032 1 10") %>%
    pull(effect.es) -> effect.es


study_effects %>%
    distinct(id2, effect.es, .keep_all=TRUE) %>%
    filter(id2 == "CD000032 1 10") %>%
    pull(effect.se) -> effect.se


metagen(effect.es, effect.se)$tau2


res %>% summary()

res %>% group_by(nr_interval) %>% summarize(n = n(),
                                            I2_median = median(I2),
                                            gamma_mean = mean(gamma_mean),
                                            gamma_hmean = mean(gamma_hmean))

res %>% ggplot(mapping = aes(x = factor(nr_interval), y = I2)) +
    geom_boxplot()

res %>% ggplot(mapping = aes(x = gamma_mean, y = I2)) +
    geom_point() +
    geom_smooth()

res %>% ggplot(mapping = aes(x = gamma_hmean, y = I2)) +
    geom_point() +
    geom_smooth()

cor(res %>% pull(I2), res %>% pull(gamma_hmean), method = "spearman")
cor(res %>% pull(I2), res %>% pull(gamma_mean), method = "spearman")



## look at a study with 7 intervals and plot p-value function
study_effects %>% select(-I2) %>%
    right_join(res, by = "id2") %>%
    filter(nr_interval == 7) %>%
    pull(id2) %>% unique()

id2_ <- "CD012318 1 24"
study_effects %>% filter(id2 == id2_) %>% pull(effect.es) -> thetahat
study_effects %>% filter(id2 == id2_) %>% pull(effect.se) -> se

level <- 0.95; alpha <- 1 - level
muSeq <- seq(-6, 0, length.out = 1000)
pValueSeq <- hMeanChiSqMu(thetahat = thetahat, se = se,
                          alternative = "none", mu = muSeq)
(CIs <- hMeanChiSqCI(thetahat = thetahat, se = se, alternative = "none"))

plot(x = muSeq, y = pValueSeq, type = "l", panel.first = grid(lty = 1),
     xlab = expression(mu), ylab = "p-value")
abline(v = thetahat, h = alpha, lty = 2)
arrows(x0 = CIs$CI[, 1], x1 = CIs$CI[, 2], y0 = alpha,
       y1 = alpha, col = "darkgreen", lwd = 3, angle = 90, code = 3)

