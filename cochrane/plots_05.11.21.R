#### Plots for Cochrane examples of Harmonic Mean ####

### Content:
## Plots as drawn in powerpoint 
## Plots comparing multiplicative with additive methods

rm(list = ls())
library(tidyverse); theme_set(theme_bw())
source("../forestplot/forestplot.R")
source("../simulation/ReplicationSuccess_extension.R")
library(doParallel)
library(meta)
registerDoParallel(3)
load("data/PubData.RData")
#load("/Users/philipheesen/Library/Mobile Documents/com~apple~CloudDocs/Documents/Forschung/Biostatistik_UZH/data_data.ma.RData") 
#load("~/data_data.ma.RData") 

# filter MAs
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


## extract study estimates + tau2 
data %>% filter(id %in% dm$id) %>%
  mutate(id2 = paste(id, comparison.nr, outcome.nr)) %>%
  right_join(dm, by = "id2") %>%
  select(id2, effect.se, effect.es, I2, tau2, effect.x, ci_start, ci_end) %>%
  group_by(id2) %>%
  summarize(id2 = id2, I2 = I2, study = 1:n(),
            effect.es = effect.es, effect.se = effect.se, tau2=tau2, 
            effect=effect.x, ci_start=ci_start, ci_end=ci_end) %>%
  ungroup() %>%
  filter(!is.na(effect.es)) ->
  study_effects

# This throws an error for i = 2720
# get number of disjunct hMeanChiSqCIs for each meta alaysis & extract MA estimates
res <- foreach(i = seq_along(unique(study_effects$id2)), .combine = rbind) %dopar% {
  id2_i <- unique(study_effects$id2)[i]
  study_effects %>% filter(id2 == id2_i) -> subs
  subs %>% pull(effect.es) -> effect.es
  subs %>% pull(effect.se) -> effect.se
  mg <- metagen(TE = effect.es, seTE = effect.se, control = list(maxiter = 1000, stepadj = 0.5))
  hm <- hMeanChiSqCI(thetahat = effect.es, se = effect.se, alternative = "none", tau2 = mg$tau2)
  data.frame(id2 = id2_i,
             nr_interval = nrow(hm$CI),
             gamma_mean = hm$gammaMean,
             gamma_hmean = hm$gammaHMean,
             gamma_min = min(hm$gamma[,2]),
             I2 = mg$I2,
             Q = mg$Q,
             pval.Q = mg$pval.Q,
             tau2 = mg$tau2,
             nr_studies = length(mg$TE),
             p.random = weights(mg, comb.random=mg$comb.random)[4], # weights of random effects
             p.random.1 = var(weights(mg, comb.random=mg$comb.random)[4])) # variance of weights of random effects
}

## Plot study with one interval (1st plot in powerpoint)
# Leflunomide for rheumatoid arthritis
id2_ <- "CD002047 9 5"
study_effects %>% filter(id2 == id2_) %>% pull(effect.es) -> thetahat
study_effects %>% filter(id2 == id2_) %>% pull(effect.se) -> se
res %>% filter(id2 == id2_) %>% pull(I2) -> I2
res %>% filter(id2 == id2_) %>% pull(tau2) -> tau2
res %>%filter(id2 == id2_) %>% pull(gamma_min) -> gamma_min
res %>% filter(id2 == id2_) %>% pull(gamma_hmean) -> gamma_hmean
res %>% filter(id2 == id2_) %>% pull(gamma_mean) -> gamma_mean
res %>% filter(id2 == id2_) %>% pull(p.random.1) -> var_weights

level <- 0.95; alpha <- 1 - level
muSeq <- seq(min(thetahat-3*se), max(thetahat+3*se), length.out = 1000)
pValueSeq <- hMeanChiSqMu(thetahat = thetahat, se = se,
                          alternative = "none", mu = muSeq, tau2 = tau2[1]) # new implementation
(CIs <- hMeanChiSqCI(thetahat = thetahat, se = se, alternative = "none", tau2 = tau2[1]))

plot(x = muSeq, y = pValueSeq, type = "l", panel.first = grid(lty = 1),
     xlab = expression(mu), ylab = "p-value", col = "red", yaxt="n")+
  axis(side = 2, las=2)+
  abline(v = thetahat, h = alpha, lty = 2)+
  arrows(x0 = CIs$CI[, 1], x1 = CIs$CI[, 2], y0 = alpha,
         y1 = alpha, col = "red", lwd = 3, angle = 90, code = 3)+
  mtext(paste("I2:", round(I2, digits=2), "|", "Tau2:", round(tau2, digits=2), "|", 
              "gamma_Hmean", round(gamma_hmean, 2), 
              "|", "gamma_min:", round(gamma_min, digits=2), "|", "sd weights:", round(sqrt(var_weights), 2)), side=3)+
  mtext(paste("id2:", id2_), side=4)

print(data.frame(thetahat=thetahat, se=se) %>% arrange(thetahat))
print(data.frame(CIs$CI))
#forest(metagen(TE = thetahat, seTE = se)) # forest plot function from meta package

studies <- tibble(y = c(0.38, -0.22, -0.09),
                  lower = c(-0.06, -0.70, -0.30),
                  upper = c(0.82, 0.26, 0.11),
                  names = paste("study", seq_along(y)))

hMean <- tibble(lower = -0.45,
                upper = 0.59)

randomEffect <- tibble(lower = -0.30,
                       upper = 0.31)

forest(studies, hMean=hMean, randomEffect = randomEffect) # Florian's forest plot function


## Plot study with multiple intervals (2nd plot in powerpoint)
# Desmopressin for male nocturia
id2_ <- "CD012059 1 16"
study_effects %>% filter(id2 == id2_) %>% pull(effect.es) -> thetahat
study_effects %>% filter(id2 == id2_) %>% pull(effect.se) -> se
res %>% filter(id2 == id2_) %>% pull(I2) -> I2
res %>% filter(id2 == id2_) %>% pull(tau2) -> tau2
res %>%filter(id2 == id2_) %>% pull(gamma_min) -> gamma_min
res %>% filter(id2 == id2_) %>% pull(gamma_hmean) -> gamma_hmean
res %>% filter(id2 == id2_) %>% pull(gamma_mean) -> gamma_mean
res %>% filter(id2 == id2_) %>% pull(p.random.1) -> var_weights

level <- 0.95; alpha <- 1 - level
muSeq <- seq(min(thetahat-3*se), max(thetahat+3*se), length.out = 1000)
pValueSeq <- hMeanChiSqMu(thetahat = thetahat, se = se,
                          alternative = "none", mu = muSeq, tau2 = tau2[1]) # new implementation
(CIs <- hMeanChiSqCI(thetahat = thetahat, se = se, alternative = "none", tau2 = tau2[1]))


plot(x = muSeq, y = pValueSeq, type = "l", panel.first = grid(lty = 1),
     xlab = expression(mu), ylab = "p-value", col = "red", yaxt="n")+
  axis(side = 2, las=2)+
  abline(v = thetahat, h = alpha, lty = 2)+
  arrows(x0 = CIs$CI[, 1], x1 = CIs$CI[, 2], y0 = alpha,
         y1 = alpha, col = "red", lwd = 3, angle = 90, code = 3)+
  mtext(paste("I2:", round(I2, digits=2), "|", "Tau2:", round(tau2, digits=2), "|", 
              "gamma_Hmean", round(gamma_hmean, 2), 
              "|", "gamma_min:", round(gamma_min, digits=2), "|", "sd weights:", round(sqrt(var_weights), 2)), side=3)+
  mtext(paste("id2:", id2_), side=4)

print(data.frame(thetahat=thetahat, se=se) %>% arrange(thetahat))
print(data.frame(CIs$CI))
#forest(metagen(TE = thetahat, seTE = se)) # forest plot function from meta package

studies <- tibble(y = c(-0.96, -0.21, -1.30, -0.45, -0.36, -1.93),
                  lower = c(-1.30, -0.36, -1.70, -0.76, -0.60, -2.97),
                  upper = c(-0.61, -0.07, -0.89, -0.15, -0.11, -0.90),
                  names = paste("study", seq_along(y)))

hMean <- tibble(lower = c(-2.16, -1.48, -0.70),
                upper = c(-1.68, -0.73, -0.04))

randomEffect <- tibble(lower = -1.10,
                       upper = -0.36)

forest(studies, hMean=hMean, randomEffect = randomEffect) # Florian's forest plot function

##
# Leflunomide for rheumatoid arthritis
id2_ <- "CD002047 9 5"
study_effects %>% filter(id2 == id2_) %>% pull(effect.es) -> thetahat
study_effects %>% filter(id2 == id2_) %>% pull(effect.se) -> se
res %>% filter(id2 == id2_) %>% pull(I2) -> I2
res %>% filter(id2 == id2_) %>% pull(tau2) -> tau2
res %>%filter(id2 == id2_) %>% pull(gamma_min) -> gamma_min
res %>% filter(id2 == id2_) %>% pull(gamma_hmean) -> gamma_hmean
res %>% filter(id2 == id2_) %>% pull(gamma_mean) -> gamma_mean
res %>% filter(id2 == id2_) %>% pull(p.random.1) -> var_weights

level <- 0.95; alpha <- 1 - level
muSeq <- seq(min(thetahat-3*se), max(thetahat+3*se), length.out = 1000)
pValueSeq <- hMeanChiSqMu(thetahat = thetahat, se = se,
                          alternative = "none", mu = muSeq, tau2 = tau2[1]) # new implementation

phi <- estimatePhi(thetahat = thetahat, se = se)
pValueSeq_m <- hMeanChiSqMuPhi(thetahat = thetahat, se = se,
                               alternative = "none", mu = muSeq, phi = phi)

(CIs <- hMeanChiSqCI(thetahat = thetahat, se = se, alternative = "none", tau2 = tau2[1]))
(CIs_m <- hMeanChiSqCIphi(thetahat = thetahat, se = se, alternative = "none"))


plot(x = muSeq, y = pValueSeq, type = "l", panel.first = grid(lty = 1),
     xlab = expression(mu), ylab = "p-value", col = "red", yaxt="n")
  lines(x = muSeq, y = pValueSeq_m, col="green")
  axis(side = 2, las=2)
  abline(v = thetahat, h = alpha, lty = 2)
  arrows(x0 = CIs$CI[, 1], x1 = CIs$CI[, 2], y0 = alpha,
         y1 = alpha, col = "red", lwd = 3, angle = 90, code = 3)
  arrows(x0 = CIs_m$CI[, 1], x1 = CIs_m$CI[, 2], y0 = alpha,
         y1 = alpha, col = "green", lwd = 3, angle = 90, code = 3)
  mtext(paste("I2:", round(I2, digits=2), "|", "Tau2:", round(tau2, digits=2), "|", 
              "gamma_Hmean", round(gamma_hmean, 2), 
              "|", "gamma_min:", round(gamma_min, digits=2), "|", "sd weights:", round(sqrt(var_weights), 2)), side=3)
  mtext(paste("id2:", id2_), side=4)
  legend(x = min(muSeq),y = 0.9, legend = c("additive", "multiplicative"),
         col=c("red", "green"), lty = 1, cex = 0.8,
         box.lty = 0)

df.pValueSeq <- data.frame(x=muSeq, y=pValueSeq, y_m=pValueSeq_m)

##
# Desmopressin for male nocturia

id2_ <- "CD012059 1 16"
study_effects %>% filter(id2 == id2_) %>% pull(effect.es) -> thetahat
study_effects %>% filter(id2 == id2_) %>% pull(effect.se) -> se
res %>% filter(id2 == id2_) %>% pull(I2) -> I2
res %>% filter(id2 == id2_) %>% pull(tau2) -> tau2
res %>%filter(id2 == id2_) %>% pull(gamma_min) -> gamma_min
res %>% filter(id2 == id2_) %>% pull(gamma_hmean) -> gamma_hmean
res %>% filter(id2 == id2_) %>% pull(gamma_mean) -> gamma_mean
res %>% filter(id2 == id2_) %>% pull(p.random.1) -> var_weights

level <- 0.95; alpha <- 1 - level
muSeq <- seq(min(thetahat-3*se), max(thetahat+3*se), length.out = 1000)

pValueSeq <- hMeanChiSqMu(thetahat = thetahat, se = se,
                          alternative = "none", mu = muSeq, tau2 = tau2[1]) # additive

phi <- estimatePhi(thetahat = thetahat, se = se)
pValueSeq_m <- hMeanChiSqMuPhi(thetahat = thetahat, se = se,
                               alternative = "none", mu = muSeq, phi = phi) # multiplicative

(CIs <- hMeanChiSqCI(thetahat = thetahat, se = se, alternative = "none", tau2 = tau2[1])) # additive
(CIs_m <- hMeanChiSqCIphi(thetahat = thetahat, se = se, alternative = "none")) # multiplicative

plot(x = muSeq, y = pValueSeq, type = "l", panel.first = grid(lty = 1),
     xlab = expression(mu), ylab = "p-value", col = "red", yaxt="n")
  lines(x = muSeq, y = pValueSeq_m, col="green")
  axis(side = 2, las=2)
  abline(v = thetahat, h = alpha, lty = 2)
  arrows(x0 = CIs$CI[, 1], x1 = CIs$CI[, 2], y0 = alpha,
         y1 = alpha, col = "red", lwd = 3, angle = 90, code = 3)
  arrows(x0 = CIs_m$CI[, 1], x1 = CIs_m$CI[, 2], y0 = alpha,
         y1 = alpha, col = "green", lwd = 3, angle = 90, code = 3)
  mtext(paste("I2:", round(I2, digits=2), "|", "Tau2:", round(tau2, digits=2), "|", 
              "gamma_Hmean", round(gamma_hmean, 2), 
              "|", "gamma_min:", round(gamma_min, digits=2), "|", "sd weights:", round(sqrt(var_weights), 2)), side=3)
  mtext(paste("id2:", id2_), side=4)
  legend(x = min(muSeq), y = 0.9, legend = c("additive", "multiplicative"),
         col=c("red", "green"), lty = 1, cex = 0.8,
         box.lty = 0)

df.pValueSeq <- data.frame(x=muSeq, y=pValueSeq, y_m=pValueSeq_m)

