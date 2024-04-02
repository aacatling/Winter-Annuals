# Data and script needed to reproduce results for JoE paper titled:
# Individual vital rates respond differently to local-scale environmental variation and neighbour removal
# 2024
# Alexandra Catling

#### Loading packages and data ####
#Data imported from data preparation sheet
source("data_preparation.R")
#vitaldata has all datasets combined, including:
# germination, survival, seed production, neighbour info, quadratic factors
# one row per subplot with seeds sown - NAs are very important since, e.g., survival info is only on germinated subplots

#Functions file
source("functions.R")

#Packages
library(ggplot2)
library(ggrepel)
library(MuMIn)
library(DHARMa)
library(glmmTMB)
library(kableExtra)
library(grid)
library(car)
library(sjPlot)
library(gridExtra)
library(corrplot)
library(scales)
library(cowplot)

#### Table S1-S4:Germination models ####

## ARCA new models
#N signif
arcagermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, arcadata)
arcagermmod1dharma <- simulateResiduals(arcagermfinalmod)
plot(arcagermmod1dharma)
vif(arcagermfinalmod)
summary(arcagermfinalmod)
testDispersion(arcagermfinalmod)

## hygl new models
#Cover signif
hyglgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, hygldata)
hyglgermmod1dharma <- simulateResiduals(hyglgermfinalmod)
plot(hyglgermmod1dharma)
vif(hyglgermfinalmod)
summary(hyglgermfinalmod)

## laro new models
#Cover signif
larogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, larodata)
larogermmod1dharma <- simulateResiduals(larogermfinalmod)
plot(larogermmod1dharma)
vif(larogermfinalmod)
summary(larogermfinalmod)

## peai new models
#Cover signif
peaigermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, peaidata)
peaigermmod1dharma <- simulateResiduals(peaigermfinalmod)
plot(peaigermmod1dharma)
vif(peaigermfinalmod)
summary(peaigermfinalmod)

## plde new models
#pH signif
pldegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, pldedata)
pldegermmod1dharma <- simulateResiduals(pldegermfinalmod)
plot(pldegermmod1dharma)
vif(pldegermfinalmod)
summary(pldegermfinalmod)

## pole new models
#N signif
polegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, poledata)
polegermmod1dharma <- simulateResiduals(polegermfinalmod)
plot(polegermmod1dharma)
vif(polegermfinalmod)
summary(polegermfinalmod)

## trcy new models
#Cover signif
trcygermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, trcydata)
trcygermmod1dharma <- simulateResiduals(trcygermfinalmod)
plot(trcygermmod1dharma)
vif(trcygermfinalmod)
summary(trcygermfinalmod)

## tror new models
#Cover signif
trorgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, trordata)
trorgermmod1dharma <- simulateResiduals(trorgermfinalmod)
plot(trorgermmod1dharma)
vif(trorgermfinalmod)
summary(trorgermfinalmod)

## vero new models
#Cover signif
verogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, verodata)
verogermmod1dharma <- simulateResiduals(verogermfinalmod)
plot(verogermmod1dharma)
vif(verogermfinalmod)
summary(verogermfinalmod)

#### Table S1-S4:Survival models ####
#arca
#N matters for survival, higher N higher survival. pH did not
#Simplifying -- watering DOES interact with cover
# water NO interaction with neighbour abundance
# cover DOES interact with neighbour abundance
#fundamentally the same as the original model

arcasurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + Cover:std_logp1_totalabund + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), arcadata)
arcasurvfinalmoddharma <- simulateResiduals(arcasurvfinalmod)
plot(arcasurvfinalmoddharma)
#good
summary(arcasurvfinalmod)

## hygl
hyglsurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, hygldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
hyglsurvmod1dharma <- simulateResiduals(hyglsurvmod1)
plot(hyglsurvmod1dharma)
#good
vif(hyglsurvmod1)
summary(hyglsurvmod1)

ggplot(hygldata, aes(y=surv_to_produce_seeds, x=Treatment, colour = Cover))+
  geom_boxplot()+
  geom_jitter(alpha = 0.3, width=0.1, height=0.2)+
  theme_classic()
#all dry plants in the shade died. HUGE estimates and SE from Dry:Shade

#N matters for survival, higher N higher survival. pH did not
#Simplyfing -- watering DOES interact with cover
# watering does interact with neighbour abundance (different from orig. model, not signif in final mod)
# cover no interaction with neighbour abundance

hyglsurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + Treatment:std_logp1_totalabund + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), hygldata)
hyglsurvfinalmoddharma <- simulateResiduals(hyglsurvfinalmod)
plot(hyglsurvfinalmoddharma)
vif(hyglsurvfinalmod)
#good
summary(hyglsurvfinalmod)
r.squaredGLMM(hyglsurvfinalmod)

## laro
larosurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, larodata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
larosurvmod1dharma <- simulateResiduals(larosurvmod1)
plot(larosurvmod1dharma)
#good
vif(larosurvmod1)
summary(larosurvmod1)

#Nothing signif
#fundamentally the same as the original model

larosurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), larodata)
larosurvfinalmoddharma <- simulateResiduals(larosurvfinalmod)
plot(larosurvfinalmoddharma)
#good
summary(larosurvfinalmod)

## peai
peaisurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, peaidata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
peaisurvmod1dharma <- simulateResiduals(peaisurvmod1)
plot(peaisurvmod1dharma)
#good
vif(peaisurvmod1)
summary(peaisurvmod1)

#Nothing signif - different from orig model, PC1:NA

peaisurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), peaidata)
peaisurvfinalmoddharma <- simulateResiduals(peaisurvfinalmod)
plot(peaisurvfinalmoddharma)
#good
summary(peaisurvfinalmod)

## plde
pldesurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, pldedata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
pldesurvmod1dharma <- simulateResiduals(pldesurvmod1)
plot(pldesurvmod1dharma)
#good
vif(pldesurvmod1)
#VIF OF 8 for pH - look into this later
summary(pldesurvmod1)

#Nothing signif - different from orig model, PC1:NA

pldesurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), pldedata)
pldesurvfinalmoddharma <- simulateResiduals(pldesurvfinalmod)
plot(pldesurvfinalmoddharma)
#good
summary(pldesurvfinalmod)

## trcy
trcysurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, trcydata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
trcysurvmod1dharma <- simulateResiduals(trcysurvmod1)
plot(trcysurvmod1dharma)
#good
vif(trcysurvmod1)
summary(trcysurvmod1)

#Simplifying -- watering DOES interact with cover
# water NO interaction with neighbour abundance
# cover DOES interact with neighbour abundance
#different from the original model

trcysurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + Cover:std_logp1_totalabund + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trcydata)
trcysurvfinalmoddharma <- simulateResiduals(trcysurvfinalmod)
plot(trcysurvfinalmoddharma)
#good
summary(trcysurvfinalmod)

## tror
trorsurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, trordata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
trorsurvmod1dharma <- simulateResiduals(trorsurvmod1)
plot(trorsurvmod1dharma)
#good
vif(trorsurvmod1)
summary(trorsurvmod1)

#Simplifying -- watering DOES interact with cover
# water NO interaction with neighbour abundance
# cover NO interaction with neighbour abundance
#fundamentally the same as the original model

trorsurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trordata)
trorsurvfinalmoddharma <- simulateResiduals(trorsurvfinalmod)
plot(trorsurvfinalmoddharma)
#good
summary(trorsurvfinalmod)

## vero
verosurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, verodata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
verosurvmod1dharma <- simulateResiduals(verosurvmod1)
plot(verosurvmod1dharma)
#good
vif(verosurvmod1)
summary(verosurvmod1)

#Nothing signif
#fundamentally the same as the original model

verosurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), verodata)
verosurvfinalmoddharma <- simulateResiduals(verosurvfinalmod)
plot(verosurvfinalmoddharma)
#good
summary(verosurvfinalmod)

##### Table S1-S4:Seed production models ####

#ARCA
#simplifying - signif Treatment:NA
#same as original
arcaseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                              Treatment:std_logp1_totalabund + (1|Site/Plot), 
                            family = nbinom2, seedarca)
arcaseedfinalmoddharma <- simulateResiduals(arcaseedfinalmod)
plot(arcaseedfinalmoddharma)
#not good
summary(arcaseedfinalmod)

#hygl
#won't converge with this new model
hyglseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedhygl)
hyglseedmod1dharma <- simulateResiduals(hyglseedmod1)
plot(hyglseedmod1dharma)
summary(hyglseedmod1)
#will converge using glmer.nb but drops TreatmentDry:CoverSun
#BECAUSE there is a complete separation issue - no dry plants in the shade produced viable seeds
#because there were no dry plants in the shade - they all died
#same estimate values between models excpt for TreatmentDry estimate
#glmmTMB optimiser stops and we don't get probabilities
#glmer.nb just drops that row from the output.
#Cover:Treatment not significant anyway so dropping it from final model
hyglseedmod2 <- glmer.nb(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                           Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), seedhygl)
hyglseedmod2dharma <- simulateResiduals(hyglseedmod2)
plot(hyglseedmod2dharma)
summary(hyglseedmod2)

# ggplot(seedhygl, aes(x=Cover, y= No_viable_seeds_grouped, colour = Treatment)) +
#   geom_jitter(alpha=0.5)+
#   theme_classic()

#simplifying - nothing signif
hyglseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedhygl)
hyglseedfinalmoddharma <- simulateResiduals(hyglseedfinalmod)
plot(hyglseedfinalmoddharma)
#good
summary(hyglseedfinalmod)
#same as original

#peai
peaiseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedpeai)
peaiseedmod1dharma <- simulateResiduals(peaiseedmod1)
plot(peaiseedmod1dharma)
summary(peaiseedmod1)

#simplifying - nothing signif, different from original
peaiseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedpeai)
peaiseedfinalmoddharma <- simulateResiduals(peaiseedfinalmod)
plot(peaiseedfinalmoddharma)
#good
summary(peaiseedfinalmod)

#laro
laroseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedlaro)
laroseedmod1dharma <- simulateResiduals(laroseedmod1)
plot(laroseedmod1dharma)
summary(laroseedmod1)

#simplifying - nothing signif (Cover:NA 0.054!), different from original
laroseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedlaro)
laroseedfinalmoddharma <- simulateResiduals(laroseedfinalmod)
plot(laroseedfinalmoddharma)
#good
summary(laroseedfinalmod)

#plde
pldeseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedplde)
pldeseedmod1dharma <- simulateResiduals(pldeseedmod1)
plot(pldeseedmod1dharma)
summary(pldeseedmod1)

#simplifying - nothing signif, same as original
pldeseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedplde)
pldeseedfinalmoddharma <- simulateResiduals(pldeseedfinalmod)
plot(pldeseedfinalmoddharma)
#good
summary(pldeseedfinalmod)

#trcy
trcyseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedtrcy)
trcyseedmod1dharma <- simulateResiduals(trcyseedmod1)
plot(trcyseedmod1dharma)
summary(trcyseedmod1)

#simplifying - signif Cover:NA (and Treatment:Cover 0.052!), differnt from original
trcyseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                              Cover:std_logp1_totalabund + (1|Site/Plot), 
                            family = nbinom2, seedtrcy)
trcyseedfinalmoddharma <- simulateResiduals(trcyseedfinalmod)
plot(trcyseedfinalmoddharma)
#good
summary(trcyseedfinalmod)

#tror
trorseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedtror)
trorseedmod1dharma <- simulateResiduals(trorseedmod1)
plot(trorseedmod1dharma)
summary(trorseedmod1)

#simplifying - signif Treatment:Cover
#same as original
trorseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                              Treatment:Cover + (1|Site/Plot), 
                            family = nbinom2, seedtror)
trorseedfinalmoddharma <- simulateResiduals(trorseedfinalmod)
plot(trorseedfinalmoddharma)
#good
summary(trorseedfinalmod)

#vero
veroseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedvero)
veroseedmod1dharma <- simulateResiduals(veroseedmod1)
plot(veroseedmod1dharma)
summary(veroseedmod1)

#simplifying - nothing signif, same as original
veroseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedvero)
veroseedfinalmoddharma <- simulateResiduals(veroseedfinalmod)
plot(veroseedfinalmoddharma)
#good
summary(veroseedfinalmod)

##### Table S1-S4:Population growth models ####

## arca
arcalambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdaarca)
arcalambdamod1dharma <- simulateResiduals(arcalambdamod1)
plot(arcalambdamod1dharma)
#not good
vif(arcalambdamod1)
summary(arcalambdamod1)
#Model simplification step - No interactions significant, removing all
#same as original
arcalambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + (1|Site/Plot), lambdaarca)
arcalambdafinalmoddharma <- simulateResiduals(arcalambdafinalmod)
plot(arcalambdafinalmoddharma)
#okay
summary(arcalambdafinalmod)

## HYGL
hygllambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdahygl)
hygllambdamod1dharma <- simulateResiduals(hygllambdamod1)
plot(hygllambdamod1dharma)
#okay
vif(hygllambdamod1)
summary(hygllambdamod1)
#Model simplification step - No interactions significant, removing all
#same as original
hygllambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + (1|Site/Plot), lambdahygl)
hygllambdafinalmoddharma <- simulateResiduals(hygllambdafinalmod)
plot(hygllambdafinalmoddharma)
#okay
summary(hygllambdafinalmod)

## laro
larolambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdalaro)
larolambdamod1dharma <- simulateResiduals(larolambdamod1)
plot(larolambdamod1dharma)
#okay
vif(larolambdamod1)
summary(larolambdamod1)
#Model simplification step - signif Cover:Neighbours, same as original
larolambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdalaro)
larolambdafinalmoddharma <- simulateResiduals(larolambdafinalmod)
plot(larolambdafinalmoddharma)
#okay
summary(larolambdafinalmod)

## peai
peailambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdapeai)
peailambdamod1dharma <- simulateResiduals(peailambdamod1)
plot(peailambdamod1dharma)
#okay
vif(peailambdamod1)
summary(peailambdamod1)
#Model simplification step - signif Treatment:Cover, Cover:Neighbours
#same as original
peailambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + Treatment:Cover + Cover:Neighbours01 + (1|Site/Plot), lambdapeai)
peailambdafinalmoddharma <- simulateResiduals(peailambdafinalmod)
plot(peailambdafinalmoddharma)
#okay
summary(peailambdafinalmod)

## plde
pldelambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdaplde)
pldelambdamod1dharma <- simulateResiduals(pldelambdamod1)
plot(pldelambdamod1dharma)
#okay
vif(pldelambdamod1)
summary(pldelambdamod1)
#Model simplification step - No interactions significant, removing all
#same as original
pldelambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + (1|Site/Plot), lambdaplde)
pldelambdafinalmoddharma <- simulateResiduals(pldelambdafinalmod)
plot(pldelambdafinalmoddharma)
#okay
summary(pldelambdafinalmod)

## trcy
trcylambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdatrcy)
trcylambdamod1dharma <- simulateResiduals(trcylambdamod1)
plot(trcylambdamod1dharma)
#okay
vif(trcylambdamod1)
summary(trcylambdamod1)
#Model simplification step - signif Treatment:Neighbours and Cover:Neighbours
#same as original
trcylambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdatrcy)
trcylambdafinalmoddharma <- simulateResiduals(trcylambdafinalmod)
plot(trcylambdafinalmoddharma)
#okay
summary(trcylambdafinalmod)

## tror
trorlambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdatror)
trorlambdamod1dharma <- simulateResiduals(trorlambdamod1)
plot(trorlambdamod1dharma)
#okay
vif(trorlambdamod1)
summary(trorlambdamod1)
#Model simplification step - No interactions significant, removing all
#same as original
trorlambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + (1|Site/Plot), lambdatror)
trorlambdafinalmoddharma <- simulateResiduals(trorlambdafinalmod)
plot(trorlambdafinalmoddharma)
#okay
summary(trorlambdafinalmod)

## vero
verolambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdavero)
verolambdamod1dharma <- simulateResiduals(verolambdamod1)
plot(verolambdamod1dharma)
#okay
vif(verolambdamod1)
summary(verolambdamod1)
#Model simplification step - signif Cover:Neighbours 
#different from original
verolambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdavero)
verolambdafinalmoddharma <- simulateResiduals(verolambdafinalmod)
plot(verolambdafinalmoddharma)
#okay
summary(verolambdafinalmod)


#### Table S1-S4: Extracting marginal and conditional r squared values ####
#Extracting values for theoretical marginal R squared
germ_model_list <- list(arcagermfinalmod, hyglgermfinalmod, larogermfinalmod, peaigermfinalmod, pldegermfinalmod, polegermfinalmod, trcygermfinalmod, trorgermfinalmod, verogermfinalmod)
rsquared = lapply(1:length(germ_model_list), function(x) {
  as.data.frame(r.squaredGLMM(germ_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
germ_rsquared_table <- do.call("rbind", rsquared)

#Rename species
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '6'] <- 'Podolepis lessonii')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '7'] <- 'Trachymene cyanopetala')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '8'] <- 'Trachymene ornata')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '9'] <- 'Goodenia rosea')

#### For Figure 3 - rates ~ cover ####
### arca ####
coef(arcagermfinalmod)
arca_germ_pred<-glmm.predict(mod=arcagermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
arca_germ_pred$Cover <- c('Shade', 'Sun')
a <- ggplot()+
  geom_jitter(data=arcadata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=arca_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(arcasurvfinalmod)
arca_surv_pred<-glmm.predict(mod=arcasurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
arca_surv_pred$Cover <- c('Shade', 'Sun')
b <- ggplot()+
  geom_jitter(data=arcadata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=arca_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

#plotting on a log scale, so need seeds+1 and +1 to negative lower limits
#exponentiating after logging, log_link=TRUE IS exponentiating it for us
coef(arcaseedfinalmod)
arca_seed_pred<-glmm.predict(mod=arcaseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*0, 0*0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
arca_seed_pred$Cover <- c('Shade', 'Sun')
#added one to the negative lower limit
arca_seed_pred$lower <- c(1.165081, 1.447482)
c <- ggplot()+
  geom_jitter(data=seedarca, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=arca_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

#in the absence of neighbours
coef(arcalambdafinalmod)
arca_lambda_pred<-glmm.predict(mod=arcalambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
arca_lambda_pred$Cover <- c('Shade', 'Sun')
d <- ggplot()+
  geom_jitter(data=lambdaarca, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=arca_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### hygl ####
coef(hyglgermfinalmod)
hygl_germ_pred<-glmm.predict(mod=hyglgermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
hygl_germ_pred$Cover <- c('Shade', 'Sun')
e <- ggplot()+
  geom_jitter(data=hygldata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=hygl_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1), n.breaks=3)+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hyglsurvfinalmod)
hygl_surv_pred<-glmm.predict(mod=hyglsurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), 0*0, 0*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
hygl_surv_pred$Cover <- c('Shade', 'Sun')
f <- ggplot()+
  geom_jitter(data=hygldata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=hygl_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hyglseedfinalmod)
hygl_seed_pred<-glmm.predict(mod=hyglseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
hygl_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
g <- ggplot()+
  geom_jitter(data=seedhygl, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=hygl_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hygllambdafinalmod)
hygl_lambda_pred<-glmm.predict(mod=hygllambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
hygl_lambda_pred$Cover <- c('Shade', 'Sun')
#negative values here are fine! Population growth rate
h <- ggplot()+
  geom_jitter(data=lambdahygl, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=hygl_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### laro ####
coef(larogermfinalmod)
laro_germ_pred<-glmm.predict(mod=larogermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
laro_germ_pred$Cover <- c('Shade', 'Sun')
i <- ggplot()+
  geom_jitter(data=larodata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=laro_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1), n.breaks=3)+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(larosurvfinalmod)
laro_surv_pred<-glmm.predict(mod=larosurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
laro_surv_pred$Cover <- c('Shade', 'Sun')
j <- ggplot()+
  geom_jitter(data=larodata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=laro_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(laroseedfinalmod)
laro_seed_pred<-glmm.predict(mod=laroseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
laro_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
k <- ggplot()+
  geom_jitter(data=seedlaro, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=laro_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(larolambdafinalmod)
laro_lambda_pred<-glmm.predict(mod=larolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
laro_lambda_pred$Cover <- c('Shade', 'Sun')
l <- ggplot()+
  geom_jitter(data=lambdalaro, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=laro_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### peai ####
coef(peaigermfinalmod)
peai_germ_pred<-glmm.predict(mod=peaigermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
peai_germ_pred$Cover <- c('Shade', 'Sun')
m <- ggplot()+
  geom_jitter(data=peaidata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=peai_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1), n.breaks=3)+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(peaisurvfinalmod)
peai_surv_pred<-glmm.predict(mod=peaisurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
peai_surv_pred$Cover <- c('Shade', 'Sun')
n <- ggplot()+
  geom_jitter(data=peaidata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=peai_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(peaiseedfinalmod)
peai_seed_pred<-glmm.predict(mod=peaiseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
peai_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
o <- ggplot()+
  geom_jitter(data=seedpeai, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=peai_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  annotate("text", x = 1.5, y= 100, label = "*", family = "", size = 20, colour = "red")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
coef(peailambdafinalmod)
peai_lambda_pred<-glmm.predict(mod=peailambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
peai_lambda_pred$Cover <- c('Shade', 'Sun')
p <- ggplot()+
  geom_jitter(data=lambdapeai, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=peai_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### plde ####
coef(pldegermfinalmod)
plde_germ_pred<-glmm.predict(mod=pldegermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plde_germ_pred$Cover <- c('Shade', 'Sun')
q <- ggplot()+
  geom_jitter(data=pldedata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=plde_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldesurvfinalmod)
plde_surv_pred<-glmm.predict(mod=pldesurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plde_surv_pred$Cover <- c('Shade', 'Sun')
r <- ggplot()+
  geom_jitter(data=pldedata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=plde_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldeseedfinalmod)
plde_seed_pred<-glmm.predict(mod=pldeseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plde_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
s <- ggplot()+
  geom_jitter(data=seedplde, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=plde_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  annotate("text", x = 1.5, y= 30, label = "*", family = "", size = 20, colour = "red")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldelambdafinalmod)
plde_lambda_pred<-glmm.predict(mod=pldelambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plde_lambda_pred$Cover <- c('Shade', 'Sun')
t <- ggplot()+
  geom_jitter(data=lambdaplde, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=plde_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### trcy ####
coef(trcygermfinalmod)
trcy_germ_pred<-glmm.predict(mod=trcygermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
trcy_germ_pred$Cover <- c('Shade', 'Sun')
u <- ggplot()+
  geom_jitter(data=trcydata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=trcy_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1), n.breaks=3)+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcysurvfinalmod)
trcy_surv_pred<-glmm.predict(mod=trcysurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
trcy_surv_pred$Cover <- c('Shade', 'Sun')
v <- ggplot()+
  geom_jitter(data=trcydata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=trcy_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcyseedfinalmod)
trcy_seed_pred<-glmm.predict(mod=trcyseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, c(1,0)*0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
trcy_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
w <- ggplot()+
  geom_jitter(data=seedtrcy, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=trcy_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcylambdafinalmod)
trcy_lambda_pred<-glmm.predict(mod=trcylambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0*0, 0*0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
trcy_lambda_pred$Cover <- c('Shade', 'Sun')
x <- ggplot()+
  geom_jitter(data=lambdatrcy, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=trcy_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### tror ####
coef(trorgermfinalmod)
tror_germ_pred<-glmm.predict(mod=trorgermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
tror_germ_pred$Cover <- c('Shade', 'Sun')
y <- ggplot()+
  geom_jitter(data=trordata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=tror_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1), n.breaks=3)+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorsurvfinalmod)
tror_surv_pred<-glmm.predict(mod=trorsurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0)),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
tror_surv_pred$Cover <- c('Shade', 'Sun')
z <- ggplot()+
  geom_jitter(data=trordata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=tror_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorseedfinalmod)
tror_seed_pred<-glmm.predict(mod=trorseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0)),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
tror_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
ab <- ggplot()+
  geom_jitter(data=seedtror, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=tror_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorlambdafinalmod)
tror_lambda_pred<-glmm.predict(mod=trorlambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
tror_lambda_pred$Cover <- c('Shade', 'Sun')
bc <- ggplot()+
  geom_jitter(data=lambdatror, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=tror_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(n.breaks=3)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### vero ####
coef(verogermfinalmod)
vero_germ_pred<-glmm.predict(mod=verogermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
vero_germ_pred$Cover <- c('Shade', 'Sun')
cd <- ggplot()+
  geom_jitter(data=verodata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=vero_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1), n.breaks=3)+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(verosurvfinalmod)
vero_surv_pred<-glmm.predict(mod=verosurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
vero_surv_pred$Cover <- c('Shade', 'Sun')
de <- ggplot()+
  geom_jitter(data=verodata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=vero_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), n.breaks=3)+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(veroseedfinalmod)
vero_seed_pred<-glmm.predict(mod=veroseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
vero_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
ef <- ggplot()+
  geom_jitter(data=seedvero, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=vero_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(verolambdafinalmod)
vero_lambda_pred<-glmm.predict(mod=verolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
vero_lambda_pred$Cover <- c('Shade', 'Sun')
fg <- ggplot()+
  geom_jitter(data=lambdavero, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=vero_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(n.breaks=3)+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())
#### Figure 3 - plotting cover_panel from above ####
arca <- ggplot()+geom_text(aes(label = "A. calendula"), size = 14, y = 0.5, x = 0.5)+theme_void()
hygl <- ggplot()+geom_text(aes(label = "H. glutinosum"), size = 14, y = 0.5, x = 0.5)+theme_void()
laro <- ggplot()+geom_text(aes(label = "L. rosea"), size = 14, y = 0.5, x = 0.5)+theme_void()
peai <- ggplot()+geom_text(aes(label = "P. airoides"), size = 14, y = 0.5, x = 0.5)+theme_void()
plde <- ggplot()+geom_text(aes(label = "P. debilis"), size = 14, y = 0.5, x = 0.5)+theme_void()
trcy <- ggplot()+geom_text(aes(label = "T. cyanopetala"), size = 14, y = 0.5, x = 0.5)+theme_void()
tror <- ggplot()+geom_text(aes(label = "T. ornata"), size = 14, y = 0.5, x = 0.5)+theme_void()
vero <- ggplot()+geom_text(aes(label = "G. rosea"), size = 14, y = 0.5, x = 0.5)+theme_void()

empty <- ggplot()+geom_text(aes(label = "", angle = 90), x = 0.5, y = 0.5)+theme_void()
germ <- ggplot()+geom_text(aes(label = "Probability of emergence", angle = 90), size = 14, y = 0.5, x = 0.5)+theme_void()
surv <- ggplot()+geom_text(aes(label = "Probability of survival", angle = 90), size = 14, y = 0.5, x = 0.5)+theme_void()
seed <- ggplot()+geom_text(aes(label = "Number of seeds produced", angle = 90), size = 14, y = 0.5, x = 0.5)+theme_void()
lambda <- ggplot()+geom_text(aes(label = "Population growth rate", angle = 90), size = 14, y = 0.5, x = 0.5)+theme_void()
species <- ggplot()+geom_text(aes(label = "Species"), size = 14, y = 0.5, x = 0.5)+theme_void()
emergence <- ggplot()+geom_text(aes(label = "Emergence"), size = 14, y = 0.5, x = 0.5)+theme_void()
survival <- ggplot()+geom_text(aes(label = "Survival"), size = 14, y = 0.5, x = 0.5)+theme_void()
seedproduction <- ggplot()+geom_text(aes(label = "Seed production"), size = 14, y = 0.5, x = 0.5)+theme_void()
population <- ggplot()+geom_text(aes(label = "Population growth rate"), size = 14, y = 0.5, x = 0.5)+theme_void()

dev.off()
pdf("Output/Figures/panel_cover_final_2024.pdf", width=21, height=21)
plot_grid(empty,empty,empty,empty,empty,empty,empty,empty,
          empty,a,empty,b,empty,c,empty,d,empty,e,empty,f,empty,g,empty,h,
          empty,i,empty,j,empty,k,empty,l,empty,m,empty,n,empty,o,empty,p,
          empty,q,empty,r,empty,s,empty,t,empty,u,empty,v,empty,w,empty,x,
          empty,y,empty,z,empty,ab,empty,bc,empty,cd,empty,de,empty,ef,empty,fg, 
          align="hv", ncol=8, rel_widths = c(3,4,1,4,1,4,1,4), rel_heights=c(1,4,4,4,4,4,4,4,4))
dev.off()
### manually adding axis and species labels in Adobe Illustrator

#### Figure 4 - Panel neighbours*cover new ####
pdf("Output/Figures/panel_NA_int_2024.pdf", width=21, height=21)
par(mfrow=c(8,3), oma = c(12, 20, 5, 1), mar =c(3,10,1,1))

####
#ARCA
#Survival - signif int
x_to_plot_sun <- seq.func(arcadata$std_logp1_totalabund[arcadata$Cover=='Sun'])
x_to_plot_shade <- seq.func(arcadata$std_logp1_totalabund[arcadata$Cover=='Shade'])
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha(ifelse(Cover=="Sun", "#CC79A7", "#0072B2"), 0.3), ylab="", xlab = NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, arcadata)
model <- arcasurvfinalmod
#Sun - blue
plotted.pred.sun <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0, x_to_plot_sun, 0, 0*0, 0*0, 0*x_to_plot_sun), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_sun, pred = plotted.pred.sun$y, upper = plotted.pred.sun$upper, lower = plotted.pred.sun$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Shade - red
plotted.pred.shade <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 1, 0, 0, x_to_plot_shade, 0, 0*1, 0*1, 1*x_to_plot_shade), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_shade, pred = plotted.pred.shade$y, upper = plotted.pred.shade$upper, lower = plotted.pred.shade$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#Fecundity
x_to_plot<-seq.func(seedarca$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 150), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedarca)
model <- arcaseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
arca_pred_nonbh<-glmm.predict(mod=arcalambdafinalmod, newdat=data.frame(1, 0, 0, 0, 0, 0, 0),
                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
arca_pred_nonbh$Neighbours01 <- 'Neighbours0'
arca_pred_nbh<-glmm.predict(mod=arcalambdafinalmod, newdat=data.frame(1, 0, 0, 0, 0, 0, 1),
                            se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
arca_pred_nbh$Neighbours01 <- 'Neighbours1'
arca_pred <- rbind(arca_pred_nonbh, arca_pred_nbh)
arca_pred$x <- c(0.5,1)
lambdaarca <- within(lambdaarca, x[Neighbours01== 'Neighbours0'] <- '0.5')
lambdaarca <- within(lambdaarca, x[Neighbours01== 'Neighbours1'] <- '1')
plot(log_lambda ~ x, pch=19, ylim=c(-0.8,3.1),xlim=c(0.25,1.25),col=alpha("grey60",0.3), ylab=NA, xlab=NA, xaxt="n", tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdaarca)
points(y ~ x, arca_pred, pch = 17, lwd=3,cex= 2.5)
arrows(x0=arca_pred$x, y0=arca_pred$lower, x1=arca_pred$x, y1=arca_pred$upper, code=3, angle=90, length=0.1, lwd=3)
text(x = 0.75, y = 1.7, "*", cex = 8, col = "red")

#Adding HYGL
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, hygldata)
model <- hyglsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0, 0*0, 0*0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedhygl$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedhygl)
model <- hyglseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
hygl_pred_nonbh<-glmm.predict(mod=hygllambdafinalmod, newdat=data.frame(1, 0, 0, 0, 0, 0, 0),
                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
hygl_pred_nonbh$Neighbours01 <- 'Neighbours0'
hygl_pred_nbh<-glmm.predict(mod=hygllambdafinalmod, newdat=data.frame(1, 0, 0, 0, 0, 0, 1),
                            se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
hygl_pred_nbh$Neighbours01 <- 'Neighbours1'
hygl_pred <- rbind(hygl_pred_nonbh, hygl_pred_nbh)
hygl_pred$x <- c(0.5,1)
lambdahygl$x <- '1'
lambdahygl <- within(lambdahygl, x[Neighbours01== 'Neighbours0'] <- '0.5')
lambdahygl <- within(lambdahygl, x[Neighbours01== 'Neighbours1'] <- '1')
plot(log_lambda ~ x, pch=19, ylim=c(-0.5,1.7),xlim=c(0.25,1.25),col=alpha("grey60",0.3), ylab=NA, xlab=NA, xaxt="n", tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdahygl)
points(y ~ x, hygl_pred, pch = 17, lwd=3,cex= 2.5)
arrows(x0=hygl_pred$x, y0=hygl_pred$lower, x1=hygl_pred$x, y1=hygl_pred$upper, code=3, angle=90, length=0.1, lwd=3)

#Adding LARO
x_to_plot<-seq.func(larodata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, larodata)
model <- larosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedlaro$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedlaro)
model <- laroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
laro_pred_nonbh<-glmm.predict(mod=larolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, c(1,0)*0),
                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
laro_pred_nonbh$Cover <- c('Shade', 'Sun')
laro_pred_nonbh$Neighbours01 <- 'Neighbours0'
laro_pred_nbh<-glmm.predict(mod=larolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 1, c(1,0)*1),
                            se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
laro_pred_nbh$Cover <- c('Shade', 'Sun')
laro_pred_nbh$Neighbours01 <- 'Neighbours1'
laro_pred <- rbind(laro_pred_nonbh,laro_pred_nbh)
laro_pred$x_label <- c("Shade:Nbh0", "Sun:Nbh0", "Shade:Nbh1", "Sun:Nbh1")
#for this figure I actually need absent sun and shade then present sun and shade
laro_pred$x <- c(1,0.5,2.5,2)
lambdalaro$x <- "1"
lambdalaro <- within(lambdalaro, x[Cover== 'Sun' & Neighbours01== 'Neighbours1'] <- '2')
lambdalaro <- within(lambdalaro, x[Cover== 'Sun' & Neighbours01== 'Neighbours0'] <- '0.5')
lambdalaro <- within(lambdalaro, x[Cover== 'Shade' & Neighbours01== 'Neighbours1'] <- '2.5')
lambdalaro <- within(lambdalaro, x[Cover== 'Shade' & Neighbours01== 'Neighbours0'] <- '1')
plot(log_lambda ~ x, pch=19, ylim=c(-0.7,2.1),xlim=c(0.25,2.75),col=alpha(ifelse(Cover=="Sun", "#CC79A7", "#0072B2"),0.2), ylab=NA, xlab=NA, xaxt="n", tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdalaro)
points(y ~ x, laro_pred, pch = 17, col=ifelse(Cover=="Sun", "#CC79A7", "#0072B2"), lwd=3,cex= 2.5)
arrows(x0=laro_pred$x, y0=laro_pred$lower, x1=laro_pred$x, y1=laro_pred$upper, col=c("#0072B2", "#CC79A7", "#0072B2", "#CC79A7"), code=3, angle=90, length=0.1, lwd=3)

#Adding PEAI
x_to_plot<-seq.func(peaidata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, peaidata)
model <- peaisurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedpeai$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedpeai)
model <- peaiseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
peai_pred_nonbh<-glmm.predict(mod=peailambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0,0*c(1,0),0*c(1,0),c(1,0)*0),
                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
peai_pred_nonbh$Cover <- c('Shade', 'Sun')
peai_pred_nonbh$Neighbours01 <- 'Neighbours0'
peai_pred_nbh<-glmm.predict(mod=peailambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 1,0*c(1,0),0*c(1,0),c(1,0)*1),
                            se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
peai_pred_nbh$Cover <- c('Shade', 'Sun')
peai_pred_nbh$Neighbours01 <- 'Neighbours1'
peai_pred <- rbind(peai_pred_nonbh,peai_pred_nbh)
peai_pred$x_label <- c("Shade:Nbh0", "Sun:Nbh0", "Shade:Nbh1", "Sun:Nbh1")
peai_pred$x <- c(1,0.5,2.5,2)
lambdapeai$x <- "1"
lambdapeai <- within(lambdapeai, x[Cover== 'Sun' & Neighbours01== 'Neighbours1'] <- '2')
lambdapeai <- within(lambdapeai, x[Cover== 'Sun' & Neighbours01== 'Neighbours0'] <- '0.5')
lambdapeai <- within(lambdapeai, x[Cover== 'Shade' & Neighbours01== 'Neighbours1'] <- '2.5')
lambdapeai <- within(lambdapeai, x[Cover== 'Shade' & Neighbours01== 'Neighbours0'] <- '1')
plot(log_lambda ~ x, pch=19, ylim=c(-1.1,4.6),xlim=c(0.25,2.75),col=alpha(ifelse(Cover=="Sun", "#CC79A7", "#0072B2"),0.2), ylab=NA, xlab=NA, xaxt="n", tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdapeai)
points(y ~ x, peai_pred, pch = 17, col=ifelse(Cover=="Sun", "#CC79A7", "#0072B2"), lwd=3,cex= 2.5)
arrows(x0=peai_pred$x, y0=peai_pred$lower, x1=peai_pred$x, y1=peai_pred$upper, col=c("#0072B2", "#CC79A7", "#0072B2", "#CC79A7"), code=3, angle=90, length=0.1, lwd=3)

#Adding PLDE
x_to_plot<-seq.func(pldedata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, pldedata)
model <- pldesurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedplde$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedplde)
model <- pldeseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
text(x = 2.3, y = 50, "*", cex = 8, col = "red")
#Lambda
plde_pred_nonbh<-glmm.predict(mod=pldelambdafinalmod, newdat=data.frame(1, 0, 0, 0, 0, 0, 0),
                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plde_pred_nonbh$Neighbours01 <- 'Neighbours0'
plde_pred_nbh<-glmm.predict(mod=pldelambdafinalmod, newdat=data.frame(1, 0, 0, 0, 0, 0, 1),
                            se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plde_pred_nbh$Neighbours01 <- 'Neighbours1'
plde_pred <- rbind(plde_pred_nonbh, plde_pred_nbh)
plde_pred$x <- c(0.5,1)
lambdaplde$x <- "1"
lambdaplde <- within(lambdaplde, x[Neighbours01== 'Neighbours0'] <- '0.5')
lambdaplde <- within(lambdaplde, x[Neighbours01== 'Neighbours1'] <- '1')
plot(log_lambda ~ x, pch=19, ylim=c(-1.2,2.5),xlim=c(0.25,1.25),col=alpha("grey60",0.3), ylab=NA, xlab=NA, xaxt="n", tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdaplde)
points(y ~ x, plde_pred, pch = 17, lwd=3,cex= 2.5)
arrows(x0=plde_pred$x, y0=plde_pred$lower, x1=plde_pred$x, y1=plde_pred$upper, code=3, angle=90, length=0.1, lwd=3)
text(x = 0.75, y = 1.8, "*", cex = 8, col = "red")

#Adding TRCY
#Survival
x_to_plot_sun <- seq.func(trcydata$std_logp1_totalabund[trcydata$Cover=='Sun'])
x_to_plot_shade <- seq.func(trcydata$std_logp1_totalabund[trcydata$Cover=='Shade'])
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha(ifelse(Cover=="Sun", "#CC79A7", "#0072B2"), 0.3), ylab="", xlab = NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, trcydata)
model <- trcysurvfinalmod
#Sun - blue
plotted.pred.sun <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_sun, pred = plotted.pred.sun$y, upper = plotted.pred.sun$upper, lower = plotted.pred.sun$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Shade - red
plotted.pred.shade <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 1, 0, 0, x_to_plot, 0, 0*1, 0*1, 1*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_shade, pred = plotted.pred.shade$y, upper = plotted.pred.shade$upper, lower = plotted.pred.shade$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)

#Fecundity
x_to_plot_sun <- seq.func(seedtrcy$std_logp1_totalabund[seedtrcy$Cover=='Sun'])
x_to_plot_shade <- seq.func(seedtrcy$std_logp1_totalabund[seedtrcy$Cover=='Shade'])
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha(ifelse(Cover=="Sun", "#CC79A7", "#0072B2"), 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedtrcy)
model <- trcyseedfinalmod
#Sun - blue
plotted.pred.sun <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0,0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_sun, pred = plotted.pred.sun$y, upper = plotted.pred.sun$upper, lower = plotted.pred.sun$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Shade - red
plotted.pred.shade <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 1, 0, 0,x_to_plot, 0,1*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_shade, pred = plotted.pred.shade$y, upper = plotted.pred.shade$upper, lower = plotted.pred.shade$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#Lambda
trcy_pred_nonbh<-glmm.predict(mod=trcylambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0,0*0,0*0,c(1,0)*0),
                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
trcy_pred_nonbh$Cover <- c('Shade', 'Sun')
trcy_pred_nonbh$Neighbours01 <- 'Neighbours0'
trcy_pred_nbh<-glmm.predict(mod=trcylambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 1,0*1,0*1,c(1,0)*1),
                            se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
trcy_pred_nbh$Cover <- c('Shade', 'Sun')
trcy_pred_nbh$Neighbours01 <- 'Neighbours1'
trcy_pred <- rbind(trcy_pred_nonbh,trcy_pred_nbh)
trcy_pred$x_label <- c("Shade:Nbh0", "Sun:Nbh0", "Shade:Nbh1", "Sun:Nbh1")
trcy_pred$x <- c(1,0.5,2.5,2)
lambdatrcy$x <- "1"
lambdatrcy <- within(lambdatrcy, x[Cover== 'Sun' & Neighbours01== 'Neighbours1'] <- '2')
lambdatrcy <- within(lambdatrcy, x[Cover== 'Sun' & Neighbours01== 'Neighbours0'] <- '0.5')
lambdatrcy <- within(lambdatrcy, x[Cover== 'Shade' & Neighbours01== 'Neighbours1'] <- '2.5')
lambdatrcy <- within(lambdatrcy, x[Cover== 'Shade' & Neighbours01== 'Neighbours0'] <- '1')
plot(log_lambda ~ x, pch=19, ylim=c(-0.7,3.1),xlim=c(0.25,2.75),col=alpha(ifelse(Cover=="Sun", "#CC79A7", "#0072B2"),0.2), ylab=NA, xlab=NA,  xaxt="n", tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdatrcy)
points(y ~ x, trcy_pred, pch = 17, col=ifelse(Cover=="Sun", "#CC79A7", "#0072B2"), lwd=3,cex= 2.5)
arrows(x0=trcy_pred$x, y0=trcy_pred$lower, x1=trcy_pred$x, y1=trcy_pred$upper, col=c("#0072B2", "#CC79A7", "#0072B2", "#CC79A7"), code=3, angle=90, length=0.1, lwd=3)

#Adding TROR
x_to_plot<-seq.func(trordata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, trordata)
model <- trorsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedtror$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedtror)
model <- trorseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
tror_pred_nonbh<-glmm.predict(mod=trorlambdafinalmod, newdat=data.frame(1, 0, 0, 0, 0, 0, 0),
                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
tror_pred_nonbh$Neighbours01 <- 'Neighbours0'
tror_pred_nbh<-glmm.predict(mod=trorlambdafinalmod, newdat=data.frame(1, 0, 0, 0, 0, 0, 1),
                            se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
tror_pred_nbh$Neighbours01 <- 'Neighbours1'
tror_pred <- rbind(tror_pred_nonbh, tror_pred_nbh)
tror_pred$x <- c(0.5,1)
lambdatror$x <- '1'
lambdatror <- within(lambdatror, x[Neighbours01== 'Neighbours0'] <- '0.5')
lambdatror <- within(lambdatror, x[Neighbours01== 'Neighbours1'] <- '1')
plot(log_lambda ~ x, pch=19, ylim=c(-1.0,1.8),xlim=c(0.25,1.25),col=alpha("grey60",0.3), ylab=NA, xlab=NA, xaxt="n", tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdatror)
points(y ~ x, tror_pred, pch = 17, lwd=3,cex= 2.5)
arrows(x0=tror_pred$x, y0=tror_pred$lower, x1=tror_pred$x, y1=tror_pred$upper, code=3, angle=90, length=0.1, lwd=3)

#Adding VERO
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, verodata)
model <- verosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedvero$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedvero)
model <- veroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, 0,x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
vero_pred_nonbh<-glmm.predict(mod=verolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0,c(1,0)*0),
                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
vero_pred_nonbh$Cover <- c('Shade', 'Sun')
vero_pred_nonbh$Neighbours01 <- 'Neighbours0'
vero_pred_nbh<-glmm.predict(mod=verolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 1,c(1,0)*1),
                            se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
vero_pred_nbh$Cover <- c('Shade', 'Sun')
vero_pred_nbh$Neighbours01 <- 'Neighbours1'
vero_pred <- rbind(vero_pred_nonbh,vero_pred_nbh)
vero_pred$x_label <- c("Shade:Nbh0", "Sun:Nbh0", "Shade:Nbh1", "Sun:Nbh1")
vero_pred$x <- c(1,0.5,2.5,2)
lambdavero$x <- "1"
lambdavero <- within(lambdavero, x[Cover== 'Sun' & Neighbours01== 'Neighbours1'] <- '2')
lambdavero <- within(lambdavero, x[Cover== 'Sun' & Neighbours01== 'Neighbours0'] <- '0.5')
lambdavero <- within(lambdavero, x[Cover== 'Shade' & Neighbours01== 'Neighbours1'] <- '2.5')
lambdavero <- within(lambdavero, x[Cover== 'Shade' & Neighbours01== 'Neighbours0'] <- '1')
plot(log_lambda ~ x, pch=19, ylim=c(-0.9,3.8),xlim=c(0.25,2.75),col=alpha(ifelse(Cover=="Sun", "#CC79A7", "#0072B2"),0.2), ylab=NA, xlab=NA,  xaxt="n", tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdavero)
points(y ~ x, vero_pred, pch = 17, col=ifelse(Cover=="Sun", "#CC79A7", "#0072B2"), lwd=3,cex= 2.5)
arrows(x0=vero_pred$x, y0=vero_pred$lower, x1=vero_pred$x, y1=vero_pred$upper, col=c("#0072B2", "#CC79A7", "#0072B2", "#CC79A7"), code=3, angle=90, length=0.1, lwd=3)
text(x = 1.5, y = 3.2, "*", cex = 8, col = "red")

mtext(side=1, "Absent", adj=0.1, line=2.5, cex =2.5)
mtext(side=1, "Present", adj=0.9, line=2.5, cex =2.5)

###Overall text
##x labels
mtext("Neighbour abundance", adj = 0.09, side = 1, line=3, cex = 3,outer = TRUE)
mtext("Neighbour abundance", adj = 0.54, side = 1, line=3, cex = 3, outer = TRUE)
mtext("Neighbour presence", adj = 0.98, side = 1, line=3, cex = 3, outer = TRUE)
##y labels
mtext("Probability of survival", side = 2, cex = 3, outer=TRUE, line=-6)
mtext("Number of viable seeds produced (log + 1)", side = 2, cex = 3, outer=TRUE, line=-50.5)
mtext("Population growth rate (log)", side = 2, cex = 3, outer=TRUE, line=-97.5)
##main labels
mtext("Survival", outer=TRUE, adj=0.15,side = 3, cex = 3)
mtext("Seed production", outer=TRUE, adj = 0.55, side = 3, cex = 3)
mtext("Population growth", outer=TRUE, adj=0.96, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.14, side = 3, cex = 3)
mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 10, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 20, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 29, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 37, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.16, padj= 35, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 52, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("G. rosea"), adj = -0.15, padj= 59, side = 3, cex = 2.5, outer = TRUE)
reset()
legend("bottom", title=NULL, horiz=T, legend=c("Open", "Shade"),
       col=c("#CC79A7","#0072B2"), pch=19, cex=3, bty="n")
dev.off()

#### For figures 5 and S4: Build dataframe with plot-level vital rates ####
#Need a dataframe with one value for each rate
#vitaldata (germ and survival)
#seedmodeldata (produced at least one seed)
#popdata (lambda) --- this dataframe has data for all rates split by neighbours

#popgrowthrate data has germination by plot, surv and seed production by neighbours by plot
# so does popdata, with lambda

meanvitalrates <- popgrowthratedata
#without pole
meanvitalrates <- meanvitalrates %>% filter(!(Species=="POLE"))
#already have plot_germ, need plot_surv, plot_fec and plot_lambda

##### plot_surv
plot_survdata <- vitaldata %>% filter(!(is.na(surv_to_produce_seeds)))
plot_surv_counts <- plot_survdata %>% group_by(Species, Site, Plot, surv_to_produce_seeds) %>% count()
Species <- rep(c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO"), each = 3, times = 8)
Site <- rep(c("1", "2", "3", "4", "5", "6", "7", "8"), each = 9, times = 3)
Plot <- rep(c("A", "B", "C"), times = 72)
plot_surv_prop <- cbind(Species, Site, Plot)
plot_surv_prop <- data.frame(plot_surv_prop)
plot_surv_prop <- plot_surv_prop %>% unite("idforjoining", Plot:Site, sep = ":", remove = FALSE)
plot_surv_counts$surv_to_produce_seeds <- as.factor(plot_surv_counts$surv_to_produce_seeds)
plot_surv_prop <- left_join(plot_surv_prop, plot_surv_counts)
plot_surv_prop <- within(plot_surv_prop, n[is.na(n)] <- 0)

#
plot_surv <- plot_surv_prop %>% group_by(Species, Site, Plot)  %>% 
  summarise(plot_surv = ifelse(n[surv_to_produce_seeds==1]==0 & n[surv_to_produce_seeds==0]==0, NA, n[surv_to_produce_seeds==1]/(n[surv_to_produce_seeds==1] + n[surv_to_produce_seeds==0])),
            total_n = sum(n))
## Merge with plot-level data
plot_surv <- plot_surv %>% select(Species, Site, Plot, plot_surv)
meanvitalrates <- left_join(meanvitalrates, plot_surv)

#Calculating mean viable seed production of logged seed values (exponentiated)
allplotfec <- viable_plot %>% mutate(log_seeds = log(No_viable_seeds_grouped)) %>%
  group_by(Species, Site, Plot) %>%
  summarise(plot_fec = exp(mean(log_seeds)))
##this is only where there were seeds produced
#Merge with plot-level data
meanvitalrates <- left_join(meanvitalrates, allplotfec)

#### Need plot_lambda (not in presence/absence of neighbours)
### Assigning all survival and fecundity NA value to 0 for lambda calculations
#(no plants germinated, none survived, none produced seeds)
one_poplongdata <- meanvitalrates
one_poplongdata <- within(one_poplongdata, plot_surv[is.na(plot_surv)] <- 0)
one_poplongdata <- within(one_poplongdata, plot_fec[is.na(plot_fec)] <- 0)

#Calculate population growth rates
# Per capita growth rate of a given population  i = seed survival*(1-germination)+number of viable seeds produced per germinant*germination
# prob of germination / germination fraction. # germinated / total germination, currently I have this as percent_germ which is a proportion (not percentage, despite the name)
one_poplongdata <- one_poplongdata %>% mutate(plot_lambda = seed_survival*(1-plot_germ)+plot_fec*plot_surv*plot_germ)
#Still get NAs where seed wasn't sown, working well

#Merge this back with plot-level data
one_poplongdata <- one_poplongdata %>% select(Species, Site, Plot, plot_lambda)
meanvitalrates <- left_join(meanvitalrates, one_poplongdata)

#Simplify main dataframe
plot_vitalrates <- meanvitalrates %>% select(Species, Site, Plot, seed_survival, plot_germ, plot_surv, plot_fec, plot_lambda)

#We need to keep our zero values for fecundity and survival, instead of NAs
#Yes! Need zero values, these are already factored into popdata too
plot_zero <- plot_vitalrates
plot_zero <- within(plot_zero, plot_surv[is.na(plot_surv)] <- 0)
plot_zero <- within(plot_zero, plot_fec[is.na(plot_fec)] <- 0)

## Change the column names, for plotting
colnames(plot_zero) <- c("Species", "Site", "Plot", "seed_survival", "emergence", "survival", "seed production", "population growth rate")

#Split the dataframe by species
plotarca <- plot_zero %>% filter(Species == "ARCA")
plothygl <- plot_zero %>% filter(Species == "HYGL")
plotlaro <- plot_zero %>% filter(Species == "LARO")
plotpeai <- plot_zero %>% filter(Species == "PEAI")
plotplde <- plot_zero %>% filter(Species == "PLDE")
plottrcy <- plot_zero %>% filter(Species == "TRCY")
plottror <- plot_zero %>% filter(Species == "TROR")
plotvero <- plot_zero %>% filter(Species == "VERO")

#Need just the rates for plotting
corrarca <- plotarca %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrhygl <- plothygl %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrlaro <- plotlaro %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrpeai <- plotpeai %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrplde <- plotplde %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrtrcy <- plottrcy %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrtror <- plottror %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrvero <- plotvero %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)

#overall (all species)
justrates <- plot_zero %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)

# logged version
# logging seed production and lambda for all species
# seed production needs to be +1 then logged
plot_zero_log <- plot_zero %>% mutate(log_seed = log(`seed production`+1), log_lambda = log(`population growth rate`))

#Split the dataframe by species
log_plotarca <- plot_zero_log %>% filter(Species == "ARCA")
log_plothygl <- plot_zero_log %>% filter(Species == "HYGL")
log_plotlaro <- plot_zero_log %>% filter(Species == "LARO")
log_plotpeai <- plot_zero_log %>% filter(Species == "PEAI")
log_plotplde <- plot_zero_log %>% filter(Species == "PLDE")
log_plottrcy <- plot_zero_log %>% filter(Species == "TRCY")
log_plottror <- plot_zero_log %>% filter(Species == "TROR")
log_plotvero <- plot_zero_log %>% filter(Species == "VERO")

#Need just the rates for plotting
log_corrarca <- log_plotarca %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_corrhygl <- log_plothygl %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_corrlaro <- log_plotlaro %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_corrpeai <- log_plotpeai %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_corrplde <- log_plotplde %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_corrtrcy <- log_plottrcy %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_corrtror <- log_plottror %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_corrvero <- log_plotvero %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)

#overall (all species)
log_justrates <- plot_zero_log %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)

#### For figures 5 and S4: Making unmanipulated plot data ####
## In the presence of neighbours, ambient plots only.
# Treatment = Ambient

#popgrowthrate data has germination by plot, surv and seed production by neighbours by plot
# so does popdata, with lambda
#already have plot_germ, need plot_surv, plot_fec and plot_lambda

unmani_plot <- popgrowthratedata
#without pole
unmani_plot <- unmani_plot %>% filter(!(Species=="POLE"))
#with neighbours only
unmani_plot <- unmani_plot %>% select(Species, Site, Plot,seed_survival,plot_germ, plot_surv=surv_prop_neighbours, plot_fec=fecundity_nbh)
#bringing in lambda
#with neighbours only
unmani_lambda <- popdata %>% filter(Neighbours01=="Neighbours1")
unmani_lambda <- unmani_lambda %>% ungroup() %>% select(Species, Site, Plot, Treatment, plot_lambda=lambda)
#joining
unmani_rates <- left_join(unmani_plot,unmani_lambda)
#watering plots only
unmani_rates <- unmani_rates %>% filter(Treatment=="Ambient")
#remove PLDE 5,6,7,8 (PLDE seeds were never sown in these four plots)
unmani_rates <- unmani_rates %>% filter(!(is.na(plot_germ)))

### Assigning all survival and fecundity NA values to 0
unmani_rates <- within(unmani_rates, plot_surv[is.na(plot_surv)] <- 0)
unmani_rates <- within(unmani_rates, plot_fec[is.na(plot_fec)] <- 0)

#Simplify main dataframe and rename columns
unmani_rates <- unmani_rates %>% select(Species, Site, Plot, seed_survival, plot_germ, plot_surv, plot_fec, plot_lambda)
colnames(unmani_rates) <- c("Species", "Site", "Plot", "seed_survival", "emergence", "survival", "seed production", "population growth rate")

#Split the dataframe by species
unmani_arca <- unmani_rates %>% filter(Species == "ARCA")
unmani_hygl <- unmani_rates %>% filter(Species == "HYGL")
unmani_laro <- unmani_rates %>% filter(Species == "LARO")
unmani_peai <- unmani_rates %>% filter(Species == "PEAI")
unmani_plde <- unmani_rates %>% filter(Species == "PLDE")
unmani_trcy <- unmani_rates %>% filter(Species == "TRCY")
unmani_tror <- unmani_rates %>% filter(Species == "TROR")
unmani_vero <- unmani_rates %>% filter(Species == "VERO")

#Need just the rates for plotting
unmani_corrarca <- unmani_arca %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
unmani_corrhygl <- unmani_hygl %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
unmani_corrlaro <- unmani_laro %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
unmani_corrpeai <- unmani_peai %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
unmani_corrplde <- unmani_plde %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
unmani_corrtrcy <- unmani_trcy %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
unmani_corrtror <- unmani_tror %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
unmani_corrvero <- unmani_vero %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)

#overall (all species)
just_unmani_rates <- unmani_rates %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)

### and for logged version
# logging seed production and lambda for all species
# seed production needs to be +1 then logged
unmani_rates_log <- unmani_rates %>% mutate(log_seed = log(`seed production`+1), log_lambda = log(`population growth rate`))

#Split the dataframe by species
log_unmani_plotarca <- unmani_rates_log %>% filter(Species == "ARCA")
log_unmani_plothygl <- unmani_rates_log %>% filter(Species == "HYGL")
log_unmani_plotlaro <- unmani_rates_log %>% filter(Species == "LARO")
log_unmani_plotpeai <- unmani_rates_log %>% filter(Species == "PEAI")
log_unmani_plotplde <- unmani_rates_log %>% filter(Species == "PLDE")
log_unmani_plottrcy <- unmani_rates_log %>% filter(Species == "TRCY")
log_unmani_plottror <- unmani_rates_log %>% filter(Species == "TROR")
log_unmani_plotvero <- unmani_rates_log %>% filter(Species == "VERO")

#Need just the rates for plotting
log_unmani_corrarca <- log_unmani_plotarca %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_unmani_corrhygl <- log_unmani_plothygl %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_unmani_corrlaro <- log_unmani_plotlaro %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_unmani_corrpeai <- log_unmani_plotpeai %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_unmani_corrplde <- log_unmani_plotplde %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_unmani_corrtrcy <- log_unmani_plottrcy %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_unmani_corrtror <- log_unmani_plottror %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)
log_unmani_corrvero <- log_unmani_plotvero %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)

#overall (all species)
log_unmani_justrates <- unmani_rates_log %>% ungroup() %>% select(emergence, survival, log_seed, log_lambda)

### Correlations for data with only no neighbours
no_nbh_plot <- popgrowthratedata
#without pole
no_nbh_plot <- no_nbh_plot %>% filter(!(Species=="POLE"))
#plots without neighbours only
no_nbh_plot <- no_nbh_plot %>% select(Species, Site, Plot,seed_survival,plot_germ, plot_surv=surv_prop_no_neighbours, plot_fec=fecundity_no_nbh)
#bringing in lambda
#plots without neighbours only
no_nbh_lambda <- popdata %>% filter(Neighbours01=="Neighbours0")
no_nbh_lambda <- no_nbh_lambda %>% ungroup() %>% select(Species, Site, Plot, Treatment, plot_lambda=lambda)
#joining
no_nbh_rates <- left_join(no_nbh_plot,no_nbh_lambda)
#remove PLDE 5,6,7,8
no_nbh_rates <- no_nbh_rates %>% filter(!(is.na(plot_germ)))

### Assigning all survival and fecundity NA values to 0
no_nbh_rates <- within(no_nbh_rates, plot_surv[is.na(plot_surv)] <- 0)
no_nbh_rates <- within(no_nbh_rates, plot_fec[is.na(plot_fec)] <- 0)

#Simplify main dataframe and rename columns
no_nbh_rates <- no_nbh_rates %>% select(Species, Site, Plot, seed_survival, plot_germ, plot_surv, plot_fec, plot_lambda)
colnames(no_nbh_rates) <- c("Species", "Site", "Plot", "seed_survival", "emergence", "survival", "seed production", "population growth rate")

#Split the dataframe by species
no_nbh_arca <- no_nbh_rates %>% filter(Species == "ARCA")
no_nbh_hygl <- no_nbh_rates %>% filter(Species == "HYGL")
no_nbh_laro <- no_nbh_rates %>% filter(Species == "LARO")
no_nbh_peai <- no_nbh_rates %>% filter(Species == "PEAI")
no_nbh_plde <- no_nbh_rates %>% filter(Species == "PLDE")
no_nbh_trcy <- no_nbh_rates %>% filter(Species == "TRCY")
no_nbh_tror <- no_nbh_rates %>% filter(Species == "TROR")
no_nbh_vero <- no_nbh_rates %>% filter(Species == "VERO")

#Need just the rates for plotting
no_nbh_corrarca <- no_nbh_arca %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
no_nbh_corrhygl <- no_nbh_hygl %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
no_nbh_corrlaro <- no_nbh_laro %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
no_nbh_corrpeai <- no_nbh_peai %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
no_nbh_corrplde <- no_nbh_plde %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
no_nbh_corrtrcy <- no_nbh_trcy %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
no_nbh_corrtror <- no_nbh_tror %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
no_nbh_corrvero <- no_nbh_vero %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)

#overall (all species)
just_no_nbh_rates <- no_nbh_rates %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)

### Figure 5 - logged plots all subplots ####
dev.off()
pdf("Output/Figures/panel_all_log_corr_2024.pdf", width = 21, height = 12)
par(mfrow=c(2,4), oma =c(1,1,1,1))

testRes = cor.mtest(log_corrarca, conf.level = 0.95)
arca_all = cor(log_corrarca, use = 'complete')
corrplot(arca_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#hygl
testRes = cor.mtest(log_corrhygl, conf.level = 0.95)
hygl_all = cor(log_corrhygl, use = 'complete')
corrplot(hygl_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#laro
testRes = cor.mtest(log_corrlaro, conf.level = 0.95)
laro_all = cor(log_corrlaro, use = 'complete')
corrplot(laro_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#peai
testRes = cor.mtest(log_corrpeai, conf.level = 0.95)
peai_all = cor(log_corrpeai, use = 'complete')
corrplot(peai_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#plde
testRes = cor.mtest(log_corrplde, conf.level = 0.95)
plde_all = cor(log_corrplde, use = 'complete')
corrplot(plde_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#trcy
testRes = cor.mtest(log_corrtrcy, conf.level = 0.95)
trcy_all = cor(log_corrtrcy, use = 'complete')
corrplot(trcy_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#tror
testRes = cor.mtest(log_corrtror, conf.level = 0.95)
tror_all = cor(log_corrtror, use = 'complete')
corrplot(tror_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#vero
testRes = cor.mtest(log_corrvero, conf.level = 0.95)
vero_all = cor(log_corrvero, use = 'complete')
corrplot(vero_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
dev.off()

### Figure S4 - corr plots no neighbours ####
dev.off()
pdf("Output/Figures/panel_corr_no_nbh_2024.pdf", width = 21, height = 12)
par(mfrow=c(2,4), oma =c(1,1,1,1))

testRes = cor.mtest(no_nbh_corrarca, conf.level = 0.95)
arca_all = cor(no_nbh_corrarca, use = 'complete')
corrplot(arca_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#hygl
testRes = cor.mtest(no_nbh_corrhygl, conf.level = 0.95)
hygl_all = cor(no_nbh_corrhygl, use = 'complete')
corrplot(hygl_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#laro
testRes = cor.mtest(no_nbh_corrlaro, conf.level = 0.95)
laro_all = cor(no_nbh_corrlaro, use = 'complete')
corrplot(laro_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#peai
testRes = cor.mtest(no_nbh_corrpeai, conf.level = 0.95)
peai_all = cor(no_nbh_corrpeai, use = 'complete')
corrplot(peai_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#plde
testRes = cor.mtest(no_nbh_corrplde, conf.level = 0.95)
plde_all = cor(no_nbh_corrplde, use = 'complete')
corrplot(plde_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#trcy
testRes = cor.mtest(no_nbh_corrtrcy, conf.level = 0.95)
trcy_all = cor(no_nbh_corrtrcy, use = 'complete')
corrplot(trcy_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#tror
testRes = cor.mtest(no_nbh_corrtror, conf.level = 0.95)
tror_all = cor(no_nbh_corrtror, use = 'complete')
corrplot(tror_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#vero
testRes = cor.mtest(no_nbh_corrvero, conf.level = 0.95)
vero_all = cor(no_nbh_corrvero, use = 'complete')
corrplot(vero_all, method = 'square', type = 'lower', p.mat = testRes$p, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
dev.off()

