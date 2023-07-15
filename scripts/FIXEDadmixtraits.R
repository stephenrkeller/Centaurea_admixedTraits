# FIXED admix traits over diffK using BLUPs
# June 20 2019
# ZMP; updated July 2023 by SRK

# load packages
library(tidyverse)
library(lmerTest)
library(sjPlot)

# load data with traits and diffK
admix <- read.csv(file = 'data/traitsAdmix.csv', header = T, sep = ',', stringsAsFactors = T)

admix$diffK1sq <- admix$diffK1^2

length(levels(admix$cross)) # 46 crosstypes
length(levels(admix$damInd)) # 37 dams
length(levels(admix$sireInd)) # 30 sires

admix$reprodStat <- ifelse(is.na(admix$dFF),"reprod","veget") # define reproductive status


################################
# making BLUPs just in case
# height
admixHt <-lmer(height ~ benchPos + planting + recorder + (1|cross), admix)
summary(admixHt)

pHt <- ranef(admixHt)
pHt <- pHt[['cross']]

pHt <- cbind(fam = rownames(pHt), pHt)
names(pHt) <- c("cross","ht")

# weight
parentWt <-lmer(weight ~ benchPos + planting + recorder + (1|cross), admix)
pWt <- ranef(parentWt)
pWt <- pWt[['cross']]

pWt <- cbind(fam = rownames(pWt), pWt)
names(pWt) <- c("cross","wt")

# stem width
parentSW <-lmer(stemWidth ~ benchPos + planting + recorder + (1|cross), admix)
pSW <- ranef(parentSW)
pSW <- pSW[['cross']]

pSW <- cbind(fam = rownames(pSW), pSW)
names(pSW) <- c("cross","stemwid")

# num capitula
parentNCap <-lmer(noCapitula ~ benchPos + planting + recorder + (1|cross), admix)
pNCap <- ranef(parentNCap)
pNCap <- pNCap[['cross']]

pNCap <- cbind(fam = rownames(pNCap), pNCap)
names(pNCap) <- c("cross","NCap")

# num flowers
parentNFl <-lmer(noFlowers ~ benchPos + planting + recorder + (1|cross), admix)
pNFl <- ranef(parentNFl)
pNFl <- pNFl[['cross']]

pNFl <- cbind(fam = rownames(pNFl), pNFl)
names(pNFl) <- c("cross","NFl")

# dFF
parentdff <-lmer(dFF ~ benchPos + planting + recorder + (1|cross), admix)
pdff <- ranef(parentdff)
pdff <- pdff[['cross']]

pdff <- cbind(fam = rownames(pdff), pdff)
names(pdff) <- c("cross","dff")

# SLA
parentsla <-lmer(SLAwt ~ benchPos + planting + recorder + (1|cross), admix)
psla <- ranef(parentsla)
psla <- psla[['cross']]

psla <- cbind(fam = rownames(psla), psla)
names(psla) <- c("cross","sla")

#### Combining into dataframe
admixBLUPs<-left_join(pHt,pWt, by = "cross")
# ,pSW,pNCap,pNFl,pdff,psla
admixBLUPs<-admixBLUPs %>%
  left_join(pSW, by = "cross")

admixBLUPs<-admixBLUPs %>%
  left_join(pNCap, by = "cross") %>%
  left_join(pNFl, by = "cross") %>%
  left_join(pdff, by = "cross") %>%
  left_join(psla, by="cross")

K <- admix %>%
  group_by(cross) %>%
  summarise(mean(diffK1), mean(diffK1sq))
colnames(K) <- c("cross", "diffK1mean", "diffK1sqmean")

diffkBLUPs <- merge(K,admixBLUPs,by="cross")
write.csv(file = "data/diffkBLUPs.csv", diffkBLUPs)
####################################
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point(size=3) +
    stat_smooth(method = "lm", col = "red")
    # labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
    #                    "Intercept =",signif(fit$coef[[1]],5 ),
    #                    " Slope =",signif(fit$coef[[2]], 5),
    #                    " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegressionQuad <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point(size=3) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), col = "red")
    # labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
    #                    "Intercept =",signif(fit$coef[[1]],5 ),
    #                    " Slope =",signif(fit$coef[[2]], 5),
    #                    " P =",signif(summary(fit)$coef[2,4], 5)))
}
####################################
# lm with control for planting, benchPos, recorder


# heightK1 (significant)
# lmquad <- lm(height ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# lm <- lm(height ~ diffK1 + benchPos + planting + recorder, data = admix)
lm <- lmer(height ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
lmquad <- lm(ht ~ diffK1mean + diffK1sqmean, data = diffkBLUPs)
anova(lm,lmquad)
summary(lm)

lmplot <- lm(ht ~ diffK1mean, data = diffkBLUPs)

ggplotRegression(lmplot) +
  #geom_text(x=0.7,y=123,label = "height ~ diffK1 + benchPos + planting + recorder") +
  xlab("Difference between K1 of Dam and Sire") +
  ylab("Height")

ggsave(file = "heightDiffK1.png", plot = last_plot())


# stem width
lmquad <- lm(stemWidth ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
lm2 <- lm(stemWidth ~ diffK1 + benchPos + planting + recorder, data = admix)
anova(lm2,lmquad)

summary(lm2)
summary(lmquad)
ggplotRegressionQuad(lmquad) +
  geom_text(x=0.3,y=8.5,label = "stemWidth ~ diffK1 + benchPos + planting + recorder") +
  xlab("Difference between K1 of Dam and Sire") +
  ylab("Stem Width")

ggplot(data = admix, mapping = aes(x = diffK1, y = stemWidth))+
  geom_point() +
  geom_smooth(method = "lm", col = "black") +
  labs(x = "Difference in K of parents (diffK)", y = "Stem width") +
  ylim(0,10) +
  theme_bw()


ggsave(file = "stemWidthDiffK1.png", plot = last_plot()) 
# biomass
# lmquad3 <- lm(weight ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
lm3 <- lm(weight ~ diffK1 + benchPos + planting + recorder, data = admix)
# anova(lm3,lmquad3)

summary(lm3)
ggplotRegression(lm3) +
  geom_text(x=0.3,y=1,label = "weight ~ diffK1 + benchPos + planting + recorder") +
  xlab("Difference between K1 of Dam and Sire") +
  ylab("Weight")

ggsave(file = "weightDiffK1.png", plot = last_plot())

# SLA
# lmquad4 <- lm(SLAwt ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
lm4 <- lm(SLAwt ~ diffK1 + benchPos + planting + recorder, data = admix)
# anova(lm4,lmquad4)

summary(lm4)
ggplotRegression(lm4) +
  geom_text(x=0.75,y=10.5,label = "SLAwt ~ diffK1 + benchPos + planting + recorder") +
  xlab("Difference between K1 of Dam and Sire") +
  ylab("SLA")

ggsave(file = "SLADiffK1.png", plot = last_plot())

# time to first flower
# lmquad5 <- lm(dFF ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
lm5 <- lm(dFF ~ diffK1 + benchPos + planting + recorder, data = admix)
# anova(lm5,lmquad5)

summary(lm5)
ggplotRegression(lm5) +
  geom_text(x=0.7,y=81,label = "dFF ~ diffK1 + benchPos + planting + recorder") +
  xlab("Difference between K1 of Dam and Sire") +
  ylab("dFF")

ggsave(file = "dFFDiffK1.png", plot = last_plot())

# noCapitula ### significant
lmquad6 <- lm(noCapitula ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
lm6 <- lm(noCapitula ~ diffK1 + benchPos + planting + recorder, data = admix)
anova(lm6,lmquad6)

summary(lmquad6)
ggplotRegressionQuad(lmquad6) +
  geom_text(x=0.7,y=45,label = "height ~ diffK1 + diffK1sq benchPos + planting + recorder") +
  xlab("Difference between K1 of Dam and Sire") +
  ylab("Number of Capitula")

ggsave(file = "noCapDiffK1.png", plot = last_plot())

## noFlowers
# lmquad7 <- lm(noFlowers ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
lm7 <- lm(noFlowers ~ diffK1 + benchPos + planting + recorder, data = admix)
# anova(lm7,lmquad7)

summary(lm7)
ggplotRegression(lm7) +
  geom_text(x=0.7,y=40,label = "noFlowers ~ diffK1 + benchPos + planting + recorder") +
  xlab("Difference between K1 of Dam and Sire") +
  ylab("Number of Flowers")

ggsave(file = "noFlDiffK1.png", plot = last_plot())

#### Total Flowers
admix$totFl <- admix$noCapitula + admix$noFlowers

# lmquad8 <- lm(totFl ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
lm8 <- lm(totFl ~ diffK1 + benchPos + planting + recorder, data = admix)
# anova(lm8,lmquad8)

summary(lm8)
ggplotRegression(lm8) +
  geom_text(x=0.7,y=47,label = "totFl ~ diffK1 + benchPos + planting + recorder") +
  xlab("Difference between K1 of Dam and Sire") +
  ylab("Total Number of Flowers")

ggsave(file = "totFlDiffK1.png", plot = last_plot())