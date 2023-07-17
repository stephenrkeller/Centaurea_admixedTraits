# FIXED admix traits over diffK
# July 2023
# code by: SRK & ZMP

# load packages
library(tidyverse)
library(lmerTest)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(forcats)

# load data with traits and diffK
admix <- read.csv(file = 'data/traitsAdmix.csv', header = T, sep = ',', stringsAsFactors = T)

admix$diffK1sq <- admix$diffK1^2
admix$totFl = admix$noCapitula + admix$noFlowers
admix$DFF=1/admix$dFF # convert to inverse of days to flowering (=flowering speed)
admix$SLA = (3*pi*1.5^2)/admix$SLAwt # 3 discs; convert to units of mm2.mg−1
admix$reprodStat <- ifelse(is.na(admix$dFF),"reprod","veget") # define reproductive status

length(levels(admix$cross)) # 46 crosstypes
length(levels(admix$damInd)) # 37 dams
length(levels(admix$sireInd)) # 30 sires

# Make heat map of crossing design

admixplot = unique(admix[c("cross","damInd","sireInd","damK1","sireK1","diffK1")])
admixplot2 = na.omit(admixplot)
admixplot2$sireInd = droplevels(admixplot2)$sireInd

# design <- ggplot(admix, aes(sireInd, damInd)) +                           # Create heatmap with ggplot2
#   geom_raster(hjust=0, vjust=0,
#             aes(fill = diffK1)) +
#             scale_fill_gradient(low = "yellow", high = "red") +
#             theme_grey(base_size=8)
# 
design <- ggplot(admixplot, aes(sireInd, damInd)) +
  geom_tile(lwd=0.25,aes(fill = diffK1)) +
    theme_grey(base_size=6) +
    scale_fill_gradient(low = "yellow", high = "red") +
    scale_x_discrete(limits=admixplot$sireInd[order(admixplot$sireK1)], guide = guide_axis(angle = 45)) +
    scale_y_discrete(limits=admixplot$damInd[order(admixplot$damK1)]) +
    #theme(element_blank()) +
    coord_fixed()

design 


  

ggsave(design, filename="figs/Crossing_design.pdf", device="pdf",height=5.5, width=10, units="in", dpi=200)
################################
# making BLUPs
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

# # num capitula
# parentNCap <-lmer(noCapitula ~ benchPos + planting + recorder + (1|cross), admix)
# pNCap <- ranef(parentNCap)
# pNCap <- pNCap[['cross']]
# 
# pNCap <- cbind(fam = rownames(pNCap), pNCap)
# names(pNCap) <- c("cross","NCap")
# 
# # num flowers
# parentNFl <-lmer(noFlowers ~ benchPos + planting + recorder + (1|cross), admix)
# pNFl <- ranef(parentNFl)
# pNFl <- pNFl[['cross']]
# 
# pNFl <- cbind(fam = rownames(pNFl), pNFl)
# names(pNFl) <- c("cross","NFl")

# Total flower heads
parenttotFl <-lmer(totFl ~ benchPos + planting + recorder + (1|cross), admix)
ptotFl <- ranef(parenttotFl)
ptotFl <- ptotFl[['cross']]

ptotFl <- cbind(fam = rownames(ptotFl), ptotFl)
names(ptotFl) <- c("cross","totFl")

# DFF
parentdff <-lmer(DFF ~ benchPos + planting + recorder + (1|cross), admix)
pdff <- ranef(parentdff)
pdff <- pdff[['cross']]

pdff <- cbind(fam = rownames(pdff), pdff)
names(pdff) <- c("cross","dff")

# SLA
parentsla <-lmer(SLA ~ benchPos + planting + recorder + (1|cross), admix)
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
  left_join(pdff, by = "cross") %>%
  left_join(psla, by="cross") %>%
  left_join(ptotFl, by = "cross") 
  #left_join(pNCap, by = "cross") %>%
  #left_join(pNFl, by = "cross") %>%


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
    geom_point(size=3, pch=21, aes(fill=diffK1)) +
    scale_fill_gradient(low = "yellow", high = "red") + 
    stat_smooth(method = "lm", col = "black") + 
    theme(legend.position="none") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       " P =",signif(summary(fit)$coef[2,4], 5))
    #                    "Intercept =",signif(fit$coef[[1]],5 ),
    #                    " Slope =",signif(fit$coef[[2]], 5),
                        )
}

# ggplotRegression2 <- function (fit) {
#   require(ggplot2)
#   ggplot(fit@frame, aes_string(x = names(fit@frame)[2], y = names(fit@frame)[1])) +
#     geom_point(size=3, pch=21, aes(fill=diffK1)) +
#     scale_fill_gradient(low = "yellow", high = "red") +
#     stat_smooth(method = "lm", col = "black") +
#     theme(legend.position="none") +
#     labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
#                        " P =",signif(summary(fit)$coef[2,4], 5)))
#   #                    "Intercept =",signif(fit$coef[[1]],5 ),
#   #                    " Slope =",signif(fit$coef[[2]], 5),
# }

# ggplotRegressionQuad <- function (fit) {
#   require(ggplot2)
#   ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
#     geom_point(size=3) +
#     stat_smooth(method = "lm", formula = y ~ x + I(x^2), col = "red")
#     # labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
#     #                    "Intercept =",signif(fit$coef[[1]],5 ),
#     #                    " Slope =",signif(fit$coef[[2]], 5),
#     #                    " P =",signif(summary(fit)$coef[2,4], 5)))
# }
####################################
# lm with control for planting, benchPos, recorder

# heightK1 (significant)
# lmquad <- lm(height ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# lm <- lm(height ~ diffK1 + benchPos + planting + recorder, data = admix)
lm <- lmer(scale(height) ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
lmquad <- lmer(scale(height) ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), data = admix)
anova(lm,lmquad)
summary(lm)
admix$s.height = scale(admix$height)
#lmplot <- lm(ht ~ diffK1mean, data = diffkBLUPs)
lmplot <- lm(s.height ~ diffK1 + benchPos + planting + recorder, data = admix)

ht <- ggplotRegression(lmplot) +
  #geom_text(x=0.7,y=123,label = "height ~ diffK1 + benchPos + planting + recorder") +
  xlab("Diff-K") +
  ylab("Height")

#ggsave(file = "heightDiffK1.png", plot = last_plot())


# stem width
# lmquad <- lm(stemWidth ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# lm2 <- lm(stemWidth ~ diffK1 + benchPos + planting + recorder, data = admix)
lm2 <- lmer(stemWidth ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
lmquad2 <- lmer(stemWidth ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), data = admix)
anova(lm2,lmquad2)
summary(lm2)
admix$s.stemWidth = scale(admix$stemWidth)
lmplot2 <- lm(s.stemWidth ~ diffK1 + benchPos + planting + recorder, data = admix)

sw <- ggplotRegression(lmplot2) +
  #geom_text(x=0.3,y=8.5,label = "stemWidth ~ diffK1 + benchPos + planting + recorder") +
  xlab("Diff-K") +
  ylab("Stem Width")

#ggsave(file = "stemWidthDiffK1.png", plot = last_plot()) 

# biomass
# lmquad3 <- lm(weight ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# lm3 <- lm(weight ~ diffK1 + benchPos + planting + recorder, data = admix)
lm3 <- lmer(weight ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
lmquad3 <- lmer(weight ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), data = admix)
anova(lm3,lmquad3)
summary(lm3)
admix$s.weight = scale(admix$weight)
lmplot3 <- lm(s.weight ~ diffK1 + benchPos + planting + recorder, data = admix)

dm <- ggplotRegression(lmplot3) +
  #geom_text(x=0.3,y=1,label = "weight ~ diffK1 + benchPos + planting + recorder") +
  xlab("Diff-K") +
  ylab("Dry Mass")

#ggsave(file = "weightDiffK1.png", plot = last_plot())

# SLA
# lmquad4 <- lm(SLAwt ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# lm4 <- lm(SLAwt ~ diffK1 + benchPos + planting + recorder, data = admix)
lm4 <- lmer(SLAwt ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
lmquad4 <- lmer(SLAwt ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), data = admix)
anova(lm4,lmquad4)
summary(lm4)
admix$s.SLAwt = scale(admix$SLAwt)
lmplot4 <- lm(s.SLAwt ~ diffK1 + benchPos + planting + recorder, data = admix)

sla <- ggplotRegression(lmplot4) +
  #geom_text(x=0.75,y=10.5,label = "SLAwt ~ diffK1 + benchPos + planting + recorder") +
  xlab("Diff-K") +
  ylab("SLA")

#ggsave(file = "SLADiffK1.png", plot = last_plot())

# time to first flower
# lmquad5 <- lm(dFF ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# lm5 <- lm(dFF ~ diffK1 + benchPos + planting + recorder, data = admix)
lm5 <- lmer(DFF ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
lmquad5 <- lmer(DFF ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), data = admix)
anova(lm5,lmquad5)
summary(lm5)
admix$s.dFFspd = scale(1/admix$dFF)
lmplot5 <- lm(s.dFFspd ~ diffK1 + benchPos + planting + recorder, data = admix)

dFF <- ggplotRegression(lmplot5) +
  #geom_text(x=0.7,y=81,label = "dFF ~ diffK1 + benchPos + planting + recorder") +
  xlab("Diff-K") +
  ylab("1 / DFF")

#ggsave(file = "dFFDiffK1.png", plot = last_plot())

# noCapitula ### significant
# lmquad6 <- lm(noCapitula ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# lm6 <- lm(noCapitula ~ diffK1 + benchPos + planting + recorder, data = admix)
# lm6 <- lmer(noCapitula ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
# lmquad6 <- lmer(noCapitula ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), data = admix)
# anova(lm6,lmquad6)
# summary(lm6)
# admix$s.noCapitula = scale(admix$noCapitula)
# lmplot6 <- lm(s.noCapitula ~ diffK1 + benchPos + planting + recorder, data = admix)
# 
# noCapitula <- ggplotRegression(lmplot6) +
#   #geom_text(x=0.7,y=45,label = "height ~ diffK1 + diffK1sq benchPos + planting + recorder") +
#   xlab("Diff-K") +
#   ylab("Num. Capitula")
# 
# #ggsave(file = "noCapDiffK1.png", plot = last_plot())
# 
# ## noFlowers
# # lmquad7 <- lm(noFlowers ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# # lm7 <- lm(noFlowers ~ diffK1 + benchPos + planting + recorder, data = admix)
# lm7 <- lmer(noFlowers ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
# lmquad7 <- lmer(noFlowers ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), data = admix)
# anova(lm7,lmquad7)
# summary(lm7)
# admix$s.noFlowers = scale(admix$noFlowers)
# lmplot7 <- lm(s.noFlowers ~ diffK1 + benchPos + planting + recorder, data = admix)
# 
# noFlowers <- ggplotRegression(lmplot7) +
#   #geom_text(x=0.7,y=40,label = "noFlowers ~ diffK1 + benchPos + planting + recorder") +
#   xlab("Diff-K") +
#   ylab("Num. Flowers")
# 
# #ggsave(file = "noFlDiffK1.png", plot = last_plot())

#### Total Flowers
# lmquad8 <- lm(totFl ~ diffK1 + diffK1sq + benchPos + planting + recorder, data = admix)
# lm8 <- lm(totFl ~ diffK1 + benchPos + planting + recorder, data = admix)
lm8 <- lmer(totFl ~ diffK1 + benchPos + planting + recorder + (1|cross), data = admix)
lmquad8 <- lmer(totFl ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), data = admix)
anova(lm8,lmquad8)
summary(lm8)
admix$s.totFl = scale(admix$totFl)
lmplot8 <- lm(s.totFl ~ diffK1 + benchPos + planting + recorder, data = admix)

totFl <- ggplotRegression(lmplot8) +
  #geom_text(x=0.7,y=47,label = "totFl ~ diffK1 + benchPos + planting + recorder") +
  xlab("Diff-K") +
  ylab("Tot. Flower Heads")

#ggsave(file = "totFlDiffK1.png", plot = last_plot())


# lm9 <- glmer(as.factor(reprodStat) ~ diffK1 + benchPos + planting + recorder + (1|cross), family="binomial", data = admix)
# lmquad9 <- glmer(as.factor(reprodStat) ~ diffK1 + I(diffK1^2) + benchPos + planting + recorder + (1|cross), family="binomial", data = admix)
# anova(lm9,lmquad9)
# summary(lm9)
# 

multpanel_diffK <- marrangeGrob(list(ht,sw,dm,sla,dFF,totFl), nrow=2,ncol=3, top=NULL, padding = unit(10, "line"))

ggsave("figs/Figure2.pdf", multpanel_diffK, width=14,height=8)

# ANOVA table
F1.Model.resultsp = list()
admix2 = admix[,c(14,15,10,26,25,24,20,2,3,22,4)]
names(admix2) <- c("Height","Stem Width","Biomass","SLA","1/DFF","Tot. Flower Heads","diffK1","benchPos","planting","recorder","cross")

for(i in 1:6) {
  admix2$trait = scale(as.numeric(admix2[,i]))
  tmp.lm = lmer(trait ~ diffK1 + benchPos + planting + (1|cross), data = admix2)
  F1.Model.resultsp[[i]] = list()
  F1.Model.resultsp[[i]][[1]] = tmp.lm
}

# Make table of model estimates across traits
tab_model(F1.Model.resultsp[[6]][[1]],
          F1.Model.resultsp[[5]][[1]],
          F1.Model.resultsp[[4]][[1]],
          F1.Model.resultsp[[3]][[1]],
          F1.Model.resultsp[[2]][[1]],
          F1.Model.resultsp[[1]][[1]],
          dv.labels=c("Height","Stem Width","Biomass","SLA","1/DFF","Tot. Flower Heads"),
#          pred.labels=c("(Intercept)","diff-K","bench","planting","recorder[a]","recorder[b]","recorder[c]","recorder[d]"),
          pred.labels=c("(Intercept)","diff-K","bench","planting"),
          string.se = "SE",
          show.ci=F, show.se=T)


