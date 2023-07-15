library(ggplot2)
library(gridExtra)
library(lmerTest)
library(sjPlot)

### Data import and merging ###
cent <- read.csv("data/Data_Greenhouse_Flowcyt_Corrected_SRK.csv", header=T)
cent$totFl = cent$numcap + cent$numfls
cent$dff=1/cent$dff # convert to inverse of days to flowering (=flowering speed)
cent$SLA = (3*pi*1.5^2)/cent$SLA # 3 discs; convert to units of mm2.mg−1
#cent$totFl = ifelse(cent$BIOMASS>0 & is.na(cent$totFl)==T, 0, cent$totFl)

gen <- read.table("data/GeneticData_180320_corrected.txt", header=T)

merged <- merge(cent,gen, by.x="fam2", by.y="Ind")

### Trait modeling and plotting ###
dat2 = merged[,c(1,4,11,12,17,18,15,19,13,14,21,22,16)]
names(dat2) = c("Ind","Pop","Height","Stem_Width","Biomass","SLA","1/DFF","Tot. Flower Heads","Num_Capitula","Num_Flowers","AnK2_K1","AnK2_K2","bench")

# Genome size ~ trait relationships
datGS = merged[,c(1,4,7,11,12,17,18,15,19,13,14,16)]
names(datGS) = c("Ind","Pop","GS","Height","Stem_Width","Biomass","SLA","1/DFF","Tot. Flower Heads","Num_Capitula","Num_Flowers","bench")

# Genome size ~ trait LMMs
GSModel.results = list()

for(i in 4:9) {
  tmp = datGS[which(datGS[i]!='NA'),]
  tmp$trait = scale(tmp[,i])
  tmp$Ind = as.factor(tmp$Ind)
  tmp$Pop = as.factor(tmp$Pop)
  tmp$bench = as.factor(tmp$bench)
  tmp.GS = lmer(trait ~ GS + bench + (1|Ind) + (1|Pop), data=tmp)
  GSModel.results[[i-3]] = tmp.GS
}

tab_model(GSModel.results[[6]],
          GSModel.results[[5]],
          GSModel.results[[4]],
          GSModel.results[[3]],
          GSModel.results[[2]],
          GSModel.results[[1]],
          dv.labels=c("Height","Stem Width","Biomass","SLA","1 / DFF","Tot. Flower Heads"),
          string.se = "SE",
          show.ci=F, show.se=T)

# Sample sizes of families and individuals by trait
N_fams = data.frame()
for(i in 3:8){
  trait = names(dat2[i])
  nfams = length(unique(dat2[which(dat2[i]!='NA'),"Ind"]))
  ninds = length(dat2[which(dat2[i]!='NA'),"Ind"])
  N_fams = rbind(N_fams, c(trait,nfams,ninds))
  names(N_fams) = c("Trait","N families","N individuals")
  }

N_fams

# Generate BLUPs for each trait and merge with ancestry values
# Not accounting for Pop rand effect in BLUP estimation
gen2 = gen[,1:3]
names(gen2) = c("Pop","Ind","K2_ancestry")
gen2$Ind = as.factor(gen2$Ind)

for(i in 3:8) {
  tmp = dat2[which(dat2[i]!='NA'),]
  trait = scale(tmp[i])
  Ind = as.factor(tmp$Ind)
  bench = as.factor(tmp$bench)
  tmp.mod = lmer(trait ~ bench + (1|Ind))
  tmp2 = coef(tmp.mod)$Ind
  names(tmp2)[1] = names(tmp[i])
  tmp2$Ind = row.names(tmp2)
  gen2 = merge(tmp2[,c(1,4)], gen2)
}

# Trait models on BLUPs ((1|Ind) model) 
# Can either include/exclude random Pop effect in 2nd round of LMMs

# without pop effect
Model.results = list()

for(i in 1:6) {
  trait = as.numeric(gen2[,i+1])
  K2_ancestry = gen2$K2_ancestry
  Pop = gen2$Pop # optional, depending on desired model
  tmp.lm2 = lm(trait ~ K2_ancestry + I(K2_ancestry^2))
  tmp.lm = lm(trait ~ K2_ancestry) 
  tmp.LRT = anova(tmp.lm, tmp.lm2)
  Model.results[[i]] = list()
  Model.results[[i]][[1]] = tmp.lm
  Model.results[[i]][[2]] = tmp.lm2
}

# Make table of model estimates across traits
tab_model(Model.results[[6]][[2]],
          Model.results[[5]][[2]],
          Model.results[[4]][[2]],
          Model.results[[3]][[2]],
          Model.results[[2]][[2]],
          Model.results[[1]][[2]],
          dv.labels=c("Height","Stem Width","Biomass","SLA","1 / DFF","Tot. Flower Heads"),
          string.se = "SE",
          show.ci=F, show.se=T)

# with pop random effect in model
LRT.resultsp = data.frame()
Model.resultsp = list()

for(i in 1:6) {
  trait = as.numeric(gen2[,i+1])
  K2_ancestry = gen2$K2_ancestry
  Pop = gen2$Pop # optional, depending on desired model
  tmp.lm2 = lmer(trait ~ K2_ancestry + I(K2_ancestry^2) + (1|Pop)) # (1|Pop) optional, depending on desired model
  tmp.lm = lmer(trait ~ K2_ancestry + (1|Pop)) # (1|Pop) optional, depending on desired model
  tmp.LRT = anova(tmp.lm, tmp.lm2)
  LRT.resultsp = rbind(LRT.results, c(names(gen2[i+1]), tmp.LRT$Chisq[2], tmp.LRT$`Pr(>Chisq)`[2])) # last 2 values will change to F and P, respectively if using 'lm' models without Pop effect
  names(LRT.resultsp) = c("trait", "ChiSq", "P") # only if using 'lmer' models with Pop effect
  Model.resultsp[[i]] = list()
  Model.resultsp[[i]][[1]] = tmp.lm
  Model.resultsp[[i]][[2]] = tmp.lm2
}

LRT.resultsp

# Make table of model estimates across traits
tab_model(Model.resultsp[[6]][[2]],
          Model.resultsp[[5]][[2]],
          Model.resultsp[[4]][[2]],
          Model.resultsp[[3]][[2]],
          Model.resultsp[[2]][[2]],
          Model.resultsp[[1]][[1]],
          dv.labels=c("Height","Stem Width","Biomass","SLA","1 / DFF","Tot. Flower Heads"),
          string.se = "SE",
          show.ci=F, show.se=T)

# Make trait~admixture plots without Pop effect
dat_hyb = gen2[,c(1,7,6,5,4,3,2,9)]
dat_hyb$class = as.factor(ifelse(dat_hyb$K2_ancestry>0.8, "nigra", "jacea"))
dat_hyb$class = ifelse((dat_hyb$K2_ancestry*(1-dat_hyb$K2_ancestry)>0.16), "hybrid", dat_hyb$class)

trait_plot_list = list()

for(i in 1:6) {
  dat_hyb$trait = scale(dat_hyb[i+1])

  p <- ggplot(dat_hyb, aes(x=K2_ancestry, y=trait, size=0.5, colour=K2_ancestry)) + geom_point(size=0.75) + 
    scale_colour_gradient(low="blue", high="red") +
    ggtitle(names(dat_hyb[i+1])) + xlab("Ancestry score (K=2)") + ylab("standardized trait") +
    theme(legend.position = "none") + 
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, colour="black") +
    theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))
  
  trait_plot_list[[i]] = p
}

traits <- marrangeGrob(trait_plot_list, nrow=1, ncol=6, top=NULL)
#ggsave("figs/trait_plots.pdf", traits)

# Make trait~admixture plots with Pop effect
gen2p = gen[,1:3]
names(gen2p) = c("Pop","Ind","K2_ancestry")
gen2p$Ind = as.factor(gen2p$Ind)
gen2p$Pop = as.factor(gen2p$Pop)

# Do so by accounting for Pop rand effect in BLUP estimation of (1|Ind)
for(i in 3:8) {
  tmp = dat2[which(dat2[i]!='NA'),]
  trait = scale(tmp[i])
  Ind = as.factor(tmp$Ind)
  Pop = as.factor(tmp$Pop)
  bench = as.factor(tmp$bench)
  tmp.mod = lmer(trait ~ bench + (1|Ind) + (1|Pop))
  tmp2 = coef(tmp.mod)$Ind
  names(tmp2)[1] = names(tmp[i])
  tmp2$Ind = row.names(tmp2)
  gen2p = merge(tmp2[,c(1,4)], gen2p)
}

dat_hybp = gen2p[,c(1,7,6,5,4,3,2,9)]
dat_hybp$class = as.factor(ifelse(dat_hybp$K2_ancestry>0.8, "nigra", "jacea"))
dat_hybp$class = ifelse((dat_hybp$K2_ancestry*(1-dat_hybp$K2_ancestry)>0.16), "hybrid", dat_hybp$class)

trait_plot_listp = list()

for(i in 1:6) {
  dat_hybp$trait = scale(dat_hybp[i+1])
  
  p <- ggplot(dat_hybp, aes(x=K2_ancestry, y=trait, size=0.5, colour=K2_ancestry)) + geom_point(size=0.75) + 
    scale_colour_gradient(low="blue", high="red") +
    ggtitle(names(dat_hybp[i+1])) + xlab("Ancestry score (K=2)") + ylab("standardized trait") +
    theme(legend.position = "none") + 
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1, colour="black") +
    theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10))
  
  trait_plot_listp[[i]] = p
}

traitsp <- marrangeGrob(trait_plot_listp, nrow=1, ncol=6, top=NULL)
# ggsave("figs/trait_plots_pop.pdf", width=14,height=4, traitsp)

# make histograms of trait values by ancestry class
# without pop effects
hist_plot_list = list()

for(j in 1:6){
  dat_hyb$trait = scale(dat_hyb[j+1])
  p <- ggplot(dat_hyb, aes(x = trait)) +
    geom_histogram(aes(color = class, fill = class), position = "identity", bins = 15, alpha = 0.4) +
    scale_color_manual(values = c("blue","red","purple")) +
    scale_fill_manual(values = c("blue","red","purple")) +
    xlab("standardized trait") +
    theme(legend.position = "none")
  
  hist_plot_list[[j]] = p
}

hist <- marrangeGrob(hist_plot_list, nrow=1, ncol=6, top=NULL)

# accounting for pop effects 
hist_plot_listp = list()

for(j in 1:6){
  dat_hybp$trait = scale(dat_hybp[j+1])
  p <- ggplot(dat_hybp, aes(x = trait)) +
      geom_histogram(aes(color = class, fill = class), position = "identity", bins = 15, alpha = 0.4) +
      scale_color_manual(values = c("blue","red","purple")) +
      scale_fill_manual(values = c("blue","red","purple")) +
      xlab("standardized trait") +
      theme(legend.position = "none")
  
  hist_plot_listp[[j]] = p
}

histp <- marrangeGrob(hist_plot_listp, nrow=1, ncol=6, top=NULL)
# ggsave("figs/hist_plots_pop.pdf", hist)

# without pop effects
multpanel <- marrangeGrob(c(traits,hist), nrow=2,ncol=1, top=NULL, padding = unit(10, "line"))
ggsave("figs/Figure1.pdf", multpanel, width=14,height=8)

# With pop effects
multpanelp <- marrangeGrob(c(traitsp,histp), nrow=2,ncol=1, top=NULL, padding = unit(10, "line"))
ggsave("figs/FigureS1.pdf", multpanelp, width=14,height=8)

