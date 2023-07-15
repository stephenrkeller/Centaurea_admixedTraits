# do hybrids (ind with high diffK1) have higher heritability (i.e. higher evolutionary potential)
# June 20 2019
# ZMP

# load packages
library(tidyverse)

# load data
mid <- read.csv(file = 'midBLUPs.csv', header = T, sep = ',')
summary(mid)
# diffK1: median = 0.3331, mean = 0.4199
# make columns for midparent values
mid$midHt <- (mid$sireHt+mid$damHt)/2
mid$midWt <- (mid$sireWt+mid$damWt)/2
mid$midSW <- (mid$sireSW+mid$damSW)/2
mid$midNCap <- (mid$sireNCap+mid$damNCap)/2
mid$midHtNFl <- (mid$sireNFl+mid$damNFl)/2
mid$middff <- (mid$siredff+mid$damdff)/2
mid$midsla <- (mid$siresla+mid$damsla)/2

## quadratic variables
mid$midHtSq <- mid$midHt^2
mid$midWtSq <- mid$midWt^2
mid$midSWSq <- mid$midSW^2
mid$midNCapSq <- mid$midNCap^2
mid$midNFlSq <- mid$midHtNFl^2
mid$middffSq <- mid$middff^2
mid$midslaSq <- mid$midsla^2

lowH <- filter(mid, diffK1mean < 0.3331)
highH <- filter(mid, diffK1mean > 0.3331)

###############################
# Ht low
htlo <- lm(data =lowH,ht ~ midHt)
# htlo2 <- lm(data = lowH, ht ~ midHt + midHtSq)
summary(htlo)
# summary(htlo2)
# anova(htlo,htlo2)

ggplot(data = lowH, aes(x = midHt, y = ht)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-1,y=20,label = "Low hybridity,
            p-value = 0.01317, Rsq = 0.3013")

ggsave(filename = "lowHHt.png",plot = last_plot())

# Ht High
hthi <- lm(data =highH,ht ~ midHt)
# hthi2 <- lm(data = highH, ht ~ midHt + midHtSq)
summary(hthi) # n.s.
# summary(hthi2)
# anova(hthi,hthi2)

ggplot(data = highH, aes(x = midHt, y = ht)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-7,y=10,label = "High hybridity,
            n.s.")

ggsave(filename = "highHHt.png",plot = last_plot())

# Wt low
wtlo <- lm(data =lowH,wt ~ midWt)
# wtlo2 <- lm(data = lowH, wt ~ midWt + midWtSq)
summary(wtlo) # n.s.
# summary(wtlo2)
# anova(wtlo,wtlo2) ### neither sig

ggplot(data = lowH, aes(x = midWt, y = wt)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-1.0,y=3,label = "Low hybridity,
            n.s.")

ggsave(filename = "lowHWt.png",plot = last_plot())

# Wt high
wthi <- lm(data =highH,wt ~ midWt)
# wthi2 <- lm(data = highH, wt ~ midWt + midWtSq)
summary(wthi) # n.s.
# summary(wthi2)
# anova(wthi,wthi2) ### neither sig

ggplot(data = highH, aes(x = midWt, y = wt)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-1.0,y=2.5,label = "High hybridity,
            n.s.")

ggsave(filename = "highHWt.png",plot = last_plot())

# SW lo
swlo <- lm(data =lowH,stemwid ~ midSW)
# swlo2 <- lm(data = lowH, stemwid ~ midSW + midSWSq)
summary(swlo)
# summary(swlo2)
# anova(swlo,swlo2)

ggplot(data = lowH, aes(x = midSW, y = stemwid)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-.3,y=0.4,label = "Low hybridity,
            p-value = 0.01599, Rsq = 0.2847")

ggsave(filename = "lowHSW.png",plot = last_plot())

# sw high
swhi <- lm(data =highH,stemwid ~ midSW)
# swhi2 <- lm(data = highH, stemwid ~ midSW + midSWSq)
summary(swhi)
# summary(swhi2)
# anova(swhi,swhi2)

ggplot(data = highH, aes(x = midSW, y = stemwid)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-.3,y=0.35,label = "High hybridity,
            n.s.")

ggsave(filename = "highHSW.png",plot = last_plot())

# NCap low
nclo <- lm(data =lowH,NCap ~ midNCap)
nclo2 <- lm(data = lowH, NCap ~ midNCap + midNCapSq)
summary(nclo)
summary(nclo2)
anova(nclo,nclo2) #### neither sig

ggplot(data = lowH, aes(x = midNCap, y = NCap)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-2,y=9,label = "Low hybridity,
            n.s.")

ggsave(filename = "lowHNCap.png",plot = last_plot())

# ncap high
nchi <- lm(data =highH,NCap ~ midNCap)
nchi2 <- lm(data = highH, NCap ~ midNCap + midNCapSq)
summary(nchi)
summary(nchi2)
anova(nchi,nchi2) #### neither sig

ggplot(data = highH, aes(x = midNCap, y = NCap)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-3,y=7,label = "High hybridity,
            n.s.")

ggsave(filename = "highHNCap.png",plot = last_plot())

# NFl low
nflo <- lm(data =lowH,NFl ~ midHtNFl)
nflo2 <- lm(data = lowH, NFl ~ midHtNFl + midNFlSq)
summary(nflo)
summary(nflo2)
anova(nflo,nflo2)

ggplot(data = lowH, aes(x = midHtNFl, y = NFl)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-1.2,y=10,label = "Low hybridity,
            p-value = 0.078 (n.s.), Rsq = 0.1391")

ggsave(filename = "lowHNFl.png",plot = last_plot())

#nFl high
nfhi <- lm(data =highH,NFl ~ midHtNFl)
# nfhi2 <- lm(data = highH, NFl ~ midHtNFl + midNFlSq)
summary(nfhi)
# summary(nfhi2)
# anova(nfhi,nfhi2)

ggplot(data = highH, aes(x = midHtNFl, y = NFl)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-1.2,y=6,label = "High hybridity,
            p-value = 0.01074, Rsq = 0.3838")

ggsave(filename = "highHNFl.png",plot = last_plot())

# dff low
dfflo <- lm(data =lowH,dff ~ middff)
# dfflo2 <- lm(data = lowH, dff ~ middff + middffSq)
summary(dfflo)
# summary(dfflo2)
# anova(dfflo,dfflo2)

ggplot(data = lowH, aes(x = middff, y = dff)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-1,y=5,label = "Low hybridity,
            p-value = 0.007884, Rsq = 0.3869")

ggsave(filename = "lowHDFF.png",plot = last_plot())

# dff high
dffhi <- lm(data =highH,dff ~ middff)
# dffhi2 <- lm(data = highH, dff ~ middff + middffSq)
summary(dffhi)
# summary(dffhi2)
# anova(dffhi,dffhi2)

ggplot(data = highH, aes(x = middff, y = dff)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-2.5,y=5,label = "High hybridity,
            p-value = 0.02789, Rsq = 0.4074")

ggsave(filename = "highHDFF.png",plot = last_plot())

# sla low
# slalo <- lm(data =lowH,sla ~ midsla)
slalo2 <- lm(data = lowH, sla ~ midsla + midslaSq)
# summary(slalo)
summary(slalo2)
anova(slalo,slalo2) ###### sig

ggplot(data = lowH, aes(x = midsla, y = sla)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  geom_text(x=0,y=3e-14,label = "Low hybridity,
            p-value = 0.0929 (n.s.), Rsq = 0.2328")

ggsave(filename = "lowHSLA.png",plot = last_plot())

# sla high
slahi <- lm(data =highH,sla ~ midsla)
slahi2 <- lm(data = highH, sla ~ midsla + midslaSq)
summary(slahi)
summary(slahi2)
anova(slahi,slahi2) ###### not sig

ggplot(data = highH, aes(x = midsla, y = sla)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=0,y=2e-14,label = "High hybridity,
            n.s.")

ggsave(filename = "highHSLA.png",plot = last_plot())

# totFl lo
lowH$TotFl <- lowH$NCap + lowH$NFl
lowH$midTotFl <- lowH$midNCap + lowH$midHtNFl
lowH$midTotFlSq <- lowH$midTotFl^2

tflo <- lm(data =lowH,TotFl ~ midTotFl)
tflo2 <- lm(data = lowH, TotFl ~ midTotFl + midTotFlSq)
summary(tflo)
summary(tflo2)
anova(tflo,tflo2) ## n.s.

ggplot(data = lowH, aes(x = midTotFl, y = TotFl)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-2.5,y=9,label = "Low Hybridity,
            n.s.")

ggsave(filename = "lowHTotFl.png",plot = last_plot())

### totFl high
highH$TotFl <- highH$NCap + highH$NFl
highH$midTotFl <- highH$midNCap + highH$midHtNFl
highH$midTotFlSq <- highH$midTotFl^2

tfhi <- lm(data =highH,TotFl ~ midTotFl)
tfhi2 <- lm(data = highH, TotFl ~ midTotFl + midTotFlSq)
summary(tfhi)
summary(tfhi2)
anova(tfhi,tfhi2) ## n.s.

ggplot(data = highH, aes(x = midTotFl, y = TotFl)) +
  geom_point() +
  stat_smooth(method = "lm", size = 1) +
  geom_text(x=-6,y=5,label = "High Hybridity,
            p-value = 0.899 (n.s.), Rsq = 0.156")

ggsave(filename = "highHTotFl.png",plot = last_plot())

#############

# trait  # lowH P  # highH P  # lowH R2  # highH R2
# height # 0.013   # n.s.     # 0.3013   # n.s.
# weight # n.s.    # n.s.     # n.s.     # n.s.
# StWid  # 0.016   # n.s.     # 0.2847   # n.s.
# NCap   # n.s.    # n.s.     # n.s.     # n.s.
# NFl    # n.s.    # 0.011    # n.s.     # 0.3838
# dFF    # 0.008   # 0.0279   # 0.3869   # 0.4074
# SLA    # n.s.    # n.s.     # n.s.     # n.s.
# totFl  # n.s.    # n.s.     # n.s.     # n.s.


