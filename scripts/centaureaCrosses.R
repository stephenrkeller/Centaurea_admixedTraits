# Centaurea crosses 2017
# March 27 2019
# ZMP

library(tidyverse)

crosses <- read.csv(file = "Centaurea_Crosses_2017ZMP.csv", header = T, sep = ',')
summary(crosses)
head(crosses)

#crosses
table(crosses$dam, crosses$sire)

#crosses attempted (# attempted), excluding NA parents: 
library(BiocManager)
# BiocManager::install("tximport")
library(tximportData)
library(tximport)
library(DESeq2)


# JF14.1xCO2.1 x JF14.1xCO2.1 (1)
# JF14.1xFC21.1 x JF14.1xFC21.1 (71)
# JF28.1xRM54.2 x JF28.1xRM54.2 (10)
# JF28.2xRM54.2 x JF28.2xRM54.2 (1)
# RM53.1xCO2.1 x RM53.1(B)xCO2.1 (73)
# SP20.2xCO2.1 x SP20.1xCO2.1 (35)

glimpse(crosses)
# crosses <- crosses[1:12]
# glimpse(crosses)
# crosses <- drop_na(crosses)
# glimpse(crosses)
# change grandparents from damFam1 -> matGDam, etc.

table(crosses$matGDam, crosses$matGSire)

##### noDevSeeds
# 0 inflated distribution
hist <- ggplot(crosses, aes(noDevSeeds)) +
  geom_histogram(aes(fill = dam),
                 binwidth = 1) +
  labs(title = "Number of Developed Seeds")
hist 
ggsave(file = "hist_devseeds.png", plot =hist)
### By matGDam
######### make sure to scale per cross individual
# total seeds
ggplot(crosses, aes(matGDam, totalSeeds)) +
  geom_boxplot(aes(fill=matGSire))
ggsave(file='box_totalseed_DambySire.png',plot = last_plot())
# dev seeds
ggplot(crosses, aes(matGDam, noDevSeeds)) +
  geom_boxplot(aes(fill=matGSire))

# small seeds
ggplot(crosses, aes(matGDam, noSmallSeeds)) +
  geom_boxplot(aes(fill=matGSire))


# Take out zeros to see if Poisson (tail goes below 0) or Gaussian
crosses_nozeros <- crosses %>%
  filter(noDevSeeds > 0)
glimpse(crosses_nozeros)

crosses_nozeros %>%
  group_by(matGDam,matGSire) %>%
  summarise(n())

ggplot(crosses_nozeros, aes(matGDam, totalSeeds)) +
  geom_boxplot(aes(fill=matGSire))
ggsave("seeds_gdambydsire_nozeros.png", plot = last_plot())

histNoZero <- ggplot(crosses_nozeros, aes(noDevSeeds)) +
  geom_histogram(aes(fill = dam),
                 binwidth = 1) +
  labs(title = "Number of Developed Seeds", subtitle = "excluding crosses with 0 developed seeds")

histNoZero +
  geom_vline(aes(xintercept = mean(noDevSeeds)))

mean(crosses$noDevSeeds) # 11.51309
var(crosses$noDevSeeds) # 366.9669
#negative binomial? do both, depends on the data


# boxplot of dam means (assuming sires are the same as dams)
# with zeros
box <- ggplot(crosses, aes(dam, noDevSeeds)) +
  geom_boxplot()
box

# without zeros
box2 <- ggplot(crosses_nozeros, aes(matGDam, noDevSeeds)) +
  geom_boxplot()
box2

#### totalSeeds
# 0 inflated distribution
histTotal <- ggplot(crosses, aes(totalSeeds)) +
  geom_histogram(aes(fill = dam),
                 binwidth = 1) +
  labs(title = "Number of Total Seeds")
histTotal

# Take out zeros to see if Poisson (tail goes below 0) or Gaussian
crosses_nozerosTotal <- crosses %>%
  filter(totalSeeds > 0)
crosses_nozerosTotal

histNoZeroTotal <- ggplot(crosses_nozerosTotal, aes(totalSeeds)) +
  geom_histogram(aes(fill = dam),
                 binwidth = 1) +
  labs(title = "Number of Total Seeds", subtitle = "excluding crosses with 0 small seeds")

histNoZeroTotal +
  geom_vline(aes(xintercept = mean(totalSeeds)))

mean(crosses_nozerosTotal$totalSeeds) # 47.70093

# boxplot of dam means (assuming sires are the same as dams)
# with zeros
box <- ggplot(crosses, aes(dam, totalSeeds)) +
  geom_boxplot()
box

# without zeros
box2 <- ggplot(crosses_nozeros, aes(dam, totalSeeds)) +
  geom_boxplot()
box2

####### What is the average number of total and developed seeds per cross?
glimpse(crosses)
damMeans<- crosses %>%
  group_by(dam) %>%
  summarise(n(),mean(noDevSeeds), mean(totalSeeds),
            mean(noDevSeeds)/n(), mean(totalSeeds)/n(),
            mean(noDevSeeds)/mean(totalSeeds))
colnames(damMeans) <- c('dam','n','meanDev','meanTotal','propDev', 'propTotal', 'percentDev')
damMeans

box2 <- ggplot(damMeans, aes(dam, percentDev)) +
  geom_boxplot()
box2
crosses
#### What is percent developed for each of the crosses (not just the dams)
crosses <- crosses %>%
  mutate(percentDev = noDevSeeds/totalSeeds)
head(crosses)

# with zeros
percentDevHist <- ggplot(crosses, aes(percentDev)) +
  geom_histogram(aes(fill = dam),
                 binwidth = .1) +
  labs(title = "Percent of Developed Seeds/Total Seeds (with 0s)")

percentDevHist +
  geom_vline(aes(xintercept = mean(percentDev)))
ggsave(file = "percdevoftotalseeds.png", plot = last_plot())
# binomial model with 2 groups?

# without zeros
percentDevHist_noZeros <- ggplot(crosses_nozeros, aes(percentDev)) +
  geom_histogram(aes(fill = dam),
                 binwidth = .1) +
  labs(title = "Number of developed/total seeds", subtitle = "excluding 0s")

# percentDevHist_noZeros
  # geom_vline(aes(xintercept = mean(damMeans$percentDev)))
# anova?

############# zero inflated poisson
# install.packages("pscl")
library(pscl)
summary(m1 <- zeroinfl(totalSeeds ~ damFam1|1, dist = "negbin",data=crosses, link = "logit"))

table(crosses$damFam1, crosses$damFam2)
table(crosses$sireFam1, crosses$sireFam2)



# Keep it simple
### noSeeds ~ crosstype (factor with n levels)
# cross types: the grandparents 
# 0 inflated neg binomial




