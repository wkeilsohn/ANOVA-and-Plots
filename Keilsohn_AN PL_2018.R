### William Keilsohn
### Anova

### Load Packages
library(dplyr)
library(nparcomp)
library(agricolae)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggsignif)


### Anova on three experimental groups
## Sample size is assumed to be equal across all groups
PA1<-aov(Data ~ Group, data = Prob.1)
summary(PA1)
#Test for normality
shapiro.test(residuals(PA1)) ### Not normally distributed.
# Test for homogeniety
Prob.1.2<-stack(Prob.1)
bartlett.test(Prob.1.2) ### Can't do this b/c it's not normal


#Transformations
Prob.1.3<-mutate(Prob.1, LogData = log10(Prob.1$Data))
PA2<-aov(LogData ~ Group, data = Prob.1.3)
shapiro.test(residuals(PA2))
Prob.1.3<-c(Prob.1.3$Group, Prob.1.3$LogData)
Prob.1.3.2<-stack(Prob.1.3)
bartlett.test(Prob.1.3)

Prob.1.4<-mutate(Prob.1, ExpData = log(Prob.1$Data))
PA3<-aov(ExpData ~ Group, data = Prob.1.4)
shapiro.test(residuals(PA3))
Prob.1.4.2<-data.frame(Prob.1.4$Group, Prob.1.4$ExpData)
bartlett.test(list(Prob.1.4$ExpData, Prob.1.4$Group))
summary(PA3)
TukeyHSD(PA3)

### Looking for added variation
## Anova across four groups
## Sample size is equal and chosen at random from each group.
P11A.1<-aov(CellAct ~ Trtmnt, data = P11.3)
shapiro.test(residuals(P11A.1))
bartlett.test(P11.3$CellAct, P11.3$Trtmnt)
summary(P11A.1)
TukeyHSD(P11A.1)

# Transformations
P11.3.2<-mutate(P11.3, LogData = log(P11.3$CellAct))
P11A.2<-aov(LogData ~ Trtmnt, data = P11.3.2)
shapiro.test(residuals(P11A.2))
bartlett.test(P11.3.2$LogData, P11.3.2$Trtmnt)
K1<-kruskal(P11.3$CellAct, P11.3$Trtmnt) ### use the kruskal command, instead of kruskal.test
summary(mctp(CellAct ~ Trtmnt, data = P11.3))


### Testing for the affect of a chemical
## Anova is used
## Equal sample size across three groups.
P11A.3<-aov(mm ~ Trtmnt, data = P11.10)
shapiro.test(residuals(P11A.3))
bartlett.test(P11.10)

# Transformations
P11.10.2<-mutate(P11.10, LogData = log(P11.10$mm))
P11A3.2<-aov(LogData ~ Trtmnt, data = P11.10.2)
shapiro.test(residuals(P11A3.2))
bartlett.test(list(P11.10.2$LogData, P11.10.2$Trtmnt))

K3<-kruskal(P11.10$mm, P11.10$Trtmnt)

### Hormones are examined in birds
## Sample size is equal amoung groups
## Individuals are chosen at random
P11A.4<-aov(TestL ~ Strain, data = P11.12)
shapiro.test(residuals(P11A.4))

## Transformations
P11.12.2<-mutate(P11.12, LogData = log(P11.12$TestL))
P11A.4.2<-aov(LogData ~ Strain, data = P11.12.2)
shapiro.test(residuals(P11A.4.2))
K2<-kruskal(P11.12$TestL, P11.12$Strain) ### This is the important line*****
summary(mctp(TestL~Strain, data = P11.12))
kruskal.test(TestL~Strain, data = P11.12)


### Urchin hormone levels are compared
## Anova is used
## Sample size amoung groups is assumed to be equal.
UA1<-aov(E2 ~ Stage, data = Parac)
shapiro.test(residuals(UA1))
bartlett.test(Parac$E2, Parac$Stage)

# Transform
Parac2<-mutate(Parac, LogData = log(Parac$E2))
UA2<-aov(LogData ~ Stage, data = Parac2)
shapiro.test(residuals(UA2))
K4<-kruskal(Parac$E2, Parac$Stage)
K4


### Series of plots
urchin_comparisons<-list(c("Mature","Growing"),c("Mature","Recovery"),
                         c("Growing", "Recovery"))
K4.2<-tapply(Parac$E2, Parac$Stage, median, na.rm = TRUE)
graph1<-ggplot(data = Parac, aes(x = Stage, y = E2))+
  geom_boxplot(fill = c("grey69","grey50","grey40"))+
  #stat_compare_means(comparisons = urchin_comparisons, position = 850)+
  stat_compare_means(label.y = 700)+
  #geom_signif(comparisons = list(c("Mature","Growing")), map_signif_level = TRUE)+
  #geom_signif(comparisons = list(c("Mature","Recovery")), map_signif_level = TRUE)+
  #geom_signif(comparisons = list(c("Recovery","Growing")), map_signif_level = TRUE)+
  plotPairwiseTests(K4$comparison, K4.5, alpha = 0.05)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(-2, "mm"), axis.title.y = element_text(size = 15),
        axis.title = element_text(size = 15), axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 12))
graph1<-graph1 + labs(x = "Urchin Life Stage", y = "Estradiol Level per Gonad Wet Mass(pg/g)")
graph1 

### New idea for a plot
graph2<-ggplot(data = Parac, aes(x = Stage, y = E2))+
  geom_tufteboxplot()
graph2 ### This is just a shit show. 