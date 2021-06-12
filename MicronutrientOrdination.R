

#----------------------------------------------------------------------------------------------------------------------------------------
#project:Shell micronutrient project 
#coder:G.W. Hopper
#----------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())
######################## 
######Step 1.##########
#######################

#Set working directory and import micronutrient data
setwd("~/Shell/Micronutrients")
micronutr<-read.csv(file="ShellMicronutrientsMass.csv", header = TRUE, sep = ",", row.names = c())
head(micronutr)
micronutr<-subset(micronutr, ID !="Fcer6" & ID !="Ouni1" & ID != "Ouni6" & ID != "Tver3" & ID != "Tver4")#no data for these individuals
#create a factor for species names, use attach to call by name

micronutr$Spp<-factor(micronutr$Spp, levels=c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor",  order=TRUE))
micronutr$Strategy<-factor(micronutr$Strategy, levels=c("Equilibrium","Periodic",  order=TRUE))
micronutr$Shell<-factor(micronutr$Shell, levels=c("UD","General", "Unsculp", order=TRUE))

attach(micronutr)

#load required packages.
library(vegan)
library(cowplot)
library(ggpubr)
library(ggfortify)
library(factoextra)
library(dplyr)
library(MVN)#for multivariate normatlity tests
library(car)
library(tidyverse)
library(dunn.test)
#============================================
#is length(age) related to nutrients?
#create a data frame of micronutrients
Micros<-(as.data.frame(micronutr[10:20]))
attach(Micros)

#Run principle components analysis
micro.pca<-prcomp(na.omit(log(Micros)), scale =TRUE)
ordiplot(micro.pca)

#check eigen values with scree plot
fviz_eig(micro.pca)
eig.val <- get_eigenvalue(micro.pca)
eig.val

screeplot(micro.pca, bstick = TRUE, type = c("barplot", "lines"),
          npcs = min(11, length(micro.pca$sdev)),
          ptype = "o", bst.col = "red", bst.lty = "solid",
          xlab = "Component", ylab = "Inertia",
          main = deparse(substitute(micro.pca)))

#Extract loadings for supplemental table and determine order for ilrs
pca.loadings <- micro.pca$rotation
pca.loadings
#extract PCA scores for axis 1, 2, 3
micronutr$pca.scores <- micro.pca$x[,1:3]


#need to do this again because more column headers have been added at this point.
attach(micronutr)
head(micronutr)

#install national parks palettes to make graph pretty as NP
#devtools::install_github("katiejolly/nationalparkcolors")
library(nationalparkcolors)
RWpal<-park_palette("Redwoods")
Evpal<-park_palette("Everglades",4)

pal<-c(RWpal,Evpal) #colorpaletts

#simple ordinatioN
Microplot<- autoplot(micro.pca, x = 1, y = 2, data = micronutr, fill = 'Spp', shape = 'Spp', size = 5,
         loadings = TRUE, loadings.colour = 'grey45', loadings.label = TRUE, loadings.label.size = 7, 
         loadings.label.colour = 'black',loadings.label.hjust = 1.2, loadings.labels.vjust=1)+
  stat_conf_ellipse(aes(color = Spp, fill = Spp, linetype = Spp), alpha = 0.05, size =1.5, geom = "polygon")+
  geom_hline(aes(yintercept=0), linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0), linetype="dashed", size=1)+
  scale_linetype_manual(values =c(1,1,1,1,1,3,3))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  scale_shape_manual(values = c(21,22,23,24,25,21,22))+
  theme(legend.position = c(.01,.85), 
        legend.title = element_blank(), 
        legend.text = element_text(face = "italic"), 
        legend.key.size =  unit(.25, 'cm'),
        legend.background = element_rect(linetype = 2, color="black", fill = "grey97"),
        panel.background = element_rect(fill="grey90"));Microplot
#adjustments to loading vectors
Microplot$layers[[2]]$aes_params$size <- 1
Microplot$layers[[2]]$geom_params$arrow$length <- unit(0, units = "points")
Microplot



#check data for assumptions of multivariation normality and heterogeneity
#prior to MANOVA
#=============================================================================
#check Mardia's test and plot chi- squared quantiles vs observed
mardia_result <- mvn(data = log10(Micros), mvnTest = "mardia", multivariatePlot = "qq")
mardia_result$multivariateNormality

#PERMANOVA esting for species differences
set.seed(1)
as.matrix(Micros)
adonis(Micros~Spp*Length, data=micronutr, permutations = 999, method = "bray",
       strata = NULL, contr.unordered= "contr.sum",
       contr.ordered = "contr.poly")
#BC distance matrix prep step to test for dispersion among species
bc_micros<-vegdist(Micros, method = "bray")
#test for dispersion among groups
disper_micros_spp<-betadisper(bc_micros, Spp, type = "centroid")
anova(disper_micros_spp)

#PERMANOVA morphology differences
set.seed(1)
adonis(Micros~Shell*Length, data=micronutr, permutations = 999, method = "bray",
       strata = NULL, contr.unordered= "contr.sum",
       contr.ordered = "contr.poly")
#BC distance matrix prep step to test for dispersion among morphologies
bc_micros<-vegdist(Micros, method = "bray")
#test for dispersion among groups
disper_micros_shell<-betadisper(bc_micros, Shell, type = "centroid")
anova(disper_micros_shell)

#PERMANOVA strategies differences
set.seed(1)
adonis(Micros~Strategy*Length, data=micronutr, permutations = 999, method = "bray",
       strata = NULL, contr.unordered= "contr.sum",
       contr.ordered = "contr.poly")

#BC distance matrix prep step to test for dispersion among stratetgies
bc_micros<-vegdist(Micros, method = "bray")
#test for dispersion among groups
disper_micros_strategy<-betadisper(bc_micros, Strategy, type = "centroid")
anova(disper_micros_strategy)

#run a multivariate analysis of variance for 10 elements testing for species effects
elem.man <- manova(log(cbind(C,N,P, K, B, Zn, Fe, Cu, Ca, Mn, Mg)) ~ Spp, data = micronutr)
summary(elem.man, test = "Wilks", tol = 0)
summary.aov(elem.man)

#run a multivariate analysis of variance for 10 elements testing for morphological effects
elem.morph<- manova(log(cbind(C,N,P, K, B, Zn, Fe, Cu, Ca, Mn, Mg)) ~ Shell, data = micronutr)
summary(elem.morph, test = "Wilks", tol = 0)
summary.aov(elem.morph)
#run a multivariate analysis of variance for 10 elements testing for life history effects
elem.LH<- manova(log(cbind(C,N,P, K, B, Zn, Fe, Cu, Ca, Mn, Mg)) ~ Strategy, data = micronutr)
summary(elem.LH, test = "Wilks", tol = 0)
summary.aov(elem.LH)
library(emmeans)
attach(micronutr)
#perform univariate tests of normality and hetergeneity; run anova, and pairwises test for species
shapiro.test(log(C))#need to exchange elements for each test
leveneTest(log(C)~Spp)
lmC<-lm(log(C)~Spp, data= micronutr)
emmeans(lmC, list(pairwise~Spp), adjust = "Tukey")
#non-parametric; Kruskal-Wallis test among species
dunn.test(C, Spp, method="bonferroni", kw = TRUE)
dunn.test(Cu, Spp, method="bonferroni", kw = TRUE)
dunn.test(N, Spp, method="bonferroni", kw = TRUE)
dunn.test(P, Spp, method="bonferroni", kw = TRUE)
#non-parametric; Kruskal-Wallis test among morphology
dunn.test(B, Shell, method="bonferroni", kw = TRUE)
dunn.test(C, Shell, method="bonferroni", kw = TRUE)
dunn.test(Cu, Shell, method="bonferroni", kw = TRUE)
dunn.test(N, Shell, method="bonferroni", kw = TRUE)
dunn.test(P, Shell, method="bonferroni", kw = TRUE)
#non-parametric; Kruskal-Wallis test among LH strategies
dunn.test(B, Strategy, method="bonferroni", kw = TRUE)
dunn.test(C, Strategy, method="bonferroni", kw = TRUE)
dunn.test(Cu, Strategy, method="bonferroni", kw = TRUE)
dunn.test(N, Strategy, method="bonferroni", kw = TRUE)
dunn.test(P, Strategy, method="bonferroni", kw = TRUE)
#

library(tidyr)
Bulk_long <- gather(micronutr, Element, Mass, C:K, factor_key=TRUE)
Bulk_long
#create box plot for all elements
elem.bulk<-ggplot()+
  geom_boxplot(data = Bulk_long, aes(x=Spp, y =Mass, fill = Tribe), width = .5, size = 1)+
  facet_wrap(~Element, ncol =1, nrow = 5, strip.position = "right", scales = "free_y")+
  ylab(expression(Element~(mg~g^{"-1"})))+
  scale_fill_manual(values =Evpal)+
  ggtitle("Bulk elements")+
  theme_bw()+
  theme(legend.position = "bottom", axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 11, margin=margin(t=-10)),
        axis.text.x = element_text(face = "italic", size = 12, angle = 45, hjust =.75, vjust = .75),
        panel.grid = element_line(color = "grey81"), panel.spacing = unit(.75, "lines"))
elem.bulk1<-elem.bulk+
  geom_point(data = Bulk_long %>% 
               group_by(Element) %>% #group by facet variable
               summarise(y.min = pretty(Mass)[1],
                         y.max = pretty(Mass)[length(pretty(Mass))]) %>%
               tidyr::gather(key, value, -Element), 
             aes(x = 1, y = value),
             inherit.aes = FALSE, alpha = 0) +
  
  # Turn off automatical scale expansion, & manually set scale breaks
  # as an evenly spaced sequence (with the "pretty" values created above
  # providing the limits for each facet). If there are many facets to
  # show, I recommend no more than 3 labels in each facet, to keep things
  # simple.
  scale_y_continuous(breaks = function(x) seq(from = x[1], 
                                              to = x[2], 
                                              length.out = 3), 
                     expand = c(0, 0))
  elem.bulk1
ggsave("Bulkplots.tiff", plot = elem.bulk1, width=8, height= 10.5, units= "in", dpi =300)



library(tidyr)
Trace_long <- gather(micronutr, Element, Mass, B:Zn, factor_key=TRUE)
Trace_long
#create box plot for all elements
elem.trace<-ggplot()+
  geom_boxplot(data = Trace_long, aes(x=Spp, y =Mass, fill = Tribe), width = .5, size = 1)+
  facet_wrap(~Element, ncol =2, nrow = 3, strip.position = "right", scales = "free_y")+
  ylab(expression(Element~(mg~g^{"-1"})))+
  scale_fill_manual(values =Evpal)+
  ggtitle("Trace elements")+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank(), 
        axis.title = element_blank(),
        axis.text.y = element_text(size = 11, margin=margin(t=-10)),
        axis.text.x = element_text(face = "italic", size = 12, angle = 45, hjust =.75, vjust = .75),
        panel.grid = element_line(color = "grey81"), panel.spacing = unit(.75, "lines"),
        strip.background = element_rect(fill = "grey96"))
elem.trace1<-elem.trace+
  geom_point(data = Trace_long %>% 
               group_by(Element) %>% #group by facet variable
               summarise(y.min = pretty(Mass)[1],
                         y.max = pretty(Mass)[length(pretty(Mass))]) %>%
               tidyr::gather(key, value, -Element), 
             aes(x = 1, y = value),
             inherit.aes = FALSE, alpha = 0) +
  
  # Turn off automatical scale expansion, & manually set scale breaks
  # as an evenly spaced sequence (with the "pretty" values created above
  # providing the limits for each facet). If there are many facets to
  # show, I recommend no more than 3 labels in each facet, to keep things
  # simple.
  scale_y_continuous(breaks = function(x) seq(from = x[1], 
                                              to = x[2], 
                                              length.out = 3), 
                     expand = c(0, 0))
elem.trace1
ggsave("Traceplots.tiff", plot = elem.trace1, width=8, height= 8, units= "in", dpi =300)


elem.plot<-plot_grid(elem.bulk1+ theme(legend.position = "none"),
                     elem.trace1, ncol = 2, align ="h", axis = "b", rel_widths = c(.6,1));elem.plot
ggsave("Elemplots.pdf", plot = elem.plot, width=10, height= 8, units= "in", dpi =300)
attach(micro_long)
library(ggthemes)

################
#lMANOVA on ILRs


# MANOVA test

ILR.man <- manova(cbind(CP, NP, All1, Bulk2,Bulk3, Trace2, Trace3) ~ Spp, data = micronutr)
summary(ILR.man, test = "Wilks", tol = 0) 
summary.aov(ILR.man)
#plot to see how Trace 3 changes with length

ILR.tribe <- manova(cbind(CP, NP, All1, Bulk2, Bulk3,Trace2, Trace3) ~ Tribe, data = micronutr)
summary(ILR.tribe, test = "Wilks", tol = 0) 
summary.aov(ILR.tribe)

ILR.morph <- manova(cbind(CP, NP, All1, Bulk2, Bulk3,Trace2, Trace3) ~ Shell, data = micronutr)
summary(ILR.morph, test = "Wilks", tol = 0) 
summary.aov(ILR.morph,)

ILR.LH <- manova(cbind(CP, NP, All1, Bulk2,Bulk3, Trace2, Trace3) ~ Strategy, data = micronutr)
summary(ILR.LH, test = "Wilks", tol = 0) 
summary.aov(ILR.LH)
dunn.test(CP, Strategy, method="bonferroni", kw = TRUE)
dunn.test(Bulk3, Strategy, method="bonferroni", kw = TRUE)
#==============================================================================================
#PERMANOVA
#==============================================================================================
ILRs<-(as.data.frame(micronutr[23:29]))
attach(ILRs)

#PERMANOVA morphology differences
set.seed(1)
adonis(ILRs~Shell, data=micronutr, permutations = 999, method = "bray",
       strata = NULL, contr.unordered= "contr.sum",
       contr.ordered = "contr.poly")
#BC distance matrix prep step to test for dispersion among species
bc_ILRs<-vegdist(ILRs, method = "bray")
#test for dispersion among groups
disper_ILR_shell<-betadisper(bc_ILRs, Shell, type = "centroid")
anova(disper_ILR_shell)
#PERMANOVA lh strategydifferences
set.seed(1)
adonis(ILRs~Strategy, data=micronutr, permutations = 999, method = "bray",
       strata = NULL, contr.unordered= "contr.sum",
       contr.ordered = "contr.poly")
#BC distance matrix prep step to test for dispersion among species
bc_ILRs<-vegdist(ILRs, method = "bray")
#test for dispersion among groups
disper_ILR_strategy<-betadisper(bc_ILRs, Strategy, type = "centroid")
anova(disper_ILR_strategy)


library(Rmisc)
SummP<-summarySE(data = micronutr, measurevar = "P", groupvars = "Spp");attach(SummP)
SummCP<-summarySE(data = micronutr, measurevar = "CP", groupvars = "Spp");attach(SummCP)
SummGR<-summarySE(data = micronutr, measurevar = "GR", groupvars = "Spp");attach(SummGR)
SummP$GR<-SummGR$GR
SummCP$P<-SummP$P
SummCP$GR<-SummGR$GR

SummP.nolorn<-subset(SummP, Spp != "L.ornata");attach(SummP.nolorn)

ShellPPlot1<-ggplot(data=SummP)+
  geom_smooth(aes(x=GR, y=P), formula = y~x, size = .75, method = "lm", se = FALSE, color = "black")+
  geom_smooth(data = SummP.nolorn, aes(x=GR, y=P), formula = y~x, size = .75, method = "lm", se = FALSE,lty=2, color = "grey40")+
  geom_point(aes(x=GR, y=P, fill = Spp, shape =Spp), color = "black", size = 3.5)+
  geom_errorbar(aes(x=GR, ymin=P-(2*se), ymax=P+(2*se)), width=0, size= .65, color="grey45") + 
  scale_y_continuous(limits = c(0, .4))+
  xlab("Growth Rate (k)")+
  ylab(expression("P (mg " ~ g^-1 *")"))+
  scale_fill_manual(values = c( "#769370", "#BDB2A7", "#F1C646", "#6E687E", "#F17236", "#91D5DE", "#2E8289"), name = "Spp", 
                    breaks = c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor", order = TRUE))+
  scale_shape_manual(values = c( 21, 22, 23, 24, 25, 21, 22), name = "Spp", 
                     breaks = c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor", order = TRUE))+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank())+
  annotate("text", x = 0.25, y = .35, label = "R²=0.13\np=0.01", size = 5)
ShellPPlot1
lmP<-lm(P~GR, data= micronutr)
summary(lmP)
#subset without L.ornata
micronutr.nolorn<-subset(micronutr, Spp != "L.ornata");attach(micronutr.nolorn)
lmP.nolorn<-lm(P~GR, data= micronutr.nolorn)
summary(lmP.nolorn)
##########################################################################
SummCP.nolorn<-subset(SummCP, Spp != "L.ornata");attach(SummCP.nolorn)
ShellCPPlot1<-ggplot(data=SummCP)+
  geom_smooth(aes(x=GR, y=CP), formula = y~x, size = .75, method = "lm", se = FALSE, color = "black")+
  geom_smooth(data = SummCP.nolorn, aes(x=GR, y=CP), formula = y~x, size = .75, method = "lm", se = FALSE,lty=2, color = "grey40")+
  geom_point(aes(x=GR, y=CP, fill = Spp, shape = Spp), color = "black", size = 3.5)+
  geom_errorbar(aes(x=GR, ymin=CP-(2*se), ymax=CP+(2*se)), width=0, size= .65, color="grey45") + 
  scale_y_continuous(limits = c(1.75, 2))+
  labs(x = "Growth Rate (k)", y = "[C|P]")+
  scale_fill_manual(values = c( "#769370", "#BDB2A7", "#F1C646", "#6E687E", "#F17236", "#91D5DE", "#2E8289"), name = "Spp", breaks = c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor", order = TRUE))+
  scale_shape_manual(values = c( 21, 22, 23, 24, 25, 21, 22), name = "Spp", breaks = c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor", order = TRUE))+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank())+
  annotate("text", x = 0.25, y = 1.95, label = "R²=0.20\np=0.003", size = 5)
ShellCPPlot1


SummBulk3<-summarySE(data = micronutr, measurevar = "Bulk3", groupvars = "Spp")
SummCP$Bulk3<-SummBulk3$Bulk3
attach(SummCP)
ShellBulk3Plot1<-ggplot(data=SummCP)+
  geom_smooth(aes(x=GR, y=Bulk3), formula = y~x, size = .75, method = "lm", se = FALSE, color = "black")+
  geom_smooth(data = SummCP.nolorn, aes(x=GR, y=Bulk3), formula = y~x, size = .75, method = "lm", se = FALSE,lty=2, color = "grey40")+
  geom_point(aes(x=GR, y=Bulk3, fill = Spp, shape=Spp), color = "black", size = 3.5)+
  geom_errorbar(aes(x=GR, ymin=Bulk3-(2*se), ymax=Bulk3+(2*se)), width=0, size= .65, color="grey45") + 
  scale_y_continuous(limits = c(11, 12.5))+
  labs(x = "Growth Rate (k)", y = "[C,Ca|P]")+
  scale_fill_manual(values = c( "#769370", "#BDB2A7", "#F1C646", "#6E687E", "#F17236", "#91D5DE", "#2E8289"), name = "Spp", breaks = c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor", order = TRUE))+
  scale_shape_manual(values = c( 21, 22, 23, 24, 25, 21, 22), name = "Spp", breaks = c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor", order = TRUE))+
    theme(legend.position = c(0.72, .25), legend.title = element_blank(), 
        legend.text = element_text(face= "italic", size =10))+
  annotate("text", x = 0.25, y = 12.2, label = "R²=0.21\np=0.002", size = 5)
ShellBulk3Plot1

GR.plot1<-plot_grid(ShellPPlot1,ShellCPPlot1,ShellBulk3Plot1, nrow = 3, align ="v", axis = "l", rel_heights = c(.9,.9, 1.1), labels = c("A","B", "C"));GR.plot1
ggsave("GRplots.tiff", plot = GR.plot1, width=5, height= 10, units= "in", dpi =300)

lmNP<-lm(NP~GR, data= micronutr)
summary(lmNP)

lmCP<-lm(CP~Spp, data= micronutr)
summary(lmCP)
emmeans(lmCP, list(pairwise~Spp), adjust = "Tukey")

lmCP.nolorn<-lm(CP~GR, data= micronutr.nolorn)
summary(lmCP.nolorn)


ShellCPPlot<-ggplot(data=micronutr)+
  geom_smooth(aes(x=GR, y=CP), formula = y~x, size = 1.5, method = "lm", se = TRUE, color = "black")+
  geom_point(aes(x=GR, y=CP, fill = Spp), shape = 21, color = "black", size = 5)+
  scale_y_continuous(limits = c(1.7, 2))+
  labs(x = "Growth Rate (k)", y = "[C|P]")+
  scale_fill_manual(values = c( "#769370", "#BDB2A7", "#F1C646", "#6E687E", "#F17236", "#91D5DE", "#2E8289"), name = "Spp", breaks = c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor", order = TRUE))+
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank())+
  annotate("text", x = 0.25, y = 1.75, label = "R²=0.20\np=0.003", size = 5)

ShellCPPlot



ShellBulk3Plot<-ggplot(data=micronutr)+
  geom_smooth(aes(x=GR, y=Bulk3), formula = y~x, size = 1.5, method = "lm", se = TRUE, color = "black")+
  geom_point(aes(x=GR, y=Bulk3, fill = Spp), shape = 21, color = "black", size = 5)+
  labs(x = "Growth Rate (k)", y = "[C,Ca|P]")+
  scale_fill_manual(values = c( "#769370", "#BDB2A7", "#F1C646", "#6E687E", "#F17236", "#91D5DE", "#2E8289"), name = "Spp", breaks = c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "L.ornata", "O.unicolor", order = TRUE))+
  theme(legend.position = c(0.75, .2), legend.title = element_blank(), legend.text = element_text(face= "italic"))+
  annotate("text", x = 0.25, y = 11, label = "R²=0.21\np=0.002", size = 5)

ShellBulk3Plot

GR.plot1<-plot_grid(ShellCPPlot,ShellBulk3Plot, nrow = 2, align ="v", axis = "l", rel_heights = c(1,1), labels = c("A","B"));GR.plot
ggsave("GRplots.tiff", plot = GR.plot, width=6, height= 10, units= "in", dpi =300)

lmBulk3<-lm(Bulk3~Spp, data = micronutr)
summary(lmBulk3)

lmB3<-lm(Bulk3~Spp, data= micronutr)
emmeans(lmBulk3, list(pairwise~Spp), adjust = "Tukey")






