

----------------------------------------------------------------------------------------------------------------------------------------
#project:Shell micronutrient project 
#coder:G.W. Hopper
----------------------------------------------------------------------------------------------------------------------------------------

######################## 
######Step 1.##########
#######################

#Set working directory and import micronutrient data
setwd("~/Shell/Micronutrients")
micronutr<-read.csv(file="ShellMicronutrientsMass.csv", header = TRUE, sep = ",", row.names = c())
head(micronutr)

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

#creat a data frame of micronutrients
Micros<-(as.data.frame(micronutr[c(6:15)]))
attach(Micros)
#Run principle components analysis
micro.pca<-prcomp(na.omit(log(Micros)), scale =TRUE)
ordiplot(micro.pca)
#check eigen values with scree plot
fviz_eig(micro.pca)
eig.val <- get_eigenvalue(micro.pca)
eig.val
screeplot(micro.pca, bstick = FALSE, type = c("barplot", "lines"),
          npcs = min(10, length(micro.pca$sdev)),
          ptype = "o", bst.col = "red", bst.lty = "solid",
          xlab = "Component", ylab = "Inertia",
          main = deparse(substitute(micro.pca)), legend = bstick)

pca.loadings <- micro.pca$rotation
pca.loadings
#extract PCA scores for axis 1, 2, 3


#need to do this again because more column headers have been added at this point.
attach(micronutr)
head(micronutr)
#instale national parks palettes to make graph pretty as NP
#devtools::install_github("katiejolly/nationalparkcolors")
library(nationalparkcolors)
RWpal<-park_palette("Redwoods")
Evpal<-park_palette("Everglades", 4)

pal<-c(RWpal,Evpal) #colorpaletts


#simple ordinatioN
Microplot<- autoplot(micro.pca, data = micronutr, fill = 'Spp', shape = 21, size = 6,
         loadings = TRUE, loadings.colour = 'grey45', loadings.label = TRUE, loadings.label.size = 7, 
         loadings.label.colour = 'black',loadings.label.hjust = 1.2, loadings.labels.vjust=1)+
  stat_conf_ellipse(aes(color = Spp, fill = Spp), alpha = 0.07, size =.0, geom = "polygon")+
  geom_hline(aes(yintercept=0), linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0), linetype="dashed", size=1)+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  theme(legend.position = c(.01,.15), legend.title = element_blank(), legend.text = element_text(face = "italic"), 
        legend.background = element_rect(linetype = 2, color="black", fill = "grey94"),
        panel.background = element_rect(fill="grey94"));Microplot
#adjustments to loading vectors
Microplot$layers[[2]]$aes_params$size <- 1
Microplot$layers[[2]]$geom_params$arrow$length <- unit(0, units = "points")
Microplot

ggsave("ShellMicrochemPCA.tiff", Microplot, width=8, height= 8, units = "in", dpi =300)


#run a multivariate analysis of variance for 10 elements testing for species effects
elem.man <- manova(cbind(C,N,P, K, B, Zn, Fe, Cu, Ca, Mn) ~ Spp, data = micronutr)
summary.aov(elem.man)
p.adjust(method ="hochberg")

#Load tidyr to use gather.
library(tidyr)


micro_long <- gather(micronutr, Element, Mass, B:Zn, factor_key=TRUE)
micro_long
#create box plot for all elements
elem.plot<-ggplot()+
  geom_boxplot(data = micro_long, aes(x=Spp, y =Mass, fill = Tribe), width = .5, size = 1)+
  facet_wrap(~Element, ncol =2, nrow = 5, strip.position = "right", scales = "free_y")+
  ylab(expression(Element~(mg~g^{"-1"})))+
  scale_fill_manual(values =Evpal)+
 
  theme_bw()+
  theme(legend.position = "bottom", axis.title.x = element_blank(),
        axis.text.y = element_text(size = 11, margin=margin(t=-10)),
        axis.text.x = element_text(face = "italic", size = 14, angle = 45, hjust =.75, vjust = .75),
        panel.grid = element_line(color = "grey81"), panel.spacing = unit(.75, "lines"))
elem.plot1<-elem.plot+
  geom_point(data = micro_long %>% 
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
  elem.plot1
ggsave("Elementplots.tiff", plot = elem.plot1, width=8, height= 10.5, units= "in", dpi =300)

attach(micro_long)
library(ggthemes)

################
#lMANOVA on ILRs
setwd("~/Shell/Micronutrients")
ILR<-read.csv(file="Shell_ILR.csv", header = TRUE, sep = ",", row.names = c())
head(ILR)

#create a factor for species names, use attach to call by name
ILR$Spp<-factor(ILR$Spp, levels=c("A.plicata","C.asperata", "T.verrucosa","F.cerina", "P.decisum", "O.unicolor", "L.ornata", order=TRUE))
ILR$Strategy<-factor(ILR$Strategy, levels=c("Equilibrium","Periodic",  order=TRUE))
ILR$Shell<-factor(ILR$Shell, levels=c("U2D2","Generalized", "Unsculptured", order=TRUE))


# MANOVA test
ILR.man <- manova(cbind(BulkTraceILR1, BulkILR1,BulkILR2,TraceILR2, CN, CP, NP) ~ Spp, data = ILR)
summary.aov(ILR.man)
library(car)
library(emmeans)
#ANOVA for Traditional Stoichiometry ILRs
#Shell CN
CNSpp<-lm(CN~Spp, data = ILR)
emmeans(CNSpp, list(pairwise~Spp), adjust ="Tukey")

CNTribe<-lm(CN~Tribe, data = ILR)
emmeans(CNTribe, list(pairwise~Tribe), adjust ="Tukey")

CNStrat<-lm(CN~Strategy, data = ILR)
emmeans(CNStrat, list(pairwise~Strategy), adjust ="Tukey")

CNShell<-lm(CN~Shell, data = ILR)
emmeans(CNShell, list(pairwise~Shell), adjust ="Tukey")
#############
#Shel CP
CPSpp<-lm(CP~Spp, data = ILR)
emmeans(CPSpp, list(pairwise~Spp), adjust ="Tukey")

CPTribe<-lm(CP~Tribe, data = ILR)
emmeans(CPTribe, list(pairwise~Tribe), adjust ="Tukey")

CPShell<-lm(CP~Shell, data = ILR)
emmeans(CPShell, list(pairwise~Shell), adjust ="Tukey")
############
#Shell NP
NPSpp<-lm(NP~Spp, data = ILR)
emmeans(NPSpp, list(pairwise~Spp), adjust ="Tukey")

NPTribe<-lm(NP~Tribe, data = ILR)
emmeans(NPTribe, list(pairwise~Tribe), adjust ="Tukey")

NPShell<-lm(NP~Shell, data = ILR)
emmeans(NPShell, list(pairwise~Shell), adjust ="Tukey")

#ANOVA for individual Bulk ILR
BulkTrace1<-lm(BulkTraceILR1~Tribe, data = ILR)
emmeans(BulkTrace1, list(pairwise~Tribe), adjust ="Tukey")

Bulk1<-lm(BulkILR1~Spp, data = ILR)
emmeans(Bulk1, list(pairwise~Spp), adjust ="Tukey")

Bulk2Spp<-lm(BulkILR2~Spp, data = ILR)
emmeans(Bulk2Spp, list(pairwise~Spp), adjust ="Tukey")

Bulk2tribe<-lm(BulkILR2~Tribe, data = ILR)
emmeans(Bulk2tribe, list(pairwise~Tribe), adjust ="Tukey")

Bulk2shell<-lm(BulkILR2~Shell, data = ILR)
emmeans(Bulk2shell, list(pairwise~Shell), adjust ="Tukey")
#ANOVA for individual Trace ILR
Trace2Spp<-lm(TraceILR2~Spp, data = ILR)
emmeans(Trace2Spp, list(pairwise~Spp), adjust ="Tukey")

Trace2tribe<-lm(TraceILR2~Tribe, data = ILR)
emmeans(Trace2tribe, list(pairwise~Tribe), adjust ="Tukey")



##############################

#ILR balance value plots
library(tidyr)
library(gghalves)
ILR_long <- gather(ILR, ilr, values, BulkTraceILR1:NP, factor_key=TRUE)
ILR_long

attach(ILR_long)
library(Rmisc)
library(ggthemes)

ILR.plot<-ggplot()+
  geom_boxplot(data = ILR_long, aes(x=Spp, y = values, fill = Tribe), width = .25, size = 1)+
  facet_wrap(~ilr, ncol =1,  strip.position = "right", scales="free_y")+
  labs(x="", y = "Nutrient Balance Value")+
  scale_fill_manual(values =Evpal)+
  theme_bw()+
  theme(legend.position = "bottom",
    axis.text.x = element_text(face = "italic", size = 14, angle = 45, hjust =.75, vjust = .75),
    axis.title.y = element_text(size = 14))

ILR.plot1<-ILR.plot+
  geom_point(data = ILR_long %>% 
               group_by(ilr) %>% #group by facet variable
               summarise(y.min = pretty(values)[-1],
                         y.max = pretty(values)[length(pretty(values))]) %>%
               tidyr::gather(key, value, -ilr), 
             aes(x = -1, y = value),
             inherit.aes = TRUE, alpha = 1) +
  
  # Turn off automatical scale expansion, & manually set scale breaks
  # as an evenly spaced sequence (with the "pretty" values created above
  # providing the limits for each facet). If there are many facets to
  # show, I recommend no more than 3 labels in each facet, to keep things
  # simple.
  scale_y_continuous(breaks = function(x) seq(from = x[1], 
                                              to = x[2], 
                                              length.out = 3), 
                     expand = c(-12, 0))
ILR.plot1
ggsave("ILRplots.tiff", plot = ILR.plot, width=6, height= 12, units= "in", dpi =300)

#################################
ilr.df<-(as.data.frame(ILR[c(6:12)]))
#Run principle components analysis
ilr.pca<-prcomp(ilr.df, scale =TRUE)
ilr.mds<-metaMDS(ilr.df, distance = "mahalanobis")
library(ggvegan)
library(ggordiplots)
gg_ordiplot(ilr.mds, groups = ILR$Spp, scaling = 1, choices = c(1, 2))
#check eigen values with scree plot
fviz_eig(ilr.pca)
eig.val <- get_eigenvalue(ilr.pca)
eig.val
screeplot(ilr.pca, bstick = FALSE, type = c("barplot", "lines"),
          npcs = min(10, length(ilr.pca$sdev)),
          ptype = "o", bst.col = "red", bst.lty = "solid",
          xlab = "Component", ylab = "Inertia",
          main = deparse(substitute(ilr.pca)), legend = bstick)

pca.loadings <- ilr.pca$rotation
pca.loadings
#extract PCA scores for axis 1, 2, 3
res.ind <- get_pca_ind(ilr.pca)

pca.scores <- res.ind$coord
as.data.frame(pca.scores)

coord.groups <- res.ind$coord %>%
  as_tibble() %>%
  select(Dim.1, Dim.2, Dim.3)



#need to do this again because more column headers have been added at this point.
attach(ILR)
head(ILR)
ILR.df<-as.data.frame(ILR)
#instale national parks palettes to make graph pretty as NP
#devtools::install_github("katiejolly/nationalparkcolors")
library(nationalparkcolors)
RWpal<-park_palette("Redwoods")
Evpal3<-park_palette("Everglades", c(5:7))

pal<-c(RWpal,Evpal) #colorpaletts



attach(ILR.df)
ILR.PCplotSpp<- autoplot(ilr.pca, data = ILR.df, fill = 'Spp', shape = 21, size = 6,
                      loadings = TRUE, loadings.colour = 'grey45', loadings.label = TRUE, loadings.label.size = 5, 
                      loadings.label.colour = 'black',loadings.labels.hjust= 2, loadings.labels.vjust=1)+
  stat_conf_ellipse(aes(color = Spp, fill = Spp), alpha = 0.1, size =.0, geom = "polygon")+
  geom_hline(aes(yintercept=0), linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0), linetype="dashed", size=1)+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  #guides(color=guide_legend(nrow=2))+
  theme(legend.position = c(.01,.8), legend.title = element_blank(), legend.text = element_text(face = "italic", size =12), 
        legend.background = element_blank(),
        panel.background = element_rect(fill="grey94"));ILR.PCplotSpp
#adjustments to loading vectors
ILR.PCplotSpp$layers[[2]]$aes_params$size <- 1
ILR.PCplotSpp$layers[[2]]$geom_params$arrow$length <- unit(0, units = "points")
ILR.PCplotSpp

ggsave("ShellILRPCSpp.tiff", ILR.PCplotSpp, width=8, height= 8, units = "in", dpi =300)



#simple ordinatioN

ILR.PCplotTribe<- autoplot(ilr.pca, data = ILR.df, fill = 'Tribe', shape = 21, size = 6,
                     loadings = TRUE, loadings.colour = 'grey45', loadings.label = TRUE, loadings.label.size = 5, 
                     loadings.label.colour = 'black',loadings.labels.hjust= 2, loadings.labels.vjust=1)+
  stat_conf_ellipse(aes(color = Tribe, fill = Tribe), alpha = 0.1, size =.0, geom = "polygon")+
  geom_hline(aes(yintercept=0), linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0), linetype="dashed", size=1)+
  scale_color_manual(values = Evpal)+
  scale_fill_manual(values = Evpal)+
  theme(legend.position = c(.025,.85), legend.title = element_blank(), legend.text = element_text(face = "italic"), 
        legend.background = element_blank(),
        panel.background = element_rect(fill="grey94"));ILR.PCplot
#adjustments to loading vectors
ILR.PCplotTribe$layers[[2]]$aes_params$size <- 1
ILR.PCplotTribe$layers[[2]]$geom_params$arrow$length <- unit(0, units = "points")
ILR.PCplotTribe

ggsave("ShellILRPCTribe.tiff", ILR.PCplotTribe, width=8, height= 8, units = "in", dpi =300)

RWpal3<-park_palette("Redwoods", 3)
ILR.PCplotShell<- autoplot(ilr.pca, data = ILR.df, fill = 'Shell', shape = 21, size = 6,
                      loadings = TRUE, loadings.colour = 'grey45', loadings.label = TRUE, loadings.label.size = 5, 
                      loadings.label.colour = 'black',loadings.labels.hjust= 2, loadings.labels.vjust=1)+
  stat_conf_ellipse(aes(color = Shell, fill = Shell), alpha = 0.1, size =.0, geom = "polygon")+
  geom_hline(aes(yintercept=0), linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0), linetype="dashed", size=1)+
  scale_color_manual(values = RWpal3)+
  scale_fill_manual(values = RWpal3)+
  theme(legend.position = c(.02,.85), legend.title = element_blank(), legend.text = element_text(face = "italic"), 
        legend.background = element_blank(),
        panel.background = element_rect(fill="grey94"));ILR.PCplotShell
#adjustments to loading vectors
ILR.PCplotShell$layers[[2]]$aes_params$size <- 1
ILR.PCplotShell$layers[[2]]$geom_params$arrow$length <- unit(0, units = "points")
ILR.PCplotShell

ggsave("ShellILRPCShell.tiff", ILR.PCShell, width=8, height= 8, units = "in", dpi =300)

#######
RWpal2<-park_palette("Redwoods")

ILR.PCplotStrat<- autoplot(ilr.pca, data = ILR.df, fill = 'Strategy', shape = 21, size = 6,
                           loadings = TRUE, loadings.colour = 'grey45', loadings.label = TRUE, loadings.label.size = 5, 
                           loadings.label.colour = 'black',loadings.labels.hjust= 2, loadings.labels.vjust=1)+
  stat_conf_ellipse(aes(color = Strategy, fill = Strategy), alpha = 0.1, size =.0, geom = "polygon")+
  geom_hline(aes(yintercept=0), linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0), linetype="dashed", size=1)+
  scale_color_manual(values = RWpal2, c(5,6))+
  scale_fill_manual(values = RWpal2, c(5,6))+
  theme(legend.position = c(.05,.85), legend.title = element_blank(), legend.text = element_text(face = "italic"), 
        legend.background = element_blank(),
        panel.background = element_rect(fill="grey94"));ILR.PCplotStrat
#adjustments to loading vectors
ILR.PCplotStrat$layers[[2]]$aes_params$size <- 1
ILR.PCplotStrat$layers[[2]]$geom_params$arrow$length <- unit(0, units = "points")
ILR.PCplotStrat

ggsave("ShellILRPCStrat.tiff", ILR.PCplotStrat, width=8, height= 8, units = "in", dpi =300)

ILR.PCplot<-plot_grid(ILR.PCplotSpp, ILR.PCplotTribe, ILR.PCplotShell, ILR.PCplotStrat, nrow = 2, ncol = 2, labels = "AUTO", align = "vh", rel_heights = c(1, 1, 1,1));ILR.PCplot
ILR.PCplot
save_plot("ILRFigure.tiff", ILR.PCplot, base_width = 10, base_height =10, dpi =300)


library(MASS)
fit<-lda(Spp~C + N + P + K + B + Zn + Fe + Cu + Ca + Mn, data = micronutr)
fit
ct <- table(micronutr$Spp, fit$class)
sum(diag(prop.table(ct)))
plot(fit)
lda.values <- predict(fit)
newdata <- data.frame(type = micronutr[,1], lda = lda.values$x)
ggplot(newdata) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 2.5)
library(caret)
training.samples <- micronutr$Spp%>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- micronutr[training.samples, ]
test.data <- micronutr[-training.samples, ]
# Estimate preprocessing parameters
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)
library(MASS)
# Fit the model
model <- lda(Spp~C + N + P + K + B + Zn + Fe + Cu + Ca + Mn, data = train.transformed)
# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy
mean(predictions$class==test.transformed$Spp)

model
plot(model)
lda.data <- cbind(train.transformed, predict(model)$x)

ggord(model, lda.data$Spp, size = 0, ellipse= FALSE, arrow = 0, veccol = "grey 40",
      vectyp = "solid", veclsz = 1, ext = 2.2, txt=8)+
  geom_hline(aes(yintercept=0), linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0), linetype="dashed", size=1)+
  stat_conf_ellipse(aes(colour = lda.data$Spp, fill = lda.data$Spp), alpha = 0.07, size =1, geom = "polygon")+
  geom_point(aes(fill = lda.data$Spp), shape = 21, size = 6, alpha = .7)+
  scale_color_manual(values = pal, c(5,6))+
  scale_fill_manual(values = pal, c(5,6))+
  theme(legend.position = c(.8,.2), legend.title = element_blank(), legend.text = element_text(face = "italic"), legend.background = element_blank(),
        panel.background = element_rect(fill="grey94"), panel.grid = element_line(color = "grey76"),
        axis.text = element_text(size =12), axis.title = element_text(size = 14))

ggplot(lda.data, aes(LD1, LD2)) +
  geom_point(aes(fill = Spp), shape = 21, size = 6)+
  stat_conf_ellipse(aes(color = Spp, fill = Spp), alpha = 0.07, size =.25, geom = "polygon")+
  geom_hline(aes(yintercept=0), linetype="dashed", size=1)+
  geom_vline(aes(xintercept=0), linetype="dashed", size=1)+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  theme(legend.position = c(.8,.15), legend.title = element_blank(), legend.text = element_text(face = "italic"), 
        legend.background = element_rect(linetype = 2, color="black", fill = "grey94"),
        panel.background = element_rect(fill="grey96"))

