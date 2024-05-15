###### 

##################################

##### clean up 
rm(list=ls())

### load packages 
library(vctrs)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(devtools)
library(zCompositions)
library(tibble)
library(ggplot2)
library(dplyr)
library(vegan)
library(picante)
library(gridExtra)
library(pairwiseAdonis)
library(lme4)
library(lmerTest)
library(microbiome)
library(DT)
library(eulerr)
library(microbiomeutilities)
library(ggrepel)
library(sjPlot)
library(sjmisc)
library(MuMIn)
library(jtools)
library(car)
library(plyr)
library(scales)
library(mgcv)
library(tidygam)
library(tidymv)

# Finally, set a random seed for reproducibility
set.seed(777)

###Upload ASV, taxonomy and meta data and compile into phyloseq object thats saved as RDS file (based on previous filtering)
ps = readRDS("~/Dropbox (Personal)/CellarBats/working_data/CellarBats_ps_filtered_100reads.rds")

### some information about the phyloseq object and make meta df

microbiome::summarize_phyloseq(ps)

sd(sample_sums(ps)) ## on average we still have more than 30 000 reads from our samples 

meta<- as(sample_data(ps), "data.frame")


######## now for the microbiome ########

ps_class<-tax_glom(ps, taxrank = "Class")

sample_data(ps)<-sample_data(ps)[,c("sample.ID")]
ps_class <- microbiome::transform(ps_class, "compositional") 
class_per_sample<-psmelt(ps_class)

ddply(class_per_sample, c("Class"), summarise, mean=round(mean(Abundance),4))


ps_genus<-tax_glom(ps, taxrank = "Genus")

ps_genus <- microbiome::transform(ps_genus, "compositional") 
genus_per_sample<-psmelt(ps_genus)

df_genus<-ddply(genus_per_sample, c("Genus"), summarise, mean=round(mean(Abundance),4))

##### find most common classes ####

top<-as_tibble(get_group_abundances(ps_class, level="Class", group="Class", transform = "compositional") %>%arrange (-mean_abundance))
print(top,n=15)

##### classes that are more prominent than 1%
top_class<-c("Bacilli",  "Actinobacteria", "Gammaproteobacteria", "Alphaproteobacteria", "Clostridia")


##### now top genera
ps_genus<-tax_glom(ps, taxrank = "Genus")

top<-as_tibble(get_group_abundances(ps_genus, level="Genus", group="Genus", transform = "compositional") %>%arrange (-mean_abundance))
print(top,n=15)

top_genus<-c("Weissella",  "Staphylococcus", "Fructobacillus", "Corynebacterium", "Streptococcus", "Actinomyces", "Gemella")


##### relative abundance plots by time point and treatment ######
ps = readRDS("~/Dropbox (Personal)/CellarBats/working_data/CellarBats_ps_filtered_100reads.rds")

treat<-subset_samples(ps, treatment =="T")
cont<-subset_samples(ps, treatment =="C")

ps_treat_time<-merge_samples(treat,"time_point", fun=mean)
ps_treat_time_class<-tax_glom(ps_treat_time, taxrank = "Class")
ps_treat_time_class <- microbiome::transform(ps_treat_time_class, "compositional")
class_by_treat_time<-psmelt(ps_treat_time_class)

class_by_treat_time$Class_plot<-as.factor(ifelse(class_by_treat_time$Class %in% top_class, class_by_treat_time$Class, "Other"))

class_by_treat_time$Class_plot<-factor(class_by_treat_time$Class_plot, 
                                       levels = c("Bacilli",  "Actinobacteria", "Gammaproteobacteria", "Alphaproteobacteria", "Clostridia", "Other"))

class_by_treat_time$Sample<-as.factor(class_by_treat_time$Sample)

summary_fig1A<-ggplot(class_by_treat_time, aes(x =fct_relevel(Sample, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                               y = Abundance, fill = Class_plot))+
  geom_col(size = 0, position = "fill")+ 
  theme_classic()+scale_fill_manual(values=c( "#B9D5C1", "#88B7B5", "#EABAED", "#C067C7", "#A886DF","grey65"))+
  theme(axis.text.x = element_text(size = 10, angle=90, vjust=0.5), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.title = element_text(colour = 'black', face="bold", size=12))+labs(fill ="Bacterial Class")+
  ylab("Relative Abundance")+xlab("day time (24hours)")+
  theme(axis.title.x = element_blank())+scale_y_continuous(expand = c(0,0))


ps_cont_time<-merge_samples(cont,"time_point", fun=mean)
ps_cont_time_class<-tax_glom(ps_cont_time, taxrank = "Class")
ps_cont_time_class <- microbiome::transform(ps_cont_time_class, "compositional")
class_by_cont_time<-psmelt(ps_cont_time_class)

class_by_cont_time$Class_plot<-as.factor(ifelse(class_by_cont_time$Class %in% top_class, class_by_cont_time$Class, "Other"))

class_by_cont_time$Class_plot<-factor(class_by_cont_time$Class_plot, 
                                      levels = c("Bacilli",  "Actinobacteria", "Gammaproteobacteria", "Alphaproteobacteria", "Clostridia", "Other"))

class_by_cont_time$Sample<-as.factor(class_by_cont_time$Sample)

summary_fig1B<-ggplot(class_by_cont_time, aes(x =fct_relevel(Sample, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                              y = Abundance, fill = Class_plot))+
  geom_col(size = 0, position = "fill")+ 
  theme_classic()+scale_fill_manual(values=c("#B9D5C1", "#88B7B5", "#EABAED", "#C067C7", "#A886DF","grey65"))+
  theme(axis.text.x = element_text(size = 10, angle=90, vjust=0.5), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.title = element_text(colour = 'black', face="bold", size=12))+labs(fill ="Bacterial Class")+
  ylab("Relative Abundance")+
  theme(axis.title.x = element_blank())+scale_y_continuous(expand = c(0,0))


ggarrange(
  summary_fig1B, summary_fig1A, 
  labels = c("A", "B"),
  nrow=2, align = c("h"),
  common.legend = TRUE, legend = "bottom"
)

#########  same for genera ######

ps_treat_time<-merge_samples(treat,"time_point", fun=mean)
ps_treat_time_genus<-tax_glom(ps_treat_time, taxrank = "Genus")
ps_treat_time_genus <- microbiome::transform(ps_treat_time_genus, "compositional")
genus_by_treat_time<-psmelt(ps_treat_time_genus)

genus_by_treat_time$Genus_plot<-as.factor(ifelse(genus_by_treat_time$Genus %in% top_genus, genus_by_treat_time$Genus, "Other"))

genus_by_treat_time$Genus_plot<-factor(genus_by_treat_time$Genus_plot, 
                                       levels = c("Weissella",  "Staphylococcus", "Fructobacillus", "Corynebacterium", "Streptococcus", "Actinomyces", "Gemella", "Other"))

genus_by_treat_time$Sample<-as.factor(genus_by_treat_time$Sample)

summary_figXA<-ggplot(genus_by_treat_time, aes(x =fct_relevel(Sample, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                               y = Abundance, fill = Genus_plot))+
  geom_col(size = 0, position = "fill")+ 
  theme_classic()+scale_fill_manual(values=c( "#A97092", "#C28CAE", "#DEC4A1", "#EBBB54","#EDEB7A", "#A3D186", "#76B19B","grey65"))+
  theme(axis.text.x = element_text(size = 10, angle=90, vjust=0.5), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.title = element_text(colour = 'black', face="bold", size=12))+labs(fill ="Bacterial Class")+
  ylab("Relative Abundance")+xlab("day time (24hours)")+
  theme(axis.title.x = element_blank())+scale_y_continuous(expand = c(0,0))


ps_cont_time<-merge_samples(cont,"time_point", fun=mean)
ps_cont_time_genus<-tax_glom(ps_cont_time, taxrank = "Genus")
ps_cont_time_genus <- microbiome::transform(ps_cont_time_genus, "compositional")
genus_by_cont_time<-psmelt(ps_cont_time_genus)

genus_by_cont_time$Genus_plot<-as.factor(ifelse(genus_by_cont_time$Genus %in% top_genus, genus_by_cont_time$Genus, "Other"))

genus_by_cont_time$Genus_plot<-factor(genus_by_cont_time$Genus_plot, 
                                      levels = c("Weissella",  "Staphylococcus", "Fructobacillus", "Corynebacterium", "Streptococcus", "Actinomyces", "Gemella", "Other"))

genus_by_cont_time$Sample<-as.factor(genus_by_cont_time$Sample)

summary_figXB<-ggplot(genus_by_cont_time, aes(x =fct_relevel(Sample, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                              y = Abundance, fill = Genus_plot))+
  geom_col(size = 0, position = "fill")+ 
  theme_classic()+scale_fill_manual(values=c("#A97092", "#C28CAE", "#DEC4A1", "#EBBB54", "#EDEB7A", "#A3D186", "#76B19B", "grey65"))+
  theme(axis.text.x = element_text(size = 10, angle=90, vjust=0.5), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.title = element_text(colour = 'black', face="bold", size=12))+labs(fill ="Bacterial Class")+
  ylab("Relative Abundance")+
  theme(axis.title.x = element_blank())+scale_y_continuous(expand = c(0,0))


ggarrange(
  summary_figXB, summary_figXA, 
  labels = c("A", "B"),
  nrow=2, align = c("h"),
  common.legend = TRUE, legend = "bottom"
)


#####  alpha diversity metrices and some stats ######
#####################################################

#first we rarefy 

total = 250

standf = function(x, t=total) round(t * (x / sum(x))) #standardise by this sampling depth 

ps.rare = transform_sample_counts(ps, standf) #normalise data accordingly 


otu.table <- as.data.frame(otu_table(ps.rare))
meta <- as(sample_data(ps.rare), "data.frame")
df.pd <- pd(t(otu.table), phy_tree(ps.rare), include.root=T)
meta$Phyogenetic_diversity <- df.pd$PD

alpha = estimate_richness(ps.rare, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")) #generate df with various diversity metices 

#add indicies to meta data (like we did with PD)
meta$Shannon <- alpha$Shannon 
meta$Chao1 <- alpha$Chao1 
meta$Observed <- alpha$Observed

#### test for correlation between unrarefied and rarefied data
alpha.unrare = estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")) #generate df with various diversity metices 

#add indicies to meta data (like we did with PD)
meta$Shannon.unrare <- alpha.unrare$Shannon 
meta$Chao1.unrare <- alpha.unrare$Chao1
meta$Observed.unrare <- alpha.unrare$Observed

ggplot(meta, aes(Shannon.unrare, Shannon))+geom_point()+geom_smooth(method="lm")
ggplot(meta, aes(Chao1.unrare, Chao1))+geom_point()+geom_smooth(method="lm")
ggplot(meta, aes(Observed.unrare, Observed))+geom_point()+geom_smooth(method="lm")

cor.test(meta$Observed.unrare, meta$Observed)

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  hjustvar = c(1.1),
  vjustvar = c(2))

SupFig3<-ggplot(meta, aes(Observed.unrare, Observed))+geom_point()+geom_smooth(method="lm")+theme_bw()+
ylab("rarefied Observed ASVs")+xlab("unrarefied Observed ASVs")+
theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), axis.title.y = element_text(colour = 'black', face="bold", size=12))+
geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=63.3, p<0.001"))


############## Alpha GAMS #############
######### OBSERVED ASVs #########

SupFig5a<-ggplot(meta, aes(fct_relevel(time_point, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                 log(Observed.unrare), fill=treatment)) + 
  geom_boxplot( outlier.alpha = 0) + geom_point(shape=21, position=position_dodge(.8))+ theme_classic() +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547")) + xlab("day time (24hours)") + ylab("Observed ASVs")+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position= "none")


##### gam observed ###
meta.gam<-meta #### because time point needs to be numerical for gams
meta.gam$time_point<-as.numeric(meta.gam$time_point)
meta.gam$Seq_depth_scaled<-as.numeric(rescale(sample_sums(ps),, to =c(-1,1)))

Observed_gam <- mgcv::gam(log(Observed.unrare)~
                            treatment+ day + s(Seq_depth_scaled) +
                            s(time_point, by = treatment, bs = "cc", k=6),
                          data=meta.gam,
                          family = gaussian, REML=T)

print(summary(Observed_gam)) # no treatment effect (great)
gam.check(Observed_gam)
plot(Observed_gam)


SupFig5b<-tidymv::plot_smooths(
  model = Observed_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","solid"))+
  ylab("Observed ASVs") + xlab("day time (24 hours)") +
  scale_x_continuous(breaks = round(seq(min(as.numeric(meta.gam$time_point)), max(as.numeric(meta.gam$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position="none")


###### SHANNON ######

Fig3a<-ggplot(meta, aes(fct_relevel(time_point, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                         (Shannon.unrare), fill=treatment)) + 
  geom_boxplot( outlier.alpha = 0) + geom_point(shape=21, position=position_dodge(.75))+ theme_classic() +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547")) + xlab("day time (24hours)") + ylab("Shannon")+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position="none")

### shannon gam ###
Shannon_gam <- mgcv::gam((Shannon.unrare)~
                           treatment+ day + s(Seq_depth_scaled) +
                           s(time_point, by = treatment, bs = "cc", k=6),
                         data=meta.gam,
                         family = gaussian, REML=T)

print(summary(Shannon_gam)) # no treatment effect (great)
gam.check(Shannon_gam)
plot(Shannon_gam)


Shannon_plot<-tidymv::plot_smooths(
  model = Shannon_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","dashed"))+
  #ylab("Shannon") + 
  xlab("day time (24 hours)") +
  scale_x_continuous(breaks = round(seq(min(as.numeric(meta.gam$time_point)), max(as.numeric(meta.gam$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none")


ggarrange(SupFig5a, SupFig5b, nrow=1, labels=c("A","B"))


######## beta diversity index ######

# calculate weighted Unifrag distance using the phyloseq package
dist_wUni = phyloseq::distance(ps.rare, method="weighted unifrac") ### abundance matters + phylogeny matters 
# calculate unweighted distance using the phyloseq package
dist_Uni = phyloseq::distance(ps.rare, method="unifrac") ### presence/absence + phylogeny matters

meta <- as(sample_data(ps.rare), "data.frame")

meta$Seq_depth_scaled<-as.numeric(rescale(sample_sums(ps), to =c(-1,1)))

adonis2(dist_Uni ~ treatment*time_point + day +Seq_depth_scaled, data = meta, method="unifrac", by="margin")
adonis2(dist_wUni ~ treatment*time_point + day + Seq_depth_scaled, data = meta, method="unifrac", by="margin")

######### for visualisation we keep the treatment and control group seperately
control<-subset_samples(ps.rare, treatment=="C")
meta_c<-sample_data(control)

dist_Uni_c = phyloseq::distance(control, method="unifrac") ### presence/absence + phylogeny matters
ordination_Uni_c<-ordinate(control, method="PCoA", distance=dist_Uni_c)

unweighted1 <- plot_ordination(control, ordination_Uni_c) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = time_point)) + # put ellipses around group centroids
  geom_point(aes(
    fill = time_point), # fill color by infection status
    size = 4, # make points size 4
    pch = 21, # Make points circular with a border 
    col = "black") + # the border colour should be black
  #labs(subtitle = "a)") + # insert subtitle
  #ggtitle("Unweighted Unifrac") + # insert title
  theme_bw() + 
  scale_fill_manual(values=c("#8A3D9C", "#A048A9", "#B552B5", "#BD64A3", "#C57691", "#CD887F",
                             "#D59A6D", "#DDAC5B", "#E5BE49", "#EDCF37", "#F3DD6B","#EDE3B4")) +
  labs(fill = "time point") +
  theme(plot.title = element_text(hjust = 0.5))+facet_wrap(~treatment)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold")) # center align the plot title
unweighted1


#Unweighted Unifrac Treatment Group
treatment<-subset_samples(ps.rare, treatment=="T")
meta_t<-sample_data(treatment)

dist_Uni_t = phyloseq::distance(treatment, method="unifrac") ### presence/absence + phylogeny matters
ordination_Uni_t<-ordinate(treatment, method="PCoA", distance=dist_Uni_t)

unweighted2<-plot_ordination(treatment, ordination_Uni_t)+ 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = time_point))+ # put ellipses around group centroids
  geom_point(aes(
    fill = time_point), #fill color by infection status
    size = 4, # make points size 4
    pch = 21, #Make points circular with a border 
    col = "black") +# the border colour should be black
  #labs(subtitle = "c)") + # insert subtitle
  #ggtitle(" ") + # insert title
  theme_bw() + 
  scale_fill_manual(values=c("#8A3D9C", "#A048A9", "#B552B5", "#BD64A3", "#C57691", "#CD887F",
                             "#D59A6D", "#DDAC5B", "#E5BE49", "#EDCF37", "#F3DD6B","#EDE3B4")) +
  labs(fill = "time point") +
  theme(plot.title = element_text(hjust = 0.5))+facet_wrap(~treatment)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold"))
unweighted2


##### weighted unifrac ######  
meta.control <- as(sample_data(control), "data.frame")
dist_wUni_c = phyloseq::distance(control, method="wunifrac") ### abundance matters + phylogeny matters
ordination_wUni_c<-ordinate(control, method="PCoA", distance=dist_wUni_c)

weighted1<-plot_ordination(control, ordination_wUni_c)+ 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = time_point))+ # put ellipses around group centroids
  geom_point(aes(
    fill = time_point), #fill color by infection status
    size = 4, # make points size 4
    pch = 21, #Make points circular with a border 
    col = "black") +# the border colour should be black
  #labs(subtitle = "b)") + # insert subtitle
  #ggtitle("Weighted Unifrac") + # insert title
  theme_bw() + 
  scale_fill_manual(values=c("#8A3D9C", "#A048A9", "#B552B5", "#BD64A3", "#C57691", "#CD887F",
                             "#D59A6D", "#DDAC5B", "#E5BE49", "#EDCF37", "#F3DD6B","#EDE3B4")) +
  labs(fill = "time point") +
  theme(plot.title = element_text(hjust = 0.5))+facet_wrap(~treatment)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold"))
weighted1



#Weighted Unifrac Treatment Group
meta.treatment <- as(sample_data(treatment), "data.frame")
dist_wUni_t = phyloseq::distance(treatment, method="wunifrac") ### presence/absence + phylogeny matters
ordination_wUni_t<-ordinate(treatment, method="PCoA", distance=dist_wUni_t)

weighted2<-plot_ordination(treatment, ordination_wUni_t)+ 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(fill = time_point))+ # put ellipses around group centroids
  geom_point(aes(
    fill = time_point), #fill color by infection status
    size = 4, # make points size 4
    pch = 21, #Make points circular with a border 
    col = "black") +# the border colour should be black
  #labs(subtitle = "d)") + # insert subtitle
  #ggtitle(" ") + # insert title
  theme_bw() + 
  scale_fill_manual(values=c("#8A3D9C", "#A048A9", "#B552B5", "#BD64A3", "#C57691", "#CD887F",
                             "#D59A6D", "#DDAC5B", "#E5BE49", "#EDCF37", "#F3DD6B","#EDE3B4")) +
  labs(fill = "time point") +
  theme(plot.title = element_text(hjust = 0.5))+facet_wrap(~treatment)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', face="bold"))
weighted2


#together 
ggarrange(unweighted1, unweighted2, weighted1, weighted2, ncol = 2, nrow=2, common.legend=TRUE, legend="right")


####### pairwise permanova #####
meta$timeBYtreatment<- paste0(meta$time_point,meta$treatment)

result_list_weighted<-pairwise.adonis2(phyloseq::distance(ps.rare, method = "weighted unifrac") ~ timeBYtreatment, 
                                       data = meta)
result_list_unweighted<-pairwise.adonis2(phyloseq::distance(ps.rare, method = "unifrac") ~ timeBYtreatment, 
                                         data = meta)


#### PC1 plot ######
PCoA_results_c_unweighted<-as.data.frame(ordination_Uni_c$vectors)
meta.control$PC1_unweighted<-PCoA_results_c_unweighted$Axis.1

PCoA_results_t_unweighted<-as.data.frame(ordination_Uni_t$vectors)
meta.treatment$PC1_unweighted<-PCoA_results_t_unweighted$Axis.1

PCoA_results_c<-as.data.frame(ordination_wUni_c$vectors)
meta.control$PC1<-PCoA_results_c$Axis.1

PCoA_results_t<-as.data.frame(ordination_wUni_t$vectors)
meta.treatment$PC1<-PCoA_results_t$Axis.1

#### fuse back together and plot
meta.PCoA<-rbind(meta.control, meta.treatment)

##### visualise using box plots
PC1_boxplot<-ggplot(meta.PCoA, aes(fct_relevel(as.factor(time_point), "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                   PC1_unweighted, fill=treatment)) + 
  geom_boxplot(outlier.alpha = 0) + geom_point(shape=21, position=position_dodge(.75))+ theme_classic() +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547")) + xlab("day time (24hours)") + ylab("PC1 scores")+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))


##### gam PC1 ###
meta.PCoA$time_point<-as.numeric(meta.PCoA$time_point)
meta.PCoA$Seq_depth_scaled<-as.numeric(rescale(sample_sums(ps), to =c(-1,1)))

PC1_gam <- mgcv::gam((PC1_unweighted)~
                       treatment+
                       s(time_point, by = treatment, bs = "cc", k=6) + # specifically interested in this variable
                       s(Seq_depth_scaled)+
                       day,
                     data=meta.PCoA,
                     family = gaussian, REML=T)

print(summary(PC1_gam)) # no treatment effect (great)
gam.check(PC1_gam)
plot(PC1_gam)


PC1_plot<-tidymv::plot_smooths(
  model = PC1_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","solid"))+
  ylab("PC1 scores") + xlab("day time (24 hours)") +
  scale_x_continuous(breaks = round(seq(min((meta.PCoA$time_point)), max((meta.PCoA$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))



##### constrained coordination analysis #######

#### weighted
final_model<-capscale(dist_wUni ~  meta$time_point * meta$treatment + meta$day + meta$Seq_depth_scaled)

library(ggord)

ggord(final_model,  
      xlims = c(min(data.frame(final_model$CCA$wa)$CAP1)-.8, 
                max(data.frame(final_model$CCA$wa)$CAP1)+.8),
      ylims = c(min(data.frame(final_model$CCA$wa)$CAP2)-.8, 
                max(data.frame(final_model$CCA$wa)$CAP2)+.8), 
      addsize = 3, size = 1, txt = 2) ### y = 24.1; x = 33.4

final_model_anova<-anova.cca(final_model, by="terms", permutations = 1000)

final_model_df<-vegan::scores(final_model)

vectors_df<-data.frame(final_model_df$sites)

vectors_df$sample.ID<-row.names(vectors_df)

sample_metadata<-data.frame(sample_data(ps.rare))[,c("sample.ID", "time_point","treatment")]
vectors_df<-merge(vectors_df, sample_metadata, by = "sample.ID")

vectors_df_treatment<-subset(vectors_df, treatment=="T")
weighted_Supp6d<-ggplot(vectors_df_treatment, aes(x = CAP1, y = CAP2))+
  geom_point(aes(fill = time_point, shape=treatment), pch=21, size = 4)+
  scale_fill_viridis_d()+
  theme_bw()

vectors_df_control<-subset(vectors_df, treatment=="C")
weighted_Supp6c<-ggplot(vectors_df_control, aes(x = CAP1, y = CAP2))+
  geom_point(aes(fill = time_point, shape=treatment), pch=21, size = 4)+
  scale_fill_viridis_d()+
  theme_bw()

#### unweighted

final_model<-capscale(dist_Uni ~  meta$time_point * meta$treatment + meta$day + meta$Seq_depth_scaled)

library(ggord)

ggord(final_model,  
      xlims = c(min(data.frame(final_model$CCA$wa)$CAP1)-.8, 
                max(data.frame(final_model$CCA$wa)$CAP1)+.8),
      ylims = c(min(data.frame(final_model$CCA$wa)$CAP2)-.8, 
                max(data.frame(final_model$CCA$wa)$CAP2)+.8), 
      addsize = 3, size = 1, txt = 2) ### y = 24.1; x = 33.4

final_model_anova<-anova.cca(final_model, by="terms", permutations = 1000)

final_model_df<-vegan::scores(final_model)

vectors_df<-data.frame(final_model_df$sites)

vectors_df$sample.ID<-row.names(vectors_df)

sample_metadata<-data.frame(sample_data(ps.rare))[,c("sample.ID", "time_point","treatment")]
vectors_df<-merge(vectors_df, sample_metadata, by = "sample.ID")

vectors_df_treatment<-subset(vectors_df, treatment=="T")
unweighted_Supp6b<-ggplot(vectors_df_treatment, aes(x = CAP1, y = CAP2))+
  geom_point(aes(fill = time_point, shape=treatment), pch=21, size = 4)+
  scale_fill_viridis_d()+
  theme_bw()

vectors_df_control<-subset(vectors_df, treatment=="C")
unweighted_Supp6a<-ggplot(vectors_df_control, aes(x = CAP1, y = CAP2))+
  geom_point(aes(fill = time_point, shape=treatment), pch=21, size = 4)+
  scale_fill_viridis_d()+
  theme_bw()

#together 
ggarrange(unweighted_Supp6a, unweighted_Supp6b, weighted_Supp6c, weighted_Supp6d, 
          ncol = 2, nrow=2, common.legend=TRUE, legend="right", labels= c("A", "B", "C", "D"))

##### GAMS ########
#continuous variables typically s(var, ...)
#categories need to be factors 
#interactions coded as s(var, by = interaction_term)
#k = number of basic functions fitted to a line 
#REML or use sp = # to adjust smoothing 
#bs = 

#edf is effective degree of freedoms
#edf= 1 is a linear fit 
#edf= 2 is a quadratic term 
#edf> 2 the more wiggly the line becomes

#a significant smooth term is one where you can not draw a horizontal line through the 95% confidence interval.


### prep for GAMMs
top_lepto<-microbiomeutilities::aggregate_top_taxa2(ps, "Genus", top = 30)
print(otu_table(top_lepto),n=41)
clr_scaled <-microbiome::transform(top_lepto, transform = "clr")

taxa_to_keep <- c("Weissella", "Staphylococcus", "Fructobacillus", 
                  "Corynebacterium", "Streptococcus", "Actinomyces", "Gemella") #### 7 top taxa

lepto_genus_scaled<-prune_taxa(taxa_to_keep, clr_scaled) 
sample_data(lepto_genus_scaled)<-sample_data(lepto_genus_scaled)[,c( "sample.ID")]
top_genera_melt<-psmelt(lepto_genus_scaled)
head(top_genera_melt)
meta.PCoA<-data.frame(sample_data(ps))
meta.PCoA$Seq_depth<-as.numeric(sample_sums(top_lepto))
meta.PCoA$Seq_depth_scaled <- as.numeric(scales::rescale(meta.PCoA$Seq_depth, to = c(-1,1)))
top_lepto_genera<- merge(top_genera_melt, meta.PCoA, by.x = "Sample", by.y = "sample.ID")
names(top_lepto_genera)
head(top_lepto_genera)
top_lepto_genera$time_point<-as.numeric(top_lepto_genera$time_point)
top_lepto_genera$time_point<-as.numeric(top_lepto_genera$time_point)

Weissella<-subset(top_lepto_genera, Genus == "Weissella")
Staphylococcus<-subset(top_lepto_genera, Genus == "Staphylococcus")
Fructobacillus<-subset(top_lepto_genera, Genus == "Fructobacillus")
Corynebacterium<-subset(top_lepto_genera, Genus == "Corynebacterium")
Streptococcus<-subset(top_lepto_genera, Genus == "Streptococcus")
Actinomyces<-subset(top_lepto_genera, Genus == "Actinomyces")
Gemella<-subset(top_lepto_genera, Genus == "Gemella")


###### GAMS #####

##### Weissella #####
Weissella_gam <- mgcv::gam(Abundance~
                             treatment+
                             s(time_point, by = treatment, bs = "cc", k=6) +
                             s(Seq_depth_scaled)+ day,
                           data=Weissella,
                           family = gaussian)

print(summary(Weissella_gam)) # no treatment effect (great)
gam.check(Weissella_gam)
plot(Weissella_gam)

Weissella_plot<-tidymv::plot_smooths(
  model = Weissella_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","dashed"))+
  ylab("Weissella") + xlab("day time (24 hours)") +
  scale_x_continuous(breaks = round(seq(min(as.numeric(Weissella$time_point)), max(as.numeric(Weissella$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))

####### Staphylococcus ########

Staphylococcus_gam <- mgcv::gam(Abundance~
                                  treatment+
                                  s(time_point, by = treatment, bs = "cc", k=6) + # specifically interested in this variable
                                  s(Seq_depth_scaled)+
                                  day,
                                data=Staphylococcus,
                                family = gaussian)

print(summary(Staphylococcus_gam)) #strong treatment effect (fewer in Staphylococcus in treatment)gam.check(Staphylococcus_gam) # plus day effect (random)
gam.check(Staphylococcus_gam)
plot(Staphylococcus_gam)


Staphylococcus_plot<-tidymv::plot_smooths(
  model = Staphylococcus_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","solid"))+
  ylab("Staphylococcus") + xlab("day time (24 hours)") +
  scale_x_continuous(breaks = round(seq(min(as.numeric(Weissella$time_point)), max(as.numeric(Weissella$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))

###########  Fructobacillus    ##############

Fructobacillus_gam <- mgcv::gam(Abundance~
                                  treatment+
                                  s(time_point, by = treatment, bs = "cc", k=6) + # specifically interested in this variable
                                  s(Seq_depth_scaled)+
                                  day,
                                data=Fructobacillus,
                                family = gaussian)

print(summary(Fructobacillus_gam)) # no effect of treatment but day (random effectgam.check(Fructobacillus_gam)
gam.check(Fructobacillus_gam)
plot(Fructobacillus_gam)


Fructobacillus_plot<-tidymv::plot_smooths(
  model = Fructobacillus_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","solid"))+
  ylab("Fructobacillus") + xlab("time points") +
  scale_x_continuous(breaks = round(seq(min(Weissella$time_point), max(Weissella$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))

###########  Corynebacterium    ##############

Corynebacterium_gam <- mgcv::gam(Abundance~
                                   treatment+
                                   s(time_point, by = treatment, bs = "cc", k=6) + # specifically interested in this variable
                                   s(Seq_depth_scaled)+
                                   day,
                                 data=Corynebacterium,
                                 family = gaussian)

print(summary(Corynebacterium_gam)) # some effect of treatment (more Corynebacterium in treatment)
gam.check(Corynebacterium_gam)
plot(Corynebacterium_gam)


Corynebacterium_plot<-tidymv::plot_smooths(
  model = Corynebacterium_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","dashed"))+
  ylab("Corynebacterium") + xlab("day time (24hrs)") +
  scale_x_continuous(breaks = round(seq(min(Weissella$time_point), max(Weissella$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))

###########  Streptococcus    ##############

Streptococcus_gam <- mgcv::gam(Abundance~
                                 treatment+
                                 s(time_point, by = treatment, bs = "cc", k=6) + # specifically interested in this variable
                                 s(Seq_depth_scaled)+
                                 day,
                               data=Streptococcus,
                               family = gaussian)

print(summary(Streptococcus_gam)) # no effect of treatment (good)
gam.check(Streptococcus_gam)
plot(Streptococcus_gam)


Streptococcus_plot<-tidymv::plot_smooths(
  model = Streptococcus_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","solid"))+
  ylab("Streptococcus") + xlab("day time (24hrs)") +
  scale_x_continuous(breaks = round(seq(min(Weissella$time_point), max(Weissella$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))



###########  Actinomyces    ##############

Actinomyces_gam <- mgcv::gam(Abundance~
                               treatment+
                               s(time_point, by = treatment, bs = "cc", k=6) + # specifically interested in this variable
                               s(Seq_depth_scaled)+
                               day,
                             data=Actinomyces,
                             family = gaussian)

print(summary(Actinomyces_gam)) ## no effect of treatment on abundance (good)
gam.check(Actinomyces_gam)
plot(Actinomyces_gam)


Actinomyces_plot<-tidymv::plot_smooths(
  model = Actinomyces_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","dashed"))+
  ylab("Actinomyces") + xlab("day time (24hrs)") +
  scale_x_continuous(breaks = round(seq(min(Actinomyces$time_point), max(Actinomyces$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))

###########  Gemella    ##############

Gemella_gam <- mgcv::gam(Abundance~
                           treatment+
                           s(time_point, by = treatment, bs = "cc", k=6) + # specifically interested in this variable
                           s(Seq_depth_scaled, bs="cr")+
                           day,
                         data=Gemella,
                         family = gaussian)

print(summary(Gemella_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(Gemella_gam) # and day effect (random)
plot(Gemella_gam)

ddply(Gemella, aes("treatment"), sum=sum(Abundance))

Gemella_plot<-tidymv::plot_smooths(
  model = Gemella_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","dashed"))+
  ylab("Gemella") + xlab("day time (24hrs)") +
  scale_x_continuous(breaks = round(seq(min(Actinomyces$time_point), max(Actinomyces$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))

# main text figure 

ggarrange(Weissella_plot,  
          Actinomyces_plot, 
          Gemella_plot,
          Corynebacterium_plot, 
          Staphylococcus_plot,
          Streptococcus_plot, 
          Fructobacillus_plot,
          ncol = 2, nrow=4, common.legend=TRUE, legend=FALSE, 
          labels = c("A", "B", "C", "D", "E", "F", "G"), 
          align="hv")

# supplementary figure 

ggarrange(Staphylococcus_plot, Fructobacillus_plot, 
          ncol = 1, common.legend=TRUE, legend="bottom", 
          labels = c("A", "B", "C"))


######### ph as functional gut activity proxy ####

#### visualise via box plot 

pH_boxplot<-ggplot(Weissella, aes(as.factor(time_point), pH_value, fill=treatment))+ geom_boxplot(outlier.alpha = 0.0)+
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_classic()+
  geom_point(shape=21, position=position_dodge(.75))+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))+
  xlab("day time (24 hours)") + ylab("faecal pH")+
  scale_x_discrete(breaks = round(seq(min(Actinomyces$time_point), max(Actinomyces$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))


#### pH GAMs ####

pH_gam <- mgcv::gam(pH_value~
                      treatment+
                      s(time_point, by = treatment, bs = "cc", k=6) + 
                      day,
                    data=Gemella,
                    family = gaussian, REML=T)

print(summary(pH_gam)) ## strong effect of treatment on abundance (implies: Gamella less frequent in control group)
gam.check(pH_gam) # and day effect (random)
#plot(Gemella_gam)

pH_gam_plot<-tidymv::plot_smooths(
  model = pH_gam,
  series = time_point,
  comparison = treatment, 
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","solid"))+
  ylab("faecal pH") + xlab("day time (24hrs)") +
  scale_x_continuous(breaks = round(seq(min(Actinomyces$time_point), max(Actinomyces$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))


ggarrange(SupFig5a, SupFig5b, 
          Fig3a, Shannon_plot, 
          PC1_boxplot, PC1_plot, 
          pH_boxplot, pH_gam_plot, ncol = 2, nrow=4, 
          common.legend=TRUE, legend="bottom", align="v",
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"))


##### functional link between pH and microbial abundances ####

anova(lm(Observed.unrare ~ pH_value * treatment+day, data=meta))

Observed_pH<-ggplot(meta, aes(pH_value, Observed.unrare, fill=treatment))+ 
  geom_point(shape=21)+
  #geom_smooth(method=lm, aes(colour=treatment), lty="dashed") + 
  #scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position = "none")+
  ylab("Observed ASVs") + xlab("fecal pH")

anova(lm(Shannon.unrare ~ pH_value * treatment+day, data=meta))

Shannon_pH<-ggplot(meta, aes(pH_value, Shannon.unrare, fill=treatment))+ 
  geom_point(shape=21)+
  #geom_smooth(method=lm, aes(colour=treatment), lty="dashed") + 
  #scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position = "none")+
  ylab("Shannon") + xlab("fecal pH")

anova(lm(PC1 ~ pH_value * treatment+day, data=meta.PCoA))

PC1_pH<-ggplot(meta.PCoA, aes(pH_value, PC1, fill=treatment))+ 
  geom_point(shape=21)+
  #geom_smooth(method=lm, aes(colour=treatment), lty="dashed") + 
  #scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position = "none")+
  ylab("PC1") + xlab("fecal pH")

anova(lm(PC1_unweighted ~ pH_value * treatment+day, data=meta.PCoA))

PC1_unweighted_pH<-ggplot(meta.PCoA, aes(pH_value, PC1_unweighted, fill=treatment))+ 
  geom_point(shape=21)+
  geom_smooth(method=lm, aes(colour=treatment), lty="dashed") + 
  scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position = "none")+
  ylab("PC1") + xlab("fecal pH")


anova(lm(Abundance ~ pH_value * treatment+day, data=Weissella))

W_pH<-ggplot(Weissella, aes(pH_value, Abundance, fill=treatment))+ 
  geom_point(shape=21)+
  geom_smooth(method=lm, aes(colour=treatment)) + 
  scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position = "none")+
  ylab("Weissella") + xlab("fecal pH")


anova(lm(Abundance ~ pH_value * treatment+day, data=Actinomyces))

A_pH<-ggplot(Actinomyces, aes(pH_value, Abundance, fill =treatment))+ 
  geom_point(shape=21)+
  geom_smooth(method=lm, aes(colour=treatment), lty="dashed") + 
  scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position="none")+
  ylab("Actinomyces") + xlab("fecal pH") 


anova(lm(Abundance ~ pH_value * treatment+day, data=Streptococcus))

S_pH<-ggplot(Streptococcus, aes(pH_value, Abundance, fill=treatment))+ 
  geom_point(shape=21)+
  #geom_smooth(method=lm, colour="black") + 
  scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position="none")+
  ylab("Streptococcus") + xlab("fecal pH") 



#### non significant interactions #####
anova(lm(Abundance ~ pH_value * treatment + day, data=Corynebacterium))
anova(lm(Abundance ~ pH_value * treatment + day, data=Fructobacillus))
anova(lm(Abundance ~ pH_value * treatment + day, data=Gemella))
anova(lm(Abundance ~ pH_value * treatment + day, data=Staphylococcus))


ggarrange(Observed_pH, Shannon_pH, PC1_unweighted_pH, 
          W_pH, A_pH, S_pH,  
          nrow=2, ncol = 3, align="hv",
          labels = c("A", "B", "C", 
                     "D", "E", "F"))


