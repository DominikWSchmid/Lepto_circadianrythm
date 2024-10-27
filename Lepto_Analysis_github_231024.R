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
ps = readRDS("~/Dropbox/CellarBats/working_data/CellarBats_ps_filtered_100reads.rds")

#### remove samples with fewer than 500 reads ####

sample_data(ps)$sample_depth<-sample_sums(ps)

ps<-subset_samples(ps, sample_data(ps)$sample_depth>=500)

### some information about the phyloseq object and make meta df

microbiome::summarize_phyloseq(ps)

sd(sample_sums(ps)) ## on average we still have more than 30 000 reads from our samples 

meta<- as(sample_data(ps), "data.frame")


##### sample size plot ####

N<-plyr::ddply(meta, c("treatment", "time_point"), summarise, n=length(time_point))

ggplot(N, aes(x =fct_relevel(time_point, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"),
              y =n))+
  geom_bar(stat="identity", fill="grey75", colour="black")+
  facet_wrap(~treatment)+
  ylim(0,15)+
  theme_classic()+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=14), 
        axis.title.y = element_text(colour = 'black', face="bold", size=14))+
  xlab("Time of Day (24hrs)")+
  ylab("Sample Size")


######## now for the microbiome ########

ps_class<-tax_glom(ps, taxrank = "Class")

sample_data(ps)<-sample_data(ps)[,c("sample.ID")]
ps_class <- microbiome::transform(ps_class, "compositional") 
class_per_sample<-psmelt(ps_class)

plyr::ddply(class_per_sample, c("Class"), summarise, mean=round(mean(Abundance),4))


ps_genus<-tax_glom(ps, taxrank = "Genus")

ps_genus <- microbiome::transform(ps_genus, "compositional") 
genus_per_sample<-psmelt(ps_genus)

df_genus<-plyr::ddply(genus_per_sample, c("Genus"), summarise, mean=round(mean(Abundance),4))

##### find most common classes ####

top<-as_tibble(microbiomeutilities::get_group_abundances(ps_class, level="Class", group="Class", transform = "compositional") %>%arrange (-mean_abundance))
print(top,n=15)

##### classes that are more prominent than 1%
top_class<-c("Bacilli",  "Actinobacteria", "Gammaproteobacteria", "Alphaproteobacteria", "Clostridia")


##### now top genera
ps_genus<-tax_glom(ps, taxrank = "Genus")

top<-as_tibble(microbiomeutilities::get_group_abundances(ps_genus, level="Genus", group="Genus", transform = "compositional") %>%arrange (-mean_abundance))
print(top,n=15)

top_genus<-c("Weissella",  "Staphylococcus", "Fructobacillus", "Corynebacterium", "Streptococcus", "Actinomyces", "Gemella")


##### relative abundance plots by time point and treatment ######
ps = readRDS("~/Dropbox/CellarBats/working_data/CellarBats_ps_filtered_100reads.rds")
sample_data(ps)$sample_depth<-sample_sums(ps)
ps<-subset_samples(ps, sample_data(ps)$sample_depth>=500)

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
        legend.title = element_text(colour = 'black', face="bold", size=12))+labs(fill ="Bacterial Genus")+
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
        legend.title = element_text(colour = 'black', face="bold", size=12))+labs(fill ="Bacterial Genus")+
  ylab("Relative Abundance")+
  theme(axis.title.x = element_blank())+scale_y_continuous(expand = c(0,0))


ggarrange(
  summary_figXB, summary_figXA, 
  labels = c("C", "D"),
  ncol=2, align = c("h"),
  common.legend = TRUE, legend = "bottom"
)


#####  alpha diversity metrices and some stats ######
#####################################################

#first we rarefy 

ps.rare <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 123)

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
cor.test(meta$Shannon.unrare, meta$Shannon)

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  hjustvar = c(1.1),
  vjustvar = c(2))

SupFig3a<-ggplot(meta, aes(Observed.unrare, Observed))+geom_point()+geom_smooth(method="lm")+theme_bw()+
ylab("rarefied Observed ASVs")+xlab("unrarefied Observed ASVs")+
theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), axis.title.y = element_text(colour = 'black', face="bold", size=12))+
geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=81.7, p<0.001"))

SupFig3b<-ggplot(meta, aes(Shannon.unrare, Shannon))+geom_point()+geom_smooth(method="lm")+theme_bw()+
  ylab("rarefied Shannon")+xlab("unrarefied Shannon Index")+
  theme(axis.title.x = element_text(colour = 'black', face="bold", size=12), axis.title.y = element_text(colour = 'black', face="bold", size=12))+
  geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label="R2=99.7, p<0.001"))

ggarrange(SupFig3a, SupFig3b, ncol = 2, nrow=1, common.legend=TRUE, legend="right", labels=c("A","B"))

#### because of these results we are happy to proceed on the unrarefied data for alpha diversity ####


##### calculation of alpha for all 187 samples ####
otu.table <- as.data.frame(otu_table(ps))
meta <- as(sample_data(ps), "data.frame")
df.pd <- pd(t(otu.table), phy_tree(ps), include.root=T)
meta$Phyogenetic_diversity <- df.pd$PD

alpha = estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")) #generate df with various diversity metices 

#add indicies to meta data (like we did with PD)
meta$Shannon.unrare <- alpha$Shannon 
meta$Chao1.unrare <- alpha$Chao1 
meta$Observed.unrare <- alpha$Observed


############## Alpha GAMS #############
######### OBSERVED ASVs #########

SupFig5a<-ggplot(meta, aes(fct_relevel(time_point, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                 (Observed.unrare), fill=treatment)) + 
  geom_boxplot( outlier.alpha = 0) + geom_point(shape=21, position=position_dodge(.8))+ theme_classic() +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547")) + xlab("day time (24hours)") + ylab("Observed ASVs")+ylim(0,70)+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position= "none")


##### gam observed ###
meta.gam<-meta #### because time point needs to be numerical for gams
meta.gam$time_point<-as.numeric(meta.gam$time_point)
meta.gam$Seq_depth_scaled<-as.numeric(rescale(sample_sums(ps),, to =c(-1,1)))

Observed_gam <- mgcv::gam(Observed.unrare~
                            treatment+ day + s(Seq_depth_scaled) +
                            s(time_point, by = treatment, bs = "cc"),
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
  theme_bw() + scale_linetype_manual(values=c("solid","dashed"))+
  ylab("Observed ASVs") + xlab("day time (24 hours)") + ylim(0,70)+
  scale_x_continuous(breaks = round(seq(min(as.numeric(meta.gam$time_point)), max(as.numeric(meta.gam$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position="none")


#### replicate on rarefied data ###
##### gam observed ###
meta.gam<-meta #### because time point needs to be numerical for gams
meta.gam$time_point<-as.numeric(meta.gam$time_point)
meta.gam$Seq_depth_scaled<-as.numeric(rescale(sample_sums(ps),, to =c(-1,1)))

Observed_gam_rarefied <- mgcv::gam(Observed~
                            treatment+ day + s(Seq_depth_scaled) +
                            s(time_point, by = treatment, bs = "cc"),
                          data=meta.gam,
                          family = gaussian, REML=T)

print(summary(Observed_gam_rarefied)) # no treatment effect (great)
gam.check(Observed_gam_rarefied)
plot(Observed_gam_rarefied)

### comparison figure ####
SupFig4a<-tidymv::plot_smooths(
  model = Observed_gam_rarefied,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","dashed"))+
  ylab("Observed ASVs") + xlab("Time of day (24 hours)") +
  scale_x_continuous(breaks = round(seq(min(as.numeric(meta.gam$time_point)), max(as.numeric(meta.gam$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(legend.position="none")


SupFig4b<-tidymv::plot_smooths(
  model = Observed_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","dashed"))+
  ylab("Observed ASVs") + xlab("day time (24 hours)") +
  scale_x_continuous(breaks = round(seq(min(as.numeric(meta.gam$time_point)), max(as.numeric(meta.gam$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(legend.position="none")

ggarrange(SupFig4a, SupFig4b, ncol = 2, nrow=1, common.legend=TRUE, legend="right", labels=c("A","B"))


###### SHANNON ######

Fig3a<-ggplot(meta, aes(fct_relevel(time_point, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                         (Shannon.unrare), fill=treatment)) + 
  geom_boxplot( outlier.alpha = 0) + geom_point(shape=21, position=position_dodge(.75))+ theme_classic() +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547")) + xlab("day time (24hours)") + ylab("Shannon")+ ylim(0,2.5)+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position="none")

### shannon gam ###
Shannon_gam <- mgcv::gam(Shannon.unrare~
                           treatment+ day + s(Seq_depth_scaled) +
                           s(time_point, by = treatment, bs = "cc"),
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
  xlab("day time (24 hours)") + ylim(0,2.5)+
  scale_x_continuous(breaks = round(seq(min(as.numeric(meta.gam$time_point)), max(as.numeric(meta.gam$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none")


######## beta diversity index ######

ps.rare <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 123)

# calculate weighted Unifrag distance using the phyloseq package
dist_wUni = phyloseq::distance(ps.rare, method="weighted unifrac") ### abundance matters + phylogeny matters 
# calculate unweighted distance using the phyloseq package
dist_Uni = phyloseq::distance(ps.rare, method="unifrac") ### presence/absence + phylogeny matters

meta <- as(sample_data(ps.rare), "data.frame")

meta$Seq_depth_scaled<-as.numeric(rescale(sample_sums(ps.rare), to =c(-1,1)))

adonis2(dist_Uni ~ treatment*time_point + day , data = meta, method="unifrac", by="margin")
adonis2(dist_wUni ~ treatment*time_point + day , data = meta, method="weighted unifrac", by="margin")


####### pairwise permanova #####
meta$timeBYtreatment<- paste0(meta$time_point,meta$treatment)

result_list_weighted<-pairwise.adonis2(phyloseq::distance(ps.rare, method = "weighted unifrac") ~ timeBYtreatment, 
                                       data = meta)
result_list_unweighted<-pairwise.adonis2(phyloseq::distance(ps.rare, method = "unifrac") ~ timeBYtreatment, 
                                         data = meta)

#####
control<-subset_samples(ps.rare, treatment=="C")
meta_c<-sample_data(control)

dist_Uni_c = phyloseq::distance(control, method="unifrac")
ordination_Uni_c<-ordinate(control, method="PCoA", distance=dist_Uni_c)
#####
treatment<-subset_samples(ps.rare, treatment=="T")
meta_t<-sample_data(treatment)

dist_Uni_t = phyloseq::distance(treatment, method="unifrac") ### presence/absence + phylogeny matters
ordination_Uni_t<-ordinate(treatment, method="PCoA", distance=dist_Uni_t)
#####
meta.control <- as(sample_data(control), "data.frame")
dist_wUni_c = phyloseq::distance(control, method="wunifrac") ### abundance matters + phylogeny matters
ordination_wUni_c<-ordinate(control, method="PCoA", distance=dist_wUni_c)
#####
meta.treatment <- as(sample_data(treatment), "data.frame")
dist_wUni_t = phyloseq::distance(treatment, method="wunifrac") ### presence/absence + phylogeny matters
ordination_wUni_t<-ordinate(treatment, method="PCoA", distance=dist_wUni_t)
#####
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
                                   PC1, fill=treatment)) + 
  geom_boxplot(outlier.alpha = 0) + geom_point(shape=21, position=position_dodge(.75))+ theme_classic() +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547")) + xlab("day time (24hours)") + ylab("PC1 scores")+ylim(-0.11,0.31)+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))


##### gam PC1 ###
meta.PCoA$time_point<-as.numeric(meta.PCoA$time_point)
#meta.PCoA$Seq_depth_scaled<-as.numeric(rescale(sample_sums(ps), to =c(-1,1))) #since we rarefied 

PC1_gam <- mgcv::gam((PC1)~
                       treatment+
                       s(time_point, by = treatment, bs = "cc") + # specifically interested in this variable
                       #s(Seq_depth_scaled)+
                       day,
                     data=meta.PCoA,
                     family = gaussian, REML=T)

print(summary(PC1_gam)) # no treatment effect (great)
mgcv::gam.check(PC1_gam)
plot(PC1_gam)


PC1_plot<-tidymv::plot_smooths(
  model = PC1_gam,
  series = time_point,
  comparison = treatment,
  exclude_terms = c("day", "Seq_depth_scaled")
) + scale_colour_manual(values = c( "#5C80BC", "#E8C547")) + 
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ 
  theme_bw() + scale_linetype_manual(values=c("solid","solid"))+
  ylab("PC1 scores") + xlab("day time (24 hours)") +ylim(-0.11, 0.31)+
  scale_x_continuous(breaks = round(seq(min((meta.PCoA$time_point)), max((meta.PCoA$time_point)), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.title = element_text(colour="black",size=12,face="bold"))


##### supplementary figure ####

# constrained ordination, weighted unifrac

control<-subset_samples(ps.rare, treatment=="C")

dist_wUni = phyloseq::distance(control, method="weighted unifrac")
otutable<-data.frame(t(control@otu_table@.Data))
metadata <- data.frame(sample_data(control))
names(metadata)

metadata$time_point2<-ifelse(metadata$time_point == "1", "1-3", NA)
metadata$time_point2<-ifelse(metadata$time_point == "3", "1-3", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "5", "5-7", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "7", "5-7", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "9", "9-11", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "11", "9-11", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "13", "13-15", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "15", "13-15", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "17", "17-19", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "19", "17-19", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "21", "21-23", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "23", "21-23", metadata$time_point2)

time_point<-factor(metadata$time_point2)



final_model<-capscale(dist_wUni ~ 
                        time_point,
                      env = metadata, 
                      comm = otutable) 
final_model$CCA$rank

anova_wunifrac<-anova.cca(final_model, by="terms")
anova_wunifrac

# extract data from model

final_model_df<-scores(final_model)

# extract CAP scores

vectors_df<-data.frame(final_model_df$sites)
vectors_df$sample.ID<-row.names(vectors_df)


# merge with info on dominant family

sample_metadata<-metadata[,c("sample.ID", "time_point2")]
vectors_df<-merge(vectors_df, sample_metadata, by = "sample.ID")


##### change order

vectors_df$time_point2<-factor(vectors_df$time_point2, levels = c("1-3", "5-7", "9-11", "13-15", "17-19", "21-23"))

vectors_control<-vectors_df

## plot

beta_control<-ggplot(vectors_control, aes(x = CAP1, y = CAP2))+
  
  stat_ellipse(geom = "polygon", aes(fill = time_point2), level = 0.9, alpha = 0.3, size = 0.5)+
  geom_point(aes(fill =time_point2), pch = 21, size = 3, alpha = 1, stroke = 1, col = "black")+
  
  theme_bw()+
  scale_fill_brewer(palette ="Paired")+
  scale_color_brewer(palette ="Paired")+
  
  theme_light(base_size = 12)+
  
  labs(fill = "Time point", col = "Time point")

##################################################################################


# treatment group 

treatment<-subset_samples(ps.rare, treatment=="T")

dist_wUni = phyloseq::distance(treatment, method="weighted unifrac")
otutable<-data.frame(t(treatment@otu_table@.Data))
metadata <- data.frame(sample_data(treatment))
names(metadata)

metadata$time_point2<-ifelse(metadata$time_point == "1", "1-3", NA)
metadata$time_point2<-ifelse(metadata$time_point == "3", "1-3", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "5", "5-7", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "7", "5-7", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "9", "9-11", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "11", "9-11", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "13", "13-15", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "15", "13-15", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "17", "17-19", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "19", "17-19", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "21", "21-23", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "23", "21-23", metadata$time_point2)



time_point<-factor(metadata$time_point2)



final_model_t<-capscale(dist_wUni ~ 
                          time_point,
                        env = metadata, 
                        comm = otutable) 


anova_wunifrac<-anova.cca(final_model_t, by="terms")
anova_wunifrac

# extract data from model

final_model_df<-scores(final_model_t)

# extract CAP scores

vectors_df<-data.frame(final_model_df$sites)
vectors_df$sample.ID<-row.names(vectors_df)


# merge with info on dominant family

sample_metadata<-metadata[,c("sample.ID", "time_point2")]
vectors_df<-merge(vectors_df, sample_metadata, by = "sample.ID")


##### change order

vectors_df$time_point2<-factor(vectors_df$time_point2, levels = c("1-3", "5-7", "9-11", "13-15", "17-19", "21-23"))


vectors_treatment <-vectors_df



## plot

beta_treatment<- ggplot(vectors_treatment, aes(x = CAP1, y = CAP2))+
  
  stat_ellipse(geom = "polygon", aes(fill = time_point2), level = 0.9, alpha = 0.3, size = 0.5)+
  geom_point(aes(fill =time_point2), pch = 21, size = 3, alpha = 1, stroke = 1, col = "black")+
  
  theme_bw()+
  scale_fill_brewer(palette ="Paired")+
  scale_color_brewer(palette ="Paired")+
  
  theme_light(base_size = 12)+
  
  labs(fill = "Time point", col = "Time point")

ggarrange(beta_control, beta_treatment, common.legend = T, 
          labels = c("Control", "Treatment"), legend = "right")



#### constrained ordination, unweighted unifrac#####

control<-subset_samples(ps.rare, treatment=="C")

dist_Uni = phyloseq::distance(control, method="unweighted unifrac")
otutable<-data.frame(t(control@otu_table@.Data))
metadata <- data.frame(sample_data(control))
names(metadata)

metadata$time_point2<-ifelse(metadata$time_point == "1", "1-3", NA)
metadata$time_point2<-ifelse(metadata$time_point == "3", "1-3", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "5", "5-7", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "7", "5-7", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "9", "9-11", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "11", "9-11", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "13", "13-15", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "15", "13-15", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "17", "17-19", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "19", "17-19", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "21", "21-23", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "23", "21-23", metadata$time_point2)

time_point<-factor(metadata$time_point2)



final_model<-capscale(dist_Uni ~ 
                        time_point,
                      env = metadata, 
                      comm = otutable) 
final_model$CCA$rank

anova_wunifrac<-anova.cca(final_model, by="terms")
anova_wunifrac

# extract data from model

final_model_df<-scores(final_model)

# extract CAP scores

vectors_df<-data.frame(final_model_df$sites)
vectors_df$sample.ID<-row.names(vectors_df)


# merge with info on dominant family

sample_metadata<-metadata[,c("sample.ID", "time_point2")]
vectors_df<-merge(vectors_df, sample_metadata, by = "sample.ID")


##### change order

vectors_df$time_point2<-factor(vectors_df$time_point2, levels = c("1-3", "5-7", "9-11", "13-15", "17-19", "21-23"))

vectors_control<-vectors_df

## plot

beta_control_unweighted<-ggplot(vectors_control, aes(x = CAP1, y = CAP2))+
  
  stat_ellipse(geom = "polygon", aes(fill = time_point2), level = 0.9, alpha = 0.3, size = 0.5)+
  geom_point(aes(fill =time_point2), pch = 21, size = 3, alpha = 1, stroke = 1, col = "black")+
  
  theme_bw()+
  scale_fill_brewer(palette ="Paired")+
  scale_color_brewer(palette ="Paired")+
  
  theme_light(base_size = 12)+
  
  labs(fill = "Time point", col = "Time point")

##################################################################################


# treatment group 

treatment<-subset_samples(ps.rare, treatment=="T")

dist_Uni = phyloseq::distance(treatment, method="unweighted unifrac")
otutable<-data.frame(t(treatment@otu_table@.Data))
metadata <- data.frame(sample_data(treatment))
names(metadata)

metadata$time_point2<-ifelse(metadata$time_point == "1", "1-3", NA)
metadata$time_point2<-ifelse(metadata$time_point == "3", "1-3", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "5", "5-7", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "7", "5-7", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "9", "9-11", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "11", "9-11", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "13", "13-15", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "15", "13-15", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "17", "17-19", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "19", "17-19", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "21", "21-23", metadata$time_point2)
metadata$time_point2<-ifelse(metadata$time_point == "23", "21-23", metadata$time_point2)



time_point<-factor(metadata$time_point2)



final_model_t<-capscale(dist_Uni ~ 
                          time_point,
                        env = metadata, 
                        comm = otutable) 


anova_wunifrac<-anova.cca(final_model_t, by="terms")
anova_wunifrac

# extract data from model

final_model_df<-scores(final_model_t)

# extract CAP scores

vectors_df<-data.frame(final_model_df$sites)
vectors_df$sample.ID<-row.names(vectors_df)


# merge with info on dominant family

sample_metadata<-metadata[,c("sample.ID", "time_point2")]
vectors_df<-merge(vectors_df, sample_metadata, by = "sample.ID")


##### change order

vectors_df$time_point2<-factor(vectors_df$time_point2, levels = c("1-3", "5-7", "9-11", "13-15", "17-19", "21-23"))


vectors_treatment <-vectors_df



## plot

beta_treatment_unweighted<- ggplot(vectors_treatment, aes(x = CAP1, y = CAP2))+
  
  stat_ellipse(geom = "polygon", aes(fill = time_point2), level = 0.9, alpha = 0.3, size = 0.5)+
  geom_point(aes(fill =time_point2), pch = 21, size = 3, alpha = 1, stroke = 1, col = "black")+
  
  theme_bw()+
  scale_fill_brewer(palette ="Paired")+
  scale_color_brewer(palette ="Paired")+
  
  theme_light(base_size = 12)+
  
  labs(fill = "Time point", col = "Time point")



SuppFig5<-ggarrange(beta_control_unweighted, beta_treatment_unweighted, beta_control, beta_treatment, 
                    common.legend = T, ncol = 2, nrow=2, legend = "right", labels=c("A","B", "C", "D"))


ggarrange(SupFig3a, SupFig3b, ncol = 2, nrow=1, common.legend=TRUE, legend="right", labels=c("A","B"))





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
  theme_bw() + scale_linetype_manual(values=c("solid","solid"))+
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
  xlab("day time (24 hours)") + ylab("faecal pH")+ylim(5.0,7.5)+
  scale_x_discrete(breaks = round(seq(min(Actinomyces$time_point), max(Actinomyces$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))


#### pH GAMs ####

pH_gam <- mgcv::gam(pH_value~
                      treatment+
                      s(time_point, by = treatment, bs = "cc") + 
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
  ylab("faecal pH") + xlab("day time (24hrs)") + ylim(5.0,7.5)+
  scale_x_continuous(breaks = round(seq(min(Actinomyces$time_point), max(Actinomyces$time_point), by = 1),1), labels=c('1', '3', '5', "7", "9", "11", "13", "15", "17", "19", "21", "23"))+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.title = element_text(colour="black",size=12,face="bold"))


ggarrange(SupFig5a, SupFig5b, 
          Fig3a, Shannon_plot, 
          PC1_boxplot, PC1_plot, 
          pH_boxplot, pH_gam_plot, ncol = 2, nrow=4, 
          common.legend=TRUE, legend="bottom", align="v",
          labels = c("A", "", "B", "", "C", "", "D", ""))


##### functional link between pH and microbial abundances ####

summary(lm(Observed.unrare ~ pH_value * treatment+day, data=meta.gam))

Observed_pH<-ggplot(meta.gam, aes(pH_value, Observed.unrare, fill=treatment))+ 
  geom_point(shape=21)+
  #geom_smooth(method=lm, aes(colour=treatment), lty="dashed") + 
  #scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position = "none")+
  ylab("Observed ASVs") + xlab("fecal pH")

summary(lm(Shannon.unrare ~ pH_value * treatment+day, data=meta.gam))

Shannon_pH<-ggplot(meta.gam, aes(pH_value, Shannon.unrare))+ 
  geom_point(shape=21, fill="black")+
  geom_smooth(method=lm, colour="black", lty="dashed") + 
  theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position = "none")+
  ylab("Shannon") + xlab("fecal pH")

summary(lm(PC1 ~ pH_value * treatment+day, data=meta.PCoA))

PC1_pH<-ggplot(meta.PCoA, aes(pH_value, PC1))+ 
  geom_point(shape=21, fill="black")+
  geom_smooth(method=lm, colour="black", lty="dashed") + 
  theme_bw()+
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position = "none")+
  ylab("PC1") + xlab("fecal pH")

summary(lm(PC1_unweighted ~ pH_value * treatment+day, data=meta.PCoA))

PC1_unweighted_pH<-ggplot(meta.PCoA, aes(pH_value, PC1_unweighted, fill=treatment))+ 
  geom_point(shape=21, fill="black")+
  #geom_smooth(method=lm, colour="black", lty="dashed") + 
  theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_blank(),
        legend.position = "none")+
  ylab("PC1") + xlab("fecal pH")


summary(lm(Abundance ~ pH_value * treatment+day, data=Weissella))

W_pH<-ggplot(Weissella, aes(pH_value, Abundance, fill=treatment))+ 
  geom_point(shape=21)+
  geom_smooth(method=lm, aes(colour=treatment)) + 
  scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position = "none")+
  ylab("        Relative abundance 
       of Weissella") + xlab("fecal pH")


summary(lm(Abundance ~ pH_value * treatment+day, data=Actinomyces))

A_pH<-ggplot(Actinomyces, aes(pH_value, Abundance, fill=treatment))+ 
  geom_point(shape=21)+
  geom_smooth(method=lm, aes(colour=treatment)) + 
  scale_colour_manual(values = c( "#5C80BC", "#E8C547")) +
  scale_fill_manual(values = c( "#5C80BC", "#E8C547"))+ theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position="none")+
  ylab("        Relative abundance 
       of Actinomyces") + xlab("fecal pH") 


summary(lm(Abundance ~ pH_value * treatment+day, data=Streptococcus))

S_pH<-ggplot(Streptococcus, aes(pH_value, Abundance))+ 
  geom_point(shape=21, fill="black")+
  geom_smooth(method=lm, colour="black") + 
  theme_bw()+
  theme(axis.title.y = element_text(colour="black",size=12,face="bold"), 
        axis.title.x = element_text(colour="black",size=12,face="bold"),
        legend.position="none")+
  ylab("        Relative abundance 
       of Streptococcus") + xlab("fecal pH") 



#### non significant interactions #####
summary(lm(Abundance ~ pH_value * treatment + day, data=Corynebacterium))
summary(lm(Abundance ~ pH_value * treatment + day, data=Fructobacillus))
summary(lm(Abundance ~ pH_value * treatment + day, data=Gemella))
summary(lm(Abundance ~ pH_value * treatment + day, data=Staphylococcus))


ggarrange(Observed_pH, Shannon_pH, PC1_unweighted_pH, 
          W_pH, A_pH, S_pH,  
          nrow=2, ncol = 3, align="hv",
          labels = c("A", "B", "C", 
                     "D", "E", "F"))


##### comp plot for each sample ordered by time of day

ps_treat_time_genus<-tax_glom(ps, taxrank = "Genus")
ps_treat_time_genus <- microbiome::transform(ps_treat_time_genus, "compositional")
genus_by_treat_time<-psmelt(ps_treat_time_genus)

genus_by_treat_time$Genus_plot<-as.factor(ifelse(genus_by_treat_time$Genus %in% top_genus, genus_by_treat_time$Genus, "Other"))

genus_by_treat_time$Genus_plot<-factor(genus_by_treat_time$Genus_plot, 
                                       levels = c("Weissella",  "Staphylococcus", "Fructobacillus", "Corynebacterium", "Streptococcus", "Actinomyces", "Gemella", "Other"))

genus_by_treat_time$Sample<-as.factor(genus_by_treat_time$Sample)

treatment<-subset(genus_by_treat_time, treatment=="T")

Supp_treatment<-ggplot(treatment, aes(x =fct_relevel(Sample, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                               y = Abundance, fill = Genus_plot))+
  geom_col(size = 0, position = "fill")+ 
  theme_classic()+scale_fill_manual(values=c( "#A97092", "#C28CAE", "#DEC4A1", "#EBBB54","#EDEB7A", "#A3D186", "#76B19B","grey65"))+
  theme(axis.text.x = element_text(size = 10, angle=90, vjust=0.5), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.title = element_text(colour = 'black', face="bold", size=12))+labs(fill ="Bacterial Genus")+
  ylab("Relative Abundance")+xlab("day time (24hours)")+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank())+scale_y_continuous(expand = c(0,0))


control<-subset(genus_by_treat_time, treatment=="C")

Supp_control<-ggplot(control, aes(x =fct_relevel(Sample, "1","3", "5", "7","9", "11", "13", "15", "17", "19", "21", "23"), 
                                      y = Abundance, fill = Genus_plot))+
  geom_col(size = 0, position = "fill")+ 
  theme_classic()+scale_fill_manual(values=c( "#A97092", "#C28CAE", "#DEC4A1", "#EBBB54","#EDEB7A", "#A3D186", "#76B19B","grey65"))+
  theme(axis.text.x = element_text(size = 10, angle=90, vjust=0.5), 
        axis.title.y = element_text(colour = 'black', face="bold", size=12), 
        legend.title = element_text(colour = 'black', face="bold", size=12))+labs(fill ="Bacterial Genus")+
  ylab("Relative Abundance")+xlab("day time (24hours)")+
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank())+scale_y_continuous(expand = c(0,0))


ggarrange(
  Supp_control, Supp_treatment, 
  labels = c("A", "B"),
  nrow=2, align = c("h"),
  common.legend = TRUE, legend = "bottom"
)


##### SAMPLE SIZE PLOT ####

N_byday<-plyr::ddply(meta, c("day", "treatment", "time_point"), summarise, n=length(time_point))
N_byday$day<-ifelse(N_byday$day==3, 2, 1)
N_byday$collection_vs_quality<-"quality"

collection_dat = read.csv("~/Dropbox/CellarBats/working_data/Sample_collection2.0.csv", sep=";")

N_total_byday<-plyr::ddply(collection_dat, c("day", "treatment", "time_point"), summarise, n=length(time_point))
N_total_byday$collection_vs_quality<-"collection"

N<-rbind(N_byday, N_total_byday)

treat.labs <- c("Control - regular feeding", "Treatment - delayed feeding")
names(treat.labs) <- c("C", "T")

day.labs <- c("first day", "second day")
names(day.labs) <- c("1", "2")

ggplot(data=N, aes(x=time_point, fill=collection_vs_quality)) +
  geom_bar(aes(y=n), stat="identity", position ="identity", alpha=.3, color='lightblue4')+
  facet_wrap(~treatment+day, labeller = labeller(treatment = treat.labs, day = day.labs))+
  geom_hline(yintercept=5, linetype="dashed")+theme_bw()+
  theme(legend.position ="none", 
        axis.title.x = element_text(colour = 'black', face="bold", size=14), 
        axis.title.y = element_text(colour = 'black', face="bold", size=14))+
  scale_fill_manual(values=c("#74D3AE", "black"))+ylab("Sample Size")+xlab("Time of Day (24 hours)")

