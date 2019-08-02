#Italy Invert Figs
#The following code is for the figures presented in the manuscript. 
#For all other code use for analysis see ItalyInvertRmarkdownJune8.rmd or ItalyInvertRmarkdownJune8.html for precomputed output
#Joe Receveur Email: josephreceveur[at] gmail.com
#######
#Import
#######

library(vegan)
library(MASS)
library(ggplot2)
library(plyr)
library(dplyr)
library(magrittr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(randomForest)
library(knitr)
library(ape)
library(ggpubr)
library(multcompView)
library(Rmisc)
library(tiff)

set.seed(27)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7","#E12D00")
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
biom=import_biom("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvert2018WTax.biom",parseFunction= parse_taxonomy_greengenes)
#tax_table(biom) <- tax_table(biom)[,-c(5:10,14)]#remove dummy ranks

metadata=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvertMetadataWDiversity3.5.2019.csv",header = TRUE)
head(metadata)

tree=read_tree("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvert11.29.18Tree.nwk")

sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$id
physeq=merge_phyloseq(biom,sampdat,tree)
physeq=rarefy_even_depth(physeq, 3000, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
physeq




#######
#Figure 3A Weights by station
#######

Trtdata <- ddply(metadata, c("Sampling_station", "FFG","Taxon_name"), summarise,
                 N    = length(id),
                 meanWeight = mean(Mass),
                 sd   = sd(Mass),
                 se   = sd / sqrt(N)
)
Trtdata

Scrapers<-subset(metadata, FFG=="Scraper")
t.test( Mass~ Sampling_station, data=Scrapers)
#t-test,t=4.219,P = 0.0047

Shredders<-subset(metadata, FFG=="Shredder")
t.test( Mass~ Sampling_station, data=Shredders)
#t=1.81,P = 0.2095

MassPlot<-ggplot(metadata,aes(x=Sampling_station,y=Mass,fill=Sampling_station))+geom_boxplot()+ #Write in labels from posthoc
  scale_fill_manual(values=cbPalette)+xlab("")+ylab("Mass (mg)")+labs(fill="FFG")+guides(fill=FALSE)+facet_wrap(~FFG,scales = "free_y")+
  theme_bw(base_size = 10)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#+theme(axis.text.x = element_text(angle = 45, 
MassPlot

tiff("Figures/MassPlot.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
MassPlot
dev.off()

#######
#Figure 3B Phylum level stacked bar graph 
#######
theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
PhylumAll=tax_glom(GPr, "Phylum")

PhylumLevel = filter_taxa(PhylumAll, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%

df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "FFG","Sampling_station"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata


RelAbuAllSamples <- ddply(df, c("Phylum"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
RelAbuAllSamples#Mean relative abundance for all samples

Trtdata<-read.csv("PhylumLevelRelAbu.csv", header=T) #Added in other category to reach 100%
Trtdata$Phylum = factor(Trtdata$Phylum, levels = c("Actinobacteria","Bacteroidetes","Firmicutes","Planctomycetes","Proteobacteria","Tenericutes","Verrucomicrobia","Other")) #fixes x-axis labels

PhylumAbu=ggplot(Trtdata, aes(x=FFG,y=mean))+geom_bar(aes(fill = Phylum),colour="black", stat="identity")+
  xlab("FFG")+ylab("Relative Abundance (%)")+theme(axis.title.x=element_blank())+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  scale_fill_manual(values=cbPalette)+facet_wrap(~Sampling_station,scale="free_x")+ theme(legend.key.size = unit(0.4, "cm"))+
  theme(legend.text = element_text(size = 8))+ theme(legend.background=element_blank())
PhylumAbu
dev.off()
tiff("Figures/PhylumLevelRelAbu.tiff", width = 3.5, height = 3.2, units = 'in', res = 300)
PhylumAbu
dev.off()


#######
#Figure 3C Phylum level facet plot
#######
#Phylum level comparisons between sites 

theme_set(theme_bw(base_size = 14)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
PhylumAll=tax_glom(GPr, "Phylum")
PhylumLevel = filter_taxa(PhylumAll, function(x) mean(x) > 1e-2, TRUE) #filter out any taxa lower tha 0.1%

df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "FFG"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

dfStatsPlot<-subset(Trtdata,Phylum!="Actinobacteria"&Phylum!= "Planctomycetes"&Phylum!="Tenericutes")
vec<-c("a","a","b","a","a","b","b","c","ab","a","a","b","ab","a","bc","c") #From plot above
dfStatsPlot # Only show significant taxa 
cdataplot=ggplot(dfStatsPlot, aes(x=FFG,y=mean))+geom_bar(aes(fill = FFG),colour="black", stat="identity")+ facet_grid(~Family)+xlab("FFG")+ylab("Relative Abundance (SEM)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.title.x=element_blank())+facet_wrap(~Phylum)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+geom_text(aes(x=FFG, y=mean+se+10,label=vec),size=3)+
  scale_fill_manual(values=cbPalette)+theme(legend.justification=c(0.05,0.95), legend.position=c(0.75,0.45))+ theme(legend.title = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=F)#+  theme(legend.text=element_text(size=rel(0.5)))#+  annotate("text", label = "KW, adj-P = 0.0014", x = 2, y = 15, size = 8)

label = c("         KW, adj-P = 0.0014", "  adj-P = 0.0014", "adj-P = 0.021","  adj-P = 0.0043")


cdataplot<-cdataplot + annotate("text",x=1.5, y=75,label=label,size=2.5)


PhylumFacet<-cdataplot

theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))


tiff("Figures/PhylumFacet.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
PhylumFacet
dev.off()
dev.off()

#######
#Figure 3D Family level facet plot
######

GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
FamilyAll=tax_glom(GPr, "Family")
FamilyLevel = filter_taxa(FamilyAll, function(x) mean(x) > 3e-2, TRUE) #filter out any taxa lower tha 0.1%

df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family", "FFG"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

dfStatsPlot<-subset(Trtdata,Family!="Enterobacteriaceae"&Family!= "Rhodobacteraceae"&Family!="Aeromonadaceae")
#dfStatsPlot

vec<-c("a","a","b","c","a","b","c","a","a","b","b","a","a","b","b","a")
#dfStatsPlot # Only show significant taxa 
cdataplot=ggplot(dfStatsPlot, aes(x=FFG,y=mean))+geom_bar(aes(fill = FFG),colour="black", stat="identity")+ facet_grid(~Family)+xlab("FFG")+ylab("Relative Abundance (SEM)") + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(axis.title.x=element_blank())+facet_wrap(~Family)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+geom_text(aes(x=FFG, y=mean+se+2,label=vec),size=3)+
  scale_fill_manual(values=cbPalette)+theme(legend.justification=c(0.05,0.95), legend.position=c(0.64,0.4))+ theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.key.size = unit(0.35, "cm"))+
  theme(legend.text = element_text(size = 8))+ theme(legend.background=element_blank())#+  annotate("text", label = "KW, adj-P = 0.0014", x = 2, y = 15, size = 8)

label = c("      KW\n adj-P < 0.001", "adj-P < 0.001", "adj-P < 0.001","           adj-P\n            < 0.001")


cdataplot<-cdataplot + annotate("text",x=3.5, y=18,label=label,size=2.5)

theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

FamilyFacet<-cdataplot
dev.off()
tiff("Figures/FamilyFacet.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
FamilyFacet
dev.off()





######
#Join together Figure 3
######

MassPlot
PhylumAbu
PhylumFacet
FamilyFacet


dev.off()
tiff("Figures/Figure3.tiff", width = 6.85, height = 6.85, units = 'in', res = 300)
ggarrange(MassPlot, PhylumAbu, PhylumFacet,FamilyFacet,
          labels = c("a", "b", "c","d"),
          ncol = 2, nrow = 2)
dev.off()





#####
#Figure 4
#####
GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
GenusAll=tax_glom(GPr,"Genus")
GenusAllStats<-GenusAll

ForestData=GenusAllStats#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
response <- as.factor(sample_data(ForestData)$FFG)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:19, ]#22
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important genera for classifying  samples\n by FFG")#\n in a string tells it to start a new line
#imp.20$MeanDecreaseGini
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
otunames
PredictorTable<-kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification
PredictorTable

GenusRandomForestSubset = subset_taxa(GenusAll, row.names(tax_table(GenusAll))%in% otunames)
GenusRandomForestSubset

df <- psmelt(GenusRandomForestSubset)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus", "FFG"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
cdataplot=ggplot(Trtdata, aes(x=FFG,y=mean))+geom_bar(aes(fill = FFG),colour="black", stat="identity")+ facet_grid(~Genus)+xlab("Feeding Group")+ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.x=element_blank())+facet_wrap(~Genus,scales = "free_y")+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+scale_fill_manual(values=cbPalette)
cdataplot

compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr",method="kruskal.test")

Means<-compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr")
Means
Means=compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr")
SigList<-length(unique(Trtdata$Genus))
for (i in levels(Means$Genus)){
  Tax<-i
  TaxAbundance<-subset(Means,Genus==i )
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  #print(Letters)
  SigList[i]<-Letters
  
}
vec<-unlist(SigList)
vec<-vec[-1]
cdataplot=ggplot(Trtdata, aes(x=FFG,y=mean))+geom_bar(aes(fill = FFG),colour="black", stat="identity")+ facet_wrap(~Genus)+xlab("FFG")+
  ylab("Relative Abundance (%)") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.x=element_blank())+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+geom_text(aes(x=FFG, y=mean+se+2,label=vec))+ scale_fill_manual(values=cbPalette)+
  theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))+
  theme(legend.justification=c(0.05,0.95), legend.position=c(0.64,0.2))+ theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme(legend.key.size = unit(0.35, "cm"))+
  theme(legend.text = element_text(size = 10))+ theme(legend.background=element_blank())+ theme(axis.title.x=element_blank())
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
cdataplot

dev.off()
tiff("Figures/Figure4.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
cdataplot
dev.off()



#######
#Figure 5A Faith's Div by FFG 
#######
FaithByFFG<-ggplot(metadata,aes(x=Sampling_station,y=faith_pd,fill=FFG))+
  geom_boxplot()+ scale_fill_manual(values=cbPalette)+xlab("")+ylab("Faith's Phylogenetic Diversity")+labs(fill="FFG")+
  guides(fill=FALSE)+facet_wrap(~FFG)+theme(axis.title.x=element_blank())#+geom_text(aes(x=FFG, y=shannon+se+1,label=Letters$Letters2), position=position_dodge(width=0.9), size=8,color="black")#+ 
FaithByFFG
theme_set(theme_bw(base_size = 11)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/FaithByFFG.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
FaithByFFG
dev.off()

#Figure 5B Shannon Div by FFG 
#######
ShannonByFFG<-ggplot(metadata,aes(x=Sampling_station,y=shannon,fill=FFG))+
  geom_boxplot()+ scale_fill_manual(values=cbPalette)+xlab("")+ylab("Shannon Diversity")+labs(fill="FFG")+
  guides(fill=FALSE)+facet_wrap(~FFG)+theme(axis.title.x=element_blank())#+geom_text(aes(x=FFG, y=shannon+se+1,label=Letters$Letters2), position=position_dodge(width=0.9), size=8,color="black")#+ 
ShannonByFFG
theme_set(theme_bw(base_size = 11)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

dev.off()
tiff("Figures/ShannonByFFG.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
ShannonByFFG
dev.off()

#####
#Figure 5C PCOA taxonomy jaccard
#####

ord=ordinate(physeq,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="FFG",shape ="FFG")+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplotTax<-ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 0, aes(fill = FFG))+theme(legend.justification=c(0.05,0.95), legend.position=c(0.71,0.36))+ theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+#+ theme(legend.key.size = unit(0.35, "cm"))
  theme(legend.text = element_text(size = 8))+ theme(legend.background=element_blank())+geom_point(size=2.5)
ordplotTax
dev.off()
tiff("Figures/PCoATax.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
ordplotTax
dev.off()

#####
#Figure 5D PCOA PICRUSt jaccard
#####
PiCRUST<-import_biom("C:\\Users\\Joe Receveur\\Downloads\\ItalyInvertPiCRUST2.20.19.biom")
PiCRUST<-merge_phyloseq(PiCRUST,sampdat)
PiCRUST=rarefy_even_depth(PiCRUST, 750000, replace = TRUE, trimOTUs = TRUE, verbose = TRUE,rngseed = TRUE)


GPdist=phyloseq::distance(PiCRUST, "jaccard")
adonis(GPdist ~ FFG*Sampling_station, as(sample_data(PiCRUST), "data.frame"))


ord=ordinate(PiCRUST,"PCoA", "jaccard")
ordplot=plot_ordination(physeq, ord,"samples", color="FFG",shape="FFG")+geom_point(size=2.5)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplotFunction<-ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 0, aes(fill = FFG))+theme(legend.justification=c(0.05,0.95), legend.position=c(0.71,0.36))+ theme(legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+#+ theme(legend.key.size = unit(0.35, "cm"))
  theme(legend.text = element_text(size = 8))+ theme(legend.background=element_blank())+geom_point(size=2.5)


ordplotFunction
dev.off()
tiff("Figures/PCoAFunction.tiff", width = 3.3, height = 3.3, units = 'in', res = 300)
ordplotFunction
dev.off()

#####
#Join together figure 5
#####


dev.off()
tiff("Figures/Figure5.tiff", width = 6.85, height = 6.85, units = 'in', res = 300)
ggarrange(FaithByFFG, ShannonByFFG, ordplotTax,ordplotFunction,
          labels = c("a", "b", "c","d"),
          ncol = 2, nrow = 2)
dev.off()






##########
#Figure Supplemental Random Forest indicator taxa Forest
############
rm(list=ls())
set.seed(27)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7","#E12D00")
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
biom=import_biom("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvert2018WTax.biom",parseFunction= parse_taxonomy_greengenes)
#tax_table(biom) <- tax_table(biom)[,-c(5:10,14)]#remove dummy ranks

metadata=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvertMetadataWDiversity3.5.2019.csv",header = TRUE)
head(metadata)

tree=read_tree("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvert11.29.18Tree.nwk")

sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$id
physeq=merge_phyloseq(biom,sampdat,tree)
physeq=rarefy_even_depth(physeq, 3000, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
GenusAll=tax_glom(GPr,"Genus")
GenusAllOstanaStats<-subset_samples(GenusAll, Sampling_station=="Forest")
GenusAllOstanaStats<-subset_samples(GenusAllOstanaStats, FFG!="Predator")

ForestData=GenusAllOstanaStats#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
response <- as.factor(sample_data(ForestData)$FFG)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:22, ]#22
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important genera for classifying  samples\n by FFG")#\n in a string tells it to start a new line
#imp.20$MeanDecreaseGini
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
otunames
PredictorTable<-kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification
PredictorTable


GenusRandomForestSubset = subset_taxa(GenusAllOstanaStats, row.names(tax_table(GenusAllOstanaStats))%in% otunames)
GenusRandomForestSubset

df <- psmelt(GenusRandomForestSubset)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus", "FFG"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
cdataplot=ggplot(Trtdata, aes(x=FFG,y=mean))+geom_bar(aes(fill = FFG),colour="black", stat="identity")+ facet_grid(~Genus)+xlab("Feeding Group")+ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.x=element_blank())+facet_wrap(~Genus,scales = "free_y")+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+scale_fill_manual(values=cbPalette)
cdataplot

compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr",method="kruskal.test")

Means<-compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr")
Means
Means=compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr")
SigList<-length(unique(Trtdata$Genus))
for (i in levels(Means$Genus)){
  Tax<-i
  TaxAbundance<-subset(Means,Genus==i )
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  #print(Letters)
  SigList[i]<-Letters
  
}
vec<-unlist(SigList)
vec<-vec[-1]
cdataplot=ggplot(Trtdata, aes(x=FFG,y=mean))+geom_bar(aes(fill = FFG),colour="black", stat="identity")+ facet_wrap(~Genus)+xlab("FFG")+ylab("Relative Abundance (%) Forest") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.title.x=element_blank())+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+geom_text(aes(x=FFG, y=mean+se+1,label=vec))+
  scale_fill_manual(values=cbPalette)+guides(fill=F)
cdataplot
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
tiff("Figures/RandomForestIndicatorsForest.tiff", width = 6.85, height = 6.85, units = 'in', res = 300)
cdataplot
dev.off()

##########
#Figure Supplemental Random Forest indicator taxa Alpine
############
rm(list=ls())
set.seed(27)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7","#E12D00")
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
biom=import_biom("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvert2018WTax.biom",parseFunction= parse_taxonomy_greengenes)
#tax_table(biom) <- tax_table(biom)[,-c(5:10,14)]#remove dummy ranks

metadata=read.csv("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvertMetadataWDiversity3.5.2019.csv",header = TRUE)
head(metadata)

tree=read_tree("C:\\Users\\Joe Receveur\\Documents\\MSU data\\ItalyInverts\\ItalyInvert11.29.18Tree.nwk")

sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$id
physeq=merge_phyloseq(biom,sampdat,tree)
physeq=rarefy_even_depth(physeq, 3000, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)



GPr  = transform_sample_counts(physeq, function(x) x / sum(x) ) #transform samples based on relative abundance
GPr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)
GenusAll=tax_glom(GPr,"Genus")

GenusAllOstanaStats<-subset_samples(GenusAll, Sampling_station!="Forest")


ForestData=GenusAllOstanaStats#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
response <- as.factor(sample_data(ForestData)$FFG)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:22, ]#22
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important genera for classifying  samples\n by FFG")#\n in a string tells it to start a new line
#imp.20$MeanDecreaseGini
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
otunames
PredictorTable<-kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification
PredictorTable


GenusRandomForestSubset = subset_taxa(GenusAllOstanaStats, row.names(tax_table(GenusAllOstanaStats))%in% otunames)
GenusRandomForestSubset

df <- psmelt(GenusRandomForestSubset)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Genus", "FFG"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
cdataplot=ggplot(Trtdata, aes(x=FFG,y=mean))+geom_bar(aes(fill = FFG),colour="black", stat="identity")+ facet_grid(~Genus)+xlab("Feeding Group")+ylab("Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.x=element_blank())+facet_wrap(~Genus,scales = "free_y")+geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+scale_fill_manual(values=cbPalette)
cdataplot

compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr",method="kruskal.test")

Means<-compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr")
Means
Means=compare_means(Abundance ~ FFG, data = df, group.by = "Genus", p.adjust.method = "fdr")
SigList<-length(unique(Trtdata$Genus))
for (i in levels(Means$Genus)){
  Tax<-i
  TaxAbundance<-subset(Means,Genus==i )
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  #print(Letters)
  SigList[i]<-Letters
  
}
vec<-unlist(SigList)
vec<-vec[-1]
cdataplot=ggplot(Trtdata, aes(x=FFG,y=mean))+geom_bar(aes(fill = FFG),colour="black", stat="identity")+ facet_wrap(~Genus)+xlab("FFG")+
  ylab("Relative Abundance (%) Alpine Prairie") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(axis.title.x=element_blank())+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se))+geom_text(aes(x=FFG, y=mean+se+1,label=vec))+ scale_fill_manual(values=cbPalette)+
  guides(fill=F)
#scale_x_discrete(labels=c("0 hrs", "24 hrs", "48 hrs","72 hrs"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))#+ scale_fill_manual(values=cbPalette)
cdataplot
tiff("Figures/RandomForestIndicatorsAlpine.tiff", width = 6.85, height = 6.85, units = 'in', res = 300)
cdataplot
dev.off()


