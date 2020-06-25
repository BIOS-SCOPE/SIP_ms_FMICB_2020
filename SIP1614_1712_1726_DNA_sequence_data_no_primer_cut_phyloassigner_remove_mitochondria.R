#install and load packages
install.packages("ggdendro")
# Load dependencies
library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")
library("dplyr")
library(gplots)
library(phyloseq)

library(ape)
library(ips)  #install package no compilation
library(vegan)
library(factoextra)
library(cluster)
source("https://bioconductor.org/biocLite.R")
biocLite("Heatplus") #update n
library(Heatplus)
library(pvclust)
library(fpc)
library(phangorn)

source("https://bioconductor.org/biocLite.R")
biocLite("DECIPHER")
library(DECIPHER)
library(grid)
library(gridExtra)
library(phytools) #no compilation
library(svglite)

library(clustsig)  #for simprof
library(scales)  #for pretty breaks function

source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
biocLite('DESeq2') 
biocLite('dendextend')
biocLite('tidyr')
biocLite('viridis')
biocLite('reshape')
install.packages("vegan") 
install.packages("ggplot2")

library("phyloseq")
library("vegan")
library("DESeq2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")

rm(list=ls())
getwd()
setwd("/Users/ShutingLiu/Desktop/Computer/F/UCSB_Research/Writing/SIP_paper/Fig/for_ms_6_3_20")
getwd()

#Create a phyloseq object out of your data files
#PRIOR TO THIS, in Excel temporarily replace "NA" with "Not_annotated"(use underscore, R doesn’t recognize space),remember to check “match case”, otherwise it will also replace ‘na’ in a taxa name. Otherwise R will remove rows with an unclassified order along w/ chloroplast data during the filter step
#Name file as "_NAchange". check sample names are all "_", no"-", will cause problems when merging for ggplot if sample names contain "-"
count_tab_temp <- read.table("SIP_table_otu_copy_withtaxa_NAchange.txt", header=T, row.names=NULL, check.names=F)
tax_tab_temp <- read.table("SIP_202_11_taxa_copy_NAchange.txt", header=T, row.names=NULL, check.names=F)

#Remove rows with chloroplast as the Class name (I found Eukaryota line and NA as kingdom in my data, remove that as well)
#Luis points out Alphaproteobacteria_Rickettsiales_Mitochondria is misannotation from silva V4, zooplankton bugs, should be excluded
OTU_nochl <- count_tab_temp %>% filter(Class !="Chloroplast" & Kingdom !="Eukaryota" & Kingdom !="Not_annotated" & Family != "Mitochondria")
TAX_nochl <- tax_tab_temp %>% filter(Class !="Chloroplast"& Kingdom !="Eukaryota" & Kingdom !="Not_annotated" & Family != "Mitochondria")


#write new csv files
write.csv(OTU_nochl, file="SIP_table_otu_copy_withtaxa_NAchange_noChl.csv", row.names=F)
write.csv(TAX_nochl, file="SIP_202_11_taxa_copy_NAchange_noChl.csv", row.names=F)

######After this just change the "Not_annotated" back to "NA" in excel, and convert files into .txt.
#convert csv file to tab delimited txt


# Force R to avoid storing strings as factors when uploading data
options(stringsAsFactors=FALSE)

#delete taxa name column in SIP_table_otu_copy_withtaxa_NAchange_noChl.txt file and save as"_otu"file
count_tab <- read.table("SIP_table_otu_copy_withtaxa_NAchange_noChl_otu.txt", header=T, row.names=1, check.names=F)
#reorder based on column names, put fraction from 1 to 10, then other ambient samples, as same order in the sample submission, which is the number after "S" in the sample name 
#extract column names
m <- colnames(count_tab)
#split the column names based on "_S"
r <- strsplit(m,"_S")
# Create a data.frame with two columns, X1 being names before _S, X2 being number after _S
q <- data.frame(matrix(unlist(r),nrow=length(m),byrow=T)) 
# Order q according to spec
class(q$X2)
#if no options stringAsFactors above, here X2 is factor, will need to convert to character first, then numeric
#q <- q[order(as.numeric(as.character(q$X2))),]
q <- q[order(as.numeric(q$X2)),]
# Reformat names
m <- paste0(q$X1,"_S",q$X2)
# Rearrange columns
count_tab_order <- count_tab[,m]

#create sample info file, put in sample names on first columnm, then experiment, treatment, density information
#save count_tab_order in a csv so can copy sample names
write.csv(count_tab_order, file="count_tab_order.csv")
#create csv then save as txt, read in txt
sample_info_tab <- read.table("Sample_info.txt", header=T, row.names=1,check.names=F, sep ="\t")


#class_order_family
tax_tab <- as.matrix(read.table("SIP_202_11_taxa_copy_NAchange_noChl.txt", header=T, row.names=1, check.names=F, na.strings="", sep="\t"))
#since a lot NA in family, include class and order name as well, we will create a new taxonomy level combining class order and family all three levels at once
# make new category with class order family levels, separated by _,insert after Family column, not at the end, because taxrank later concatenate from "upstream", don't put behind genus, otherwise will not concatenate correctly later
tax_tab<-data.frame(tax_tab) #convert to data frame first, otherwise next step won't work
tax_tab$multi_Taxonomy <- paste(tax_tab$Class,tax_tab$Order,tax_tab$Family, sep ="_")
tax_tab<-tax_tab[,c("Kingdom","Phylum","Class","Order","Family","multi_Taxonomy","Genus")]
#convert back to matrix for phyloseq object
tax_tab<-as.matrix(tax_tab)

#first look at how many reads get for each sample
colSums(count_tab_order)  #all samples and mock >10K except SIP1614 IJ5_F10 reads low (6354), all control clean <335 reads,except Blank_F5 contaminated by Moraxella, but this OTU not in other samples, so fine
#average of all samples except controls and mock
mean(colSums(count_tab_order[,grepl("SIP",colnames(count_tab_order))]))
min(colSums(count_tab_order[,grepl("SIP",colnames(count_tab_order))]))
max(colSums(count_tab_order[,grepl("SIP",colnames(count_tab_order))]))

#now ready to create phyloseq object
OTU = otu_table(count_tab_order, taxa_are_rows = TRUE)
TAX = tax_table(tax_tab)
SAM = sample_data(sample_info_tab)

physeqhet = phyloseq(OTU,TAX,SAM) #create phyloseq object with three components, you can also include fourth object of phylogenetic tree
physeqhet


# now we'll generate a proportions table for summarizing:
taxa_proportions_tab <- apply(count_tab_order, 2, function(x) x/sum(x)*100)
colSums(taxa_proportions_tab) #all should equal 100% except blank_whole all 0, NaN

# if we check the dimensions of this table at this point
dim(taxa_proportions_tab)

#save proportion file
write.csv(taxa_proportions_tab, file="SIP_table_otu_copy_NAchange_noChl_proportion.csv", row.names=T)

#Then add back taxa columns to the end of above proportion file, name as "_withtaxa"

#Proceed with all phyloseq analysis

####barplot using taxanomy
#We will first create a phyloseq object for only samples, leave mock, control behind, if just use whole phyloseq object, barplot using family will be too many rows due to mock and control community totally different and account for lots >5% rows
count_tab_order_Sample<-count_tab_order[,grepl("SIP",colnames(count_tab_order))]
sample_info_tab_Sample<-sample_info_tab[grepl("SIP",rownames(sample_info_tab)),]
OTU_Sample = otu_table(count_tab_order_Sample, taxa_are_rows = TRUE)
#TAX same as before
SAM_Sample = sample_data(sample_info_tab_Sample)

physeqhet_Sample = phyloseq(OTU_Sample,TAX,SAM_Sample) #create phyloseq object with three components, you can also include fourth object of phylogenetic tree
physeqhet_Sample

#barplot at the family level
family_counts_tab <- otu_table(tax_glom(physeqhet_Sample, taxrank="multi_Taxonomy")) 
dim(family_counts_tab)

family_tax_vec <- data.frame(tax_table(tax_glom(physeqhet_Sample, taxrank="multi_Taxonomy")))$multi_Taxonomy
family_tax_vec 
  
rownames(family_counts_tab) <- as.vector(family_tax_vec)
head(family_counts_tab)

#But NA_NA_NA as multiple rows here still in when using multi_Taxonomy, we need to manually filter it out
#filter doesn't work for phyloseq otu_tables, can use subset function if need to filter out some data from otu_tables.
annotated_family_counts_tab <- subset(family_counts_tab,rownames(family_counts_tab) != "NA_NA_NA")
unannotated_family_tax_counts <- colSums(count_tab_order_Sample) - colSums(annotated_family_counts_tab)
unannotated_family_tax_counts
family_and_unidentified_counts_tab <- rbind(annotated_family_counts_tab, "Unannotated"=unannotated_family_tax_counts)
View(family_and_unidentified_counts_tab)
rownames(family_and_unidentified_counts_tab)
  
family_taxa_proportions_tab <- apply(family_and_unidentified_counts_tab, 2, function(x) x/sum(x)*100)
colSums(family_taxa_proportions_tab)
dim(family_taxa_proportions_tab) #too many rows

# here, we'll only keep rows (taxa) that make up greater than 5% in any individual sample, try 0.1%, 0.2%,1%,2%,3%,4%,5%
temp_filt_family_taxa_proportions_tab <- data.frame(family_taxa_proportions_tab[apply(family_taxa_proportions_tab, 1, max,na.rm=TRUE) > 5, ]) 
dim(temp_filt_family_taxa_proportions_tab) # now we have 31 rows
View(temp_filt_family_taxa_proportions_tab)
  
filtered_family_proportions <- colSums(family_taxa_proportions_tab) - colSums(temp_filt_family_taxa_proportions_tab)
filt_family_taxa_proportions_tab <- rbind(temp_filt_family_taxa_proportions_tab, "Other(<=5%)"=filtered_family_proportions)  
View(filt_family_taxa_proportions_tab)

# first let's make a copy of our table that's safe for manipulating
filt_family_taxa_proportions_tab_for_plot <- filt_family_taxa_proportions_tab
# and add a column of the taxa names so that it is within the table, rather than just as row names
filt_family_taxa_proportions_tab_for_plot$Family_Taxa <- row.names(filt_family_taxa_proportions_tab_for_plot)

# now we'll transform the table into narrow, or long, format
filt_family_taxa_proportions_tab_for_plot.g <- gather(filt_family_taxa_proportions_tab_for_plot, Sample, Proportion, -Family_Taxa) 
#-Family_Taxa so this column keep at original position as a column, otherwise will be merged with all other sample info as a single column 

# take a look at the new table and compare it with the old one
head(filt_family_taxa_proportions_tab_for_plot.g)
head(filt_family_taxa_proportions_tab_for_plot)

# now we want a table with "color" and "characteristics" of each sample to merge into our plotting table so we can use that more easily in our plotting function
# here we're making a new table by pulling what we want from the sample information table
sample_info_for_merge<-data.frame("Sample"=row.names(sample_info_tab_Sample), "Experiment"=sample_info_tab_Sample$Experiment, "Treatment"=sample_info_tab_Sample$Treatment, "Density"=sample_info_tab_Sample$Density, stringsAsFactors=F)
sample_info_for_merge
# and here we are merging this table with the plotting table we just made (this is an awesome function!)
filt_family_taxa_proportions_tab_for_plot.g2 <- merge(filt_family_taxa_proportions_tab_for_plot.g, sample_info_for_merge,sort=F) #sort=F to keep original row order

Rhgcols_family<-c("darkolivegreen3","indianred4","dark blue","darkkhaki","red","darkseagreen","darkmagenta","deepskyblue","darkorange1","blue","forestgreen","brown","aquamarine3","salmon1","green","olivedrab1","plum","darkorchid1","rosybrown1","darkslategray4","cyan","lightskyblue1","honeydew2","tan3","goldenrod2","burlywood4","lightgoldenrod2","lavender","lightyellow","wheat1","wheat4","grey")

#convert sample and family name to factor, so it can keep at the original order in plot
filt_family_taxa_proportions_tab_for_plot.g2$Sample <- factor(filt_family_taxa_proportions_tab_for_plot.g2$Sample, levels=unique(filt_family_taxa_proportions_tab_for_plot.g2$Sample))
filt_family_taxa_proportions_tab_for_plot.g2$Family_Taxa <- factor(filt_family_taxa_proportions_tab_for_plot.g2$Family_Taxa, levels=unique(filt_family_taxa_proportions_tab_for_plot.g2$Family_Taxa))

View(filt_family_taxa_proportions_tab_for_plot.g2)
pdf("SIP1614_AB_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1614_A|SIP1614_B",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_family)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
    facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()

pdf("SIP1614_CD_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1614_C|SIP1614_D",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family)+
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
  facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()

pdf("SIP1614_EF_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1614_E|SIP1614_F",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_family)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
    facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()

pdf("SIP1614_GH_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1614_G|SIP1614_H",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family)+
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
  facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()

pdf("SIP1614_IJ_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1614_I|SIP1614_J",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_family)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
    facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()

pdf("SIP1614_KL_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1614_K|SIP1614_L",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family)+
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
  facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()
 
pdf("SIP1712_AB_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1712_A|SIP1712_B",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_family)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
    facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()
  
pdf("SIP1712_CD_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1712_C|SIP1712_D",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_family)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
    facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()
  
pdf("SIP1712_EF_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1712_E|SIP1712_F",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_family)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
    facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()
  
pdf("SIP1726_AB_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1726_A|SIP1726_B|SIP1726_Omnipore_A|SIP1726_Omnipore_B|SIP1726_GF75_A|SIP1726_GF75_B",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_family)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
    facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()
  
pdf("SIP1726_CD_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1726_C|SIP1726_D|SIP1726_Omnipore_C|SIP1726_Omnipore_D|SIP1726_GF75_C|SIP1726_GF75_D",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family)+
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
  facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()
  
pdf("SIP1726_EF_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1726_E|SIP1726_F|SIP1726_Omnipore_E|SIP1726_Omnipore_F|SIP1726_GF75_E|SIP1726_GF75_F",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family)+
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
  facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()

pdf("SIP1726_GH_Family_barplot.pdf",width=20,height=8)
ggplot(filt_family_taxa_proportions_tab_for_plot.g2[grepl("SIP1726_GH|SIP1726_Omnipore_G|SIP1726_Omnipore_H|SIP1726_GF75_G|SIP1726_GF75_H",filt_family_taxa_proportions_tab_for_plot.g2$Sample),], aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family)+
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
  facet_wrap(Experiment~Treatment,scales = "free") + theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()



####Use NMDS plots to check fractions in every sample and decide which fractions are outliers to exclude for SIP density comparison
#I will follow the online tutorial of phyloseq ordination to plot NMDS
#first create a new phyloseq object, use same otu table and tax table as for physeqhet_Sample above, but need to include sample name in the sample info for plotting label in NMDS plots
#copy sample_info_tab_Sample_order to a new data frame for putting sample name as a new column
sample_info_tab_Sample_NMDS <- sample_info_tab_Sample
sample_info_tab_Sample_NMDS$Sample_Name <- sapply(strsplit(row.names(sample_info_tab_Sample_NMDS),"_S"),`[`,1)
SAM_Sample_NMDS = sample_data(sample_info_tab_Sample_NMDS)
physeqhet_Sample_NMDS = phyloseq(OTU_Sample,TAX,SAM_Sample_NMDS)
physeqhet_Sample_NMDS
#Transform to even sampling depth:
physeqhet_Sample_NMDS_norm = transform_sample_counts(physeqhet_Sample_NMDS, function(x) 1E6 * x/sum(x))
#plot NMDS
#Look at SIP1614_AB
sample_names(physeqhet_Sample_NMDS_norm)
physeqhet_Sample_NMDS_norm_SIP1614_AB = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1614_A|SIP1614_B",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1614_AB)
physeqhet_Sample_NMDS_norm_SIP1614_AB.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1614_AB, "NMDS", "bray")
p1 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1614_AB, physeqhet_Sample_NMDS_norm_SIP1614_AB.ord)
pdf("SIP1614_AB_NMDS.pdf",width=10)
p1+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

#repeat for other samples
physeqhet_Sample_NMDS_norm_SIP1614_CD = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1614_C|SIP1614_D",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1614_CD)
physeqhet_Sample_NMDS_norm_SIP1614_CD.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1614_CD, "NMDS", "bray")
p2 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1614_CD, physeqhet_Sample_NMDS_norm_SIP1614_CD.ord)
pdf("SIP1614_CD_NMDS.pdf",width=10)
p2+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1614_EF = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1614_E|SIP1614_F",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1614_EF)
physeqhet_Sample_NMDS_norm_SIP1614_EF.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1614_EF, "NMDS", "bray")
p3 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1614_EF, physeqhet_Sample_NMDS_norm_SIP1614_EF.ord)
pdf("SIP1614_EF_NMDS.pdf",width=10)
p3+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1614_GH = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1614_G|SIP1614_H",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1614_GH)
physeqhet_Sample_NMDS_norm_SIP1614_GH.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1614_GH, "NMDS", "bray")
p4 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1614_GH, physeqhet_Sample_NMDS_norm_SIP1614_GH.ord)
pdf("SIP1614_GH_NMDS.pdf",width=10)
p4+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1614_IJ = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1614_I|SIP1614_J",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1614_IJ)
physeqhet_Sample_NMDS_norm_SIP1614_IJ.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1614_IJ, "NMDS", "bray")
p5 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1614_IJ, physeqhet_Sample_NMDS_norm_SIP1614_IJ.ord)
pdf("SIP1614_IJ_NMDS.pdf",width=10)
p5+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1614_KL = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1614_K|SIP1614_L",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1614_KL)
physeqhet_Sample_NMDS_norm_SIP1614_KL.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1614_KL, "NMDS", "bray")
p6 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1614_KL, physeqhet_Sample_NMDS_norm_SIP1614_KL.ord)
pdf("SIP1614_KL_NMDS.pdf",width=10)
p6+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1712_AB = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1712_A|SIP1712_B",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1712_AB)
physeqhet_Sample_NMDS_norm_SIP1712_AB.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1712_AB, "NMDS", "bray")
p7 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1712_AB, physeqhet_Sample_NMDS_norm_SIP1712_AB.ord)
pdf("SIP1712_AB_NMDS.pdf",width=10)
p7+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1712_CD = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1712_C|SIP1712_D",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1712_CD)
physeqhet_Sample_NMDS_norm_SIP1712_CD.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1712_CD, "NMDS", "bray")
p8 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1712_CD, physeqhet_Sample_NMDS_norm_SIP1712_CD.ord)
pdf("SIP1712_CD_NMDS.pdf",width=10)
p8+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1712_EF = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1712_E|SIP1712_F",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1712_EF)
physeqhet_Sample_NMDS_norm_SIP1712_EF.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1712_EF, "NMDS", "bray")#insufficient data
p9 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1712_EF, physeqhet_Sample_NMDS_norm_SIP1712_EF.ord)
pdf("SIP1712_EF_NMDS.pdf",width=10)
p9+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1726_AB = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1726_A|SIP1726_B|SIP1726_Omnipore_A|SIP1726_Omnipore_B|SIP1726_GF75_A|SIP1726_GF75_B",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1726_AB)
physeqhet_Sample_NMDS_norm_SIP1726_AB.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1726_AB, "NMDS", "bray")
p10 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1726_AB, physeqhet_Sample_NMDS_norm_SIP1726_AB.ord)
pdf("SIP1726_AB_NMDS.pdf",width=10)
p10+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1726_CD = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1726_C|SIP1726_D|SIP1726_Omnipore_C|SIP1726_Omnipore_D|SIP1726_GF75_C|SIP1726_GF75_D",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1726_CD)
physeqhet_Sample_NMDS_norm_SIP1726_CD.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1726_CD, "NMDS", "bray")
p11 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1726_CD, physeqhet_Sample_NMDS_norm_SIP1726_CD.ord)
pdf("SIP1726_CD_NMDS.pdf",width=10)
p11+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1726_EF = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1726_E|SIP1726_F|SIP1726_Omnipore_E|SIP1726_Omnipore_F|SIP1726_GF75_E|SIP1726_GF75_F",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1726_EF)
physeqhet_Sample_NMDS_norm_SIP1726_EF.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1726_EF, "NMDS", "bray")
p12 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1726_EF, physeqhet_Sample_NMDS_norm_SIP1726_EF.ord)
pdf("SIP1726_EF_NMDS.pdf",width=10)
p12+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()

physeqhet_Sample_NMDS_norm_SIP1726_GH = subset_samples(physeqhet_Sample_NMDS_norm,grepl("SIP1726_GH|SIP1726_Omnipore_G|SIP1726_Omnipore_H|SIP1726_GF75_G|SIP1726_GF75_H",sample_names(physeqhet_Sample_NMDS_norm)))
sample_names(physeqhet_Sample_NMDS_norm_SIP1726_GH)
physeqhet_Sample_NMDS_norm_SIP1726_GH.ord <- ordinate(physeqhet_Sample_NMDS_norm_SIP1726_GH, "NMDS", "bray")
p13 = plot_ordination(physeqhet_Sample_NMDS_norm_SIP1726_GH, physeqhet_Sample_NMDS_norm_SIP1726_GH.ord)
pdf("SIP1726_GH_NMDS.pdf",width=10)
p13+geom_point(size=3)+geom_text(mapping = aes(label = Sample_Name), size = 2.5,vjust=2) #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
dev.off()


#Based on NMDS, barplots and density range, I am planning to exlude F1 and F10.




####NMDS of non-fractioned samples

#######NMDS of non-fractioned samples#####
View(count_tab_order)
View(sample_info_tab)
#take out non-fractioned T0 and Harvesting time point samples, only take out 13C samples and no lignin
count_tab_order_nonSIP<-count_tab_order[,c(70,72,73,75,76,78,79,81,135,137,138,140,193:200)]
View(count_tab_order_nonSIP)
sample_info_tab_nonSIP<-sample_info_tab[c(70,72,73,75,76,78,79,81,135,137,138,140,193:200),]
View(sample_info_tab_nonSIP)
OTU_nonSIP = otu_table(count_tab_order_nonSIP, taxa_are_rows = TRUE)
sample_info_tab_nonSIP_NMDS <- sample_info_tab_nonSIP
sample_info_tab_nonSIP_NMDS$Sample_Name <- c(rep(c("S_Control","S_TW_lysate_PPL","M_Control","M_TW_lysate_PPL"),2),rep(c("M_Control","M_Amino_acid"),2),rep(c("M_Control","M_Syn_lysate","M_Syn_lysate_PPL","M_Syn_exudate_PPL"),2))
sample_info_tab_nonSIP_NMDS$Time_point <- c(rep("TF",4),rep("T0",4),rep("TF",2),rep("T0",2),rep("TF",4),rep("T0",4))
sample_info_tab_nonSIP_NMDS$Sample_Name<-factor(sample_info_tab_nonSIP_NMDS$Sample_Name,levels=unique(sample_info_tab_nonSIP_NMDS$Sample_Name)) #for color coding follow own order later
sample_info_tab_nonSIP_NMDS$Substrate<-c(rep(c("Control","TW_PPL"),4),rep(c("Control","Amino_acid"),2),rep(c("Control","Syn_lysate","Syn_PPL","Syn_PPL"),2))
sample_info_tab_nonSIP_NMDS$Substrate<-factor(sample_info_tab_nonSIP_NMDS$Substrate,levels=unique(sample_info_tab_nonSIP_NMDS$Substrate)) #for color coding follow own order later


View(sample_info_tab_nonSIP_NMDS)

SAM_nonSIP_NMDS = sample_data(sample_info_tab_nonSIP_NMDS)
physeqhet_nonSIP_NMDS = phyloseq(OTU_nonSIP,TAX,SAM_nonSIP_NMDS)
physeqhet_nonSIP_NMDS
#Transform to even sampling depth:
physeqhet_nonSIP_NMDS_norm = transform_sample_counts(physeqhet_nonSIP_NMDS, function(x) 1E6 * x/sum(x))
#plot NMDS
sample_names(physeqhet_nonSIP_NMDS_norm)
physeqhet_nonSIP_NMDS_norm.ord <- ordinate(physeqhet_nonSIP_NMDS_norm, "NMDS", "bray")
p_nonSIP = plot_ordination(physeqhet_nonSIP_NMDS_norm, physeqhet_nonSIP_NMDS_norm.ord)

pdf("NonSIP_NMDS_for_ms.pdf",width=10)
p_nonSIP+geom_point(size=5, aes(color=Sample_Name,shape=Time_point)) + #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(size = 14)) + theme(axis.title = element_text(size = 14), 
                                                       axis.text.y = element_text(size = 14)) +
  scale_x_continuous(limits=c(-2.7,1.7))
dev.off()
#stress 0.09496325,<0.1, fair fit

#try simprof
#extract phyloseq otu table as matrix
OTU1 = as(otu_table(physeqhet_nonSIP_NMDS_norm), "matrix")
View(OTU1)
# transpose if necessary
if(taxa_are_rows(physeqhet_nonSIP_NMDS_norm)){OTU1 <- t(OTU1)}
# Run simprof on the data
res <- simprof(data= OTU1, 
               method.distance="braycurtis")  #default alpha 0.05
# Graph the result
pl.color <- simprof.plot(res)
pdf("simprof_bray_cluster.pdf")
plot(pl.color)
dev.off()
#11 clusters


#to determine best cluster number, use same pam1 func as above
pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}

gskmn_bray = clusGap(physeqhet_nonSIP_NMDS_norm.ord$points, FUN=pam1, K.max = 10, B = 500,d.power=2)
gskmn_bray

pdf("Gap statistics_bray.pdf")
plot(gskmn_bray, main = "Gap statistic")
dev.off()
with(gskmn_bray, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
#two clusters

###remove Alteromonadaceae, then NMDS bray curtis, then cluster and simprof
physeqhet_nonSIP_NMDS_noAlt=subset_taxa(physeqhet_nonSIP_NMDS,multi_Taxonomy!="Gammaproteobacteria_Alteromonadales_Alteromonadaceae")
#Transform to even sampling depth:
physeqhet_nonSIP_NMDS_noAlt_norm = transform_sample_counts(physeqhet_nonSIP_NMDS_noAlt, function(x) 1E6 * x/sum(x))
#plot NMDS
physeqhet_nonSIP_NMDS_noAlt_norm.ord <- ordinate(physeqhet_nonSIP_NMDS_noAlt_norm, "NMDS", "bray")
p_nonSIP_noAlt = plot_ordination(physeqhet_nonSIP_NMDS_noAlt_norm, physeqhet_nonSIP_NMDS_noAlt_norm.ord)
sample_names(physeqhet_nonSIP_NMDS_noAlt_norm)
pdf("NonSIP_NMDS_noAlt_for_ms.pdf",width=10)
p_nonSIP_noAlt+geom_point(size=5, aes(color=Sample_Name,shape=Time_point)) + #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(size = 14)) + theme(axis.title = element_text(size = 14), 
                                                       axis.text.y = element_text(size = 14)) 
dev.off()
#stress 0.07632086,<0.1, fair fit


#color not very constrast, self define colors
rhg_cols_NMDS<-c("grey35","gold","grey","orange","plum","forestgreen","cyan","blue")
pdf("NonSIP_NMDS_noAlt_for_ms2.pdf",width=10)
p_nonSIP_noAlt+geom_point(size=5, aes(color=Sample_Name,shape=Time_point)) + #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
  scale_color_manual(values = rhg_cols_NMDS)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(size = 14)) + theme(axis.title = element_text(size = 14), 
                                                       axis.text.y = element_text(size = 14))  + theme(axis.ticks = element_line(colour = "black"), 
    panel.background = element_rect(fill = NA)) + theme(axis.text = element_text(colour = "black"), 
    panel.background = element_rect(colour = "black")) + theme(panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    panel.background = element_rect(size = 0.9))
dev.off()


#####cluster analysis to determine clustered groups
y_noAlt_bray= phyloseq::distance(physeqhet_nonSIP_NMDS_noAlt_norm, method="bray")
tiff(file = "Bray_Cluster_noAlt_Tree.tiff",units="in",width=5,height=5,res=150)
colorScale_noAlt_bray <- rainbow(length(levels(get_variable(physeqhet_nonSIP_NMDS_noAlt, "Sample_Name"))))
cols_noAlt_bray <- colorScale_noAlt_bray[get_variable(physeqhet_nonSIP_NMDS_noAlt, "Sample_Name")]
GP.tip.labels_noAlt_bray <- paste(as(get_variable(physeqhet_nonSIP_NMDS_noAlt, "Sample_Name"), "character"),as(get_variable(physeqhet_nonSIP_NMDS_noAlt, "Time_point"), "character"),sep="_")
GP.hclust_noAlt_bray <- hclust(y_noAlt_bray, method = "ward.D2")
plot(as.phylo(GP.hclust_noAlt_bray), show.tip.label = TRUE, tip.color = "white")
tiplabels(GP.tip.labels_noAlt_bray, col = cols_noAlt_bray, frame = "none", adj = -0.02, cex = 0.6, font = 0.5)
title("Bray-Curtis Distances")
dev.off()

#to determine best cluster number, use same pam1 func as above
pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}

gskmn_noAlt_bray = clusGap(physeqhet_nonSIP_NMDS_noAlt_norm.ord$points, FUN=pam1, K.max = 10, B = 500,d.power=2)
gskmn_noAlt_bray

pdf("Gap statistics_noAlt_bray.pdf")
plot(gskmn_noAlt_bray, main = "Gap statistic")
dev.off()
with(gskmn_noAlt_bray, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method="firstSEmax"))
#two clusters

OTU1_noAlt = as(otu_table(physeqhet_nonSIP_NMDS_noAlt_norm), "matrix")
# transpose if necessary
if(taxa_are_rows(physeqhet_nonSIP_NMDS_noAlt_norm)){OTU1_noAlt <- t(OTU1_noAlt)}
#Run simprof on the data
res_noAlt_bray <- simprof(data= OTU1_noAlt, 
                          method.distance="braycurtis")  #alpha 0.05
# Graph the result
pl.color_noAlt <- simprof.plot(res_noAlt_bray)
pdf("simprof_bray_cluster_noAlt.pdf")
plot(pl.color_noAlt)
dev.off()# 13 clusters, but big cluster makes more sense, 70% dissimilarity level 5 clusters
write.csv(OTU1_noAlt,"OTU1_noAlt.csv")





####if only plot NMDS for mesopelagic samples(remove Alteromonas)
OTU_nonSIP_M = otu_table(count_tab_order_nonSIP[,-c(1,2,5,6)], taxa_are_rows = TRUE)
View(sample_info_tab_nonSIP_NMDS)
SAM_nonSIP_NMDS_M = sample_data(sample_info_tab_nonSIP_NMDS[-c(1,2,5,6),])
physeqhet_nonSIP_NMDS_M = phyloseq(OTU_nonSIP_M,TAX,SAM_nonSIP_NMDS_M)
physeqhet_nonSIP_NMDS_M
physeqhet_nonSIP_NMDS_M_noAlt=subset_taxa(physeqhet_nonSIP_NMDS_M,multi_Taxonomy!="Gammaproteobacteria_Alteromonadales_Alteromonadaceae")
#Transform to even sampling depth:
physeqhet_nonSIP_NMDS_M_noAlt_norm = transform_sample_counts(physeqhet_nonSIP_NMDS_M_noAlt, function(x) 1E6 * x/sum(x))
#plot NMDS
physeqhet_nonSIP_NMDS_M_noAlt_norm.ord <- ordinate(physeqhet_nonSIP_NMDS_M_noAlt_norm, "NMDS", "bray")
p_nonSIP_M_noAlt = plot_ordination(physeqhet_nonSIP_NMDS_M_noAlt_norm, physeqhet_nonSIP_NMDS_M_noAlt_norm.ord)
sample_names(physeqhet_nonSIP_NMDS_M_noAlt_norm)
rhg_cols_NMDS_M<-c("grey","orange","plum","forestgreen","cyan","blue")
pdf("NonSIP_NMDS_M_noAlt_for_ms2.pdf",width=10)
p_nonSIP_M_noAlt+geom_point(size=5, aes(color=Sample_Name,shape=Time_point)) + #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
  scale_color_manual(values = rhg_cols_NMDS_M)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(size = 14)) + theme(axis.title = element_text(size = 14), 
                                                       axis.text.y = element_text(size = 14))  + theme(axis.ticks = element_line(colour = "black"), 
                                                                                                       panel.background = element_rect(fill = NA)) + theme(axis.text = element_text(colour = "black"), 
                                                                                                                                                           panel.background = element_rect(colour = "black")) + theme(panel.grid.major = element_line(linetype = "blank"), 



#barplot for non-fractioned samples to see which taxa driving the difference in NMDS plots
#barplot at the family level
family_counts_tab_nonSIP <- otu_table(tax_glom(physeqhet_nonSIP_NMDS, taxrank="multi_Taxonomy")) 
dim(family_counts_tab_nonSIP)

family_tax_vec_nonSIP <- data.frame(tax_table(tax_glom(physeqhet_nonSIP_NMDS, taxrank="multi_Taxonomy")))$multi_Taxonomy
family_tax_vec_nonSIP 

rownames(family_counts_tab_nonSIP) <- as.vector(family_tax_vec_nonSIP)
head(family_counts_tab_nonSIP)

#But NA_NA_NA as multiple rows here still in when using multi_Taxonomy, we need to manually filter it out
#filter doesn't work for phyloseq otu_tables, can use subset function if need to filter out some data from otu_tables.
annotated_family_counts_tab_nonSIP <- subset(family_counts_tab_nonSIP,rownames(family_counts_tab_nonSIP) != "NA_NA_NA")
unannotated_family_tax_counts_nonSIP <- colSums(family_counts_tab_nonSIP) - colSums(annotated_family_counts_tab_nonSIP)
unannotated_family_tax_counts_nonSIP
family_and_unidentified_counts_tab_nonSIP <- rbind(annotated_family_counts_tab_nonSIP, "Unannotated"=unannotated_family_tax_counts_nonSIP)
View(family_and_unidentified_counts_tab_nonSIP)
rownames(family_and_unidentified_counts_tab_nonSIP)

family_taxa_proportions_tab_nonSIP <- apply(family_and_unidentified_counts_tab_nonSIP, 2, function(x) x/sum(x)*100)
colSums(family_taxa_proportions_tab_nonSIP)
dim(family_taxa_proportions_tab_nonSIP) #too many rows

# here, we'll only keep rows (taxa) that make up greater than 5% in any individual sample, try 0.1%, 0.2%,1%,2%,5%
temp_filt_family_taxa_proportions_tab_nonSIP <- data.frame(family_taxa_proportions_tab_nonSIP[apply(family_taxa_proportions_tab_nonSIP, 1, max,na.rm=TRUE) > 5, ]) 
dim(temp_filt_family_taxa_proportions_tab_nonSIP) # now we have 19 rows
View(temp_filt_family_taxa_proportions_tab_nonSIP)

filtered_family_proportions_nonSIP <- colSums(family_taxa_proportions_tab_nonSIP) - colSums(temp_filt_family_taxa_proportions_tab_nonSIP)
filt_family_taxa_proportions_tab_nonSIP <- rbind(temp_filt_family_taxa_proportions_tab_nonSIP, "Other(<=5%)"=filtered_family_proportions_nonSIP)  
View(filt_family_taxa_proportions_tab_nonSIP)

# first let's make a copy of our table that's safe for manipulating
filt_family_taxa_proportions_tab_nonSIP_for_plot <- filt_family_taxa_proportions_tab_nonSIP
# and add a column of the taxa names so that it is within the table, rather than just as row names
filt_family_taxa_proportions_tab_nonSIP_for_plot$Family_Taxa <- row.names(filt_family_taxa_proportions_tab_nonSIP_for_plot)

# now we'll transform the table into narrow, or long, format
filt_family_taxa_proportions_tab_nonSIP_for_plot.g <- gather(filt_family_taxa_proportions_tab_nonSIP_for_plot, Sample, Proportion, -Family_Taxa) 
#-Family_Taxa so this column keep at original position as a column, otherwise will be merged with all other sample info as a single column 

# take a look at the new table and compare it with the old one
head(filt_family_taxa_proportions_tab_nonSIP_for_plot.g)
head(filt_family_taxa_proportions_tab_nonSIP_for_plot)

# now we want a table with "color" and "characteristics" of each sample to merge into our plotting table so we can use that more easily in our plotting function
# here we're making a new table by pulling what we want from the sample information table
sample_info_for_merge_nonSIP<-data.frame("Sample"=row.names(sample_info_tab_nonSIP), "Experiment"=sample_info_tab_nonSIP$Experiment, "Treatment"=sample_info_tab_nonSIP$Treatment, stringsAsFactors=F)
sample_info_for_merge_nonSIP
# and here we are merging this table with the plotting table we just made (this is an awesome function!)
filt_family_taxa_proportions_tab_nonSIP_for_plot.g2 <- merge(filt_family_taxa_proportions_tab_nonSIP_for_plot.g, sample_info_for_merge_nonSIP,sort=F) #sort=F to keep original row order

Rhgcols_family_nonSIP<-c("darkolivegreen3","indianred4","dark blue","darkkhaki","red","darkseagreen","darkmagenta","deepskyblue","darkorange1","blue","forestgreen","brown","aquamarine3","salmon1","green","olivedrab1","plum","darkorchid1","rosybrown1","grey")

#convert sample and family name to factor, so it can keep at the original order in plot
filt_family_taxa_proportions_tab_nonSIP_for_plot.g2$Sample <- factor(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2$Sample, levels=unique(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2$Sample))
filt_family_taxa_proportions_tab_nonSIP_for_plot.g2$Family_Taxa <- factor(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2$Family_Taxa, levels=unique(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2$Family_Taxa))

View(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2)
pdf("NonSIP_Family_barplot.pdf",width=30,height=15)
ggplot(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2, aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family_nonSIP)+
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance")+
  facet_wrap(~Experiment,scales = "free") +theme(axis.text.x = element_text(angle = 15))  + theme(axis.text.x = element_text(size = 10))
dev.off()

#barplot with order of samples and separate figures for manuscript
View(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2)
pdf("NonSIP_Family_barplot_SIP1614.pdf",width=20,height=10)
ggplot(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2[c(1:160),], aes(x=factor(Sample,levels = c("SIP1614_AB0_S76","SIP1614_A5B5_S70","SIP1614_EF0_S78","SIP1614_E5F5_S72","SIP1614_GH0_S79","SIP1614_G5H5_S73","SIP1614_KL0_S81","SIP1614_K5L5_S75")), y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family_nonSIP,name="Family")+
  scale_y_discrete(limits=c(0,25,50,75,100)) +
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies") + theme(panel.grid.major = element_line(linetype = "blank"), 
                                                                      panel.grid.minor = element_line(linetype = "blank"))+
  scale_x_discrete(labels=c("S_Control_T0", "S_Control_TF", "S_TW_lysate_SPE_T0","S_TW_lysate_SPE_TF","M_Control_T0","M_Control_TF","M_TW_lysate_SPE_T0","M_TW_lysate_SPE_TF")) +
  theme(axis.title = element_text(size = 16, 
                                  face = "bold"), axis.text = element_text(size = 16, 
                                                                           face = "bold"), axis.text.x = element_text(angle = 30,vjust=0.5,size = 16), 
        axis.text.y = element_text(size = 16), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) +labs(x = NULL)
dev.off()
pdf("NonSIP_Family_barplot_SIP1712.pdf",width=12,height=10)
ggplot(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2[c(161:240),], aes(x=factor(Sample,levels = c("SIP1712_AB0_S138","SIP1712_A6B6_S135","SIP1712_EF0_S140","SIP1712_E6F6_S137")), y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family_nonSIP,name="Family")+
  scale_y_discrete(limits=c(0,25,50,75,100)) +
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies") + theme(panel.grid.major = element_line(linetype = "blank"), 
                                                                      panel.grid.minor = element_line(linetype = "blank"))+
  scale_x_discrete(labels=c("M_Control_T0", "M_Control_TF", "M_Amino_acid_T0","M_Amino_acid_TF")) +
  theme(axis.title = element_text(size = 16, 
                                  face = "bold"), axis.text = element_text(size = 16, 
                                                                           face = "bold"), axis.text.x = element_text(angle = 30,vjust=0.5,size = 16), 
        axis.text.y = element_text(size = 16), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) +labs(x = NULL)
dev.off()
pdf("NonSIP_Family_barplot_SIP1726.pdf",width=20,height=10)
ggplot(filt_family_taxa_proportions_tab_nonSIP_for_plot.g2[c(241:400),], aes(x=factor(Sample,levels = c("SIP1726_Omnipore_AB0_S197","SIP1726_Omnipore_A5B5_S193","SIP1726_Omnipore_CD0_S198","SIP1726_Omnipore_C5D5_S194","SIP1726_Omnipore_EF0_S199","SIP1726_Omnipore_E5F5_S195","SIP1726_Omnipore_GH0_S200","SIP1726_Omnipore_G5H5_S196")), y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family_nonSIP,name="Family")+
  scale_y_discrete(limits=c(0,25,50,75,100)) +
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies") + theme(panel.grid.major = element_line(linetype = "blank"), 
                                                                      panel.grid.minor = element_line(linetype = "blank"))+
  scale_x_discrete(labels=c("M_Control_T0", "M_Control_TF", "M_Syn_lysate_T0","M_Syn_lysate_TF","M_Syn_lysate_SPE_T0","M_Syn_lysate_SPE_TF","M_Syn_exudate_SPE_T0","M_Syn_exudate_SPE_TF")) +
  theme(axis.title = element_text(size = 16, 
                                  face = "bold"), axis.text = element_text(size = 16, 
                                                                           face = "bold"), axis.text.x = element_text(angle = 30,vjust=0.5,size = 16), 
        axis.text.y = element_text(size = 16), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) +labs(x = NULL)
dev.off()

###Filter out copiotrophs, only look at oligotrophs for NMDS
#define copiotrophs (family):Alteromonadaceae, Vibrionaceae,Rhodobacteraceae,Flavobacteriaceae,Oceanospirillaceae,Erythrobacteraceae,Moraxellaceae,Pseudoalteromonadaceae,Oleiphilaceae,Shewanellaceae,Bacteroidaceae
physeqhet_nonSIP_oligotroph_NMDS=subset_taxa(physeqhet_nonSIP_NMDS,Family != "Alteromonadaceae" & Family != "Vibrionaceae"& Family != "Rhodobacteraceae"& Family != "Flavobacteriaceae"& Family != "Oceanospirillaceae"& Family != "Erythrobacteraceae"& Family != "Moraxellaceae"& Family != "Pseudoalteromonadaceae"& Family != "Oleiphilaceae"& Family != "Shewanellaceae"& Family != "Bacteroidaceae")
physeqhet_nonSIP_oligotroph_NMDS
physeqhet_nonSIP_oligotroph_NMDS_norm = transform_sample_counts(physeqhet_nonSIP_oligotroph_NMDS, function(x) 1E6 * x/sum(x))
sample_names(physeqhet_nonSIP_oligotroph_NMDS_norm)
physeqhet_nonSIP_oligotroph_NMDS_norm.ord <- ordinate(physeqhet_nonSIP_oligotroph_NMDS_norm, "NMDS", "bray")
p_nonSIP_oligotroph = plot_ordination(physeqhet_nonSIP_oligotroph_NMDS_norm, physeqhet_nonSIP_oligotroph_NMDS_norm.ord)
pdf("NonSIP_NMDS_oligotroph.pdf",width=10)
p_nonSIP_oligotroph+geom_point(size=5, aes(color=Sample_Name,shape=Time_point)) + #if I put label="Sample_Name" above in the plot_ordination, the size is too small, so separate layer here in ggplot
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text.x = element_text(size = 14)) + theme(axis.title = element_text(size = 14), 
                                                       axis.text.y = element_text(size = 14))
dev.off()


#Try ways to deal with high abundance Altermonas issue in Non-SIP barplot (suggestions get from Craig Nelson 4/7/2020)

#remove Alteromonas, redistribute other organisms
View(family_and_unidentified_counts_tab_nonSIP)
family_and_unidentified_counts_tab_nonSIP_no_Altero<-family_and_unidentified_counts_tab_nonSIP[-1,]
View(family_and_unidentified_counts_tab_nonSIP_no_Altero)
family_taxa_proportions_tab_nonSIP_no_Altero <- apply(family_and_unidentified_counts_tab_nonSIP_no_Altero, 2, function(x) x/sum(x)*100)
colSums(family_taxa_proportions_tab_nonSIP_no_Altero)
dim(family_taxa_proportions_tab_nonSIP_no_Altero) #too many rows

# here, we'll only keep rows (taxa) that make up greater than 5% in any individual sample, try 0.1%, 0.2%,1%,2%,3%
temp_filt_family_taxa_proportions_tab_nonSIP_no_Altero <- data.frame(family_taxa_proportions_tab_nonSIP_no_Altero[apply(family_taxa_proportions_tab_nonSIP_no_Altero, 1, max,na.rm=TRUE) > 3, ]) 
dim(temp_filt_family_taxa_proportions_tab_nonSIP_no_Altero) # now we have 30 rows
View(temp_filt_family_taxa_proportions_tab_nonSIP_no_Altero)

filtered_family_proportions_nonSIP_no_Altero <- colSums(family_taxa_proportions_tab_nonSIP_no_Altero) - colSums(temp_filt_family_taxa_proportions_tab_nonSIP_no_Altero)
filt_family_taxa_proportions_tab_nonSIP_no_Altero <- rbind(temp_filt_family_taxa_proportions_tab_nonSIP_no_Altero, "Other(<=3%)"=filtered_family_proportions_nonSIP_no_Altero)  
View(filt_family_taxa_proportions_tab_nonSIP_no_Altero)

# first let's make a copy of our table that's safe for manipulating
filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot <- filt_family_taxa_proportions_tab_nonSIP_no_Altero
# and add a column of the taxa names so that it is within the table, rather than just as row names
filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot$Family_Taxa <- row.names(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot)

# now we'll transform the table into narrow, or long, format
filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g <- gather(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot, Sample, Proportion, -Family_Taxa) 
#-Family_Taxa so this column keep at original position as a column, otherwise will be merged with all other sample info as a single column 

# take a look at the new table and compare it with the old one
head(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g)
head(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot)

#use same sample info as before
View(sample_info_for_merge_nonSIP)
# and here we are merging this table with the plotting table we just made (this is an awesome function!)
filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2 <- merge(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g, sample_info_for_merge_nonSIP,sort=F) #sort=F to keep original row order

Rhgcols_family_nonSIP_no_Altero<-c("darkolivegreen3","indianred4","dark blue","darkkhaki","red","darkseagreen","darkmagenta","deepskyblue","darkorange1","blue","forestgreen","brown","aquamarine3","salmon1","green","olivedrab1","plum","darkorchid1","rosybrown1","darkslategray4","cyan","lightskyblue1","honeydew2","tan3","goldenrod2","burlywood4","lightgoldenrod2","lavender","lightyellow","wheat1","grey")

#convert sample and family name to factor, so it can keep at the original order in plot
filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2$Sample <- factor(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2$Sample, levels=unique(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2$Sample))
filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2$Family_Taxa <- factor(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2$Family_Taxa, levels=unique(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2$Family_Taxa))

View(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2)
pdf("NonSIP_no_Altero_Family_barplot_SIP1614.pdf",width=20,height=10)
ggplot(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2[c(1:248),], aes(x=factor(Sample,levels = c("SIP1614_AB0_S76","SIP1614_A5B5_S70","SIP1614_EF0_S78","SIP1614_E5F5_S72","SIP1614_GH0_S79","SIP1614_G5H5_S73","SIP1614_KL0_S81","SIP1614_K5L5_S75")), y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family_nonSIP_no_Altero,name="Family")+
  scale_y_discrete(limits=c(0,25,50,75,100)) +
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies") + theme(panel.grid.major = element_line(linetype = "blank"), 
                                                                      panel.grid.minor = element_line(linetype = "blank"))+
  scale_x_discrete(labels=c("S Control T0", "S Control TF", "S TW lysate PPL T0","S TW lysate PPL TF","M Control T0","M Control TF","M TW lysate PPL T0","M TW lysate PPL TF")) +
  theme(axis.title = element_text(size = 16, 
                                  face = "bold"), axis.text = element_text(size = 16, 
                                                                           face = "bold"), axis.text.x = element_text(angle = 30,vjust=0.5,size = 16), 
        axis.text.y = element_text(size = 16), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) +labs(x = NULL)
dev.off()
pdf("NonSIP_no_Altero_Family_barplot_SIP1712.pdf",width=20,height=10)
ggplot(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2[c(249:372),], aes(x=factor(Sample,levels = c("SIP1712_AB0_S138","SIP1712_A6B6_S135","SIP1712_EF0_S140","SIP1712_E6F6_S137")), y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family_nonSIP_no_Altero,name="Family")+
  scale_y_discrete(limits=c(0,25,50,75,100)) +
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies") + theme(panel.grid.major = element_line(linetype = "blank"), 
                                                                      panel.grid.minor = element_line(linetype = "blank"))+
  scale_x_discrete(labels=c("M Control T0", "M Control TF", "M Amino acid T0","M Amino acid TF")) +
  theme(axis.title = element_text(size = 16, 
                                  face = "bold"), axis.text = element_text(size = 16, 
                                                                           face = "bold"), axis.text.x = element_text(angle = 30,vjust=0.5,size = 16), 
        axis.text.y = element_text(size = 16), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) +labs(x = NULL)
dev.off()
pdf("NonSIP_no_Altero_Family_barplot_SIP1726.pdf",width=20,height=10)
ggplot(filt_family_taxa_proportions_tab_nonSIP_no_Altero_for_plot.g2[c(373:620),], aes(x=factor(Sample,levels = c("SIP1726_Omnipore_AB0_S197","SIP1726_Omnipore_A5B5_S193","SIP1726_Omnipore_CD0_S198","SIP1726_Omnipore_C5D5_S194","SIP1726_Omnipore_EF0_S199","SIP1726_Omnipore_E5F5_S195","SIP1726_Omnipore_GH0_S200","SIP1726_Omnipore_G5H5_S196")), y=Proportion, fill=Family_Taxa))  +
  geom_bar(width=0.6, stat="identity") +
  scale_fill_manual(values =Rhgcols_family_nonSIP_no_Altero,name="Family")+
  scale_y_discrete(limits=c(0,25,50,75,100)) +
  theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies") + theme(panel.grid.major = element_line(linetype = "blank"), 
                                                                      panel.grid.minor = element_line(linetype = "blank"))+
  scale_x_discrete(labels=c("M Control T0", "M Control TF", "M Syn lysate T0","M Syn lysate TF","M Syn lysate PPL T0","M Syn lysate PPL TF","M Syn exudate PPL T0","M Syn exudate PPL TF")) +
  theme(axis.title = element_text(size = 16, 
                                  face = "bold"), axis.text = element_text(size = 16, 
                                                                           face = "bold"), axis.text.x = element_text(angle = 30,vjust=0.5,size = 16), 
        axis.text.y = element_text(size = 16), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12)) +labs(x = NULL)
dev.off()
#seems showing more difference among treatments now

#z scoring and view in heatmap
View(family_taxa_proportions_tab_nonSIP)
write.csv(family_taxa_proportions_tab_nonSIP,"family_taxa_proportions_tab_nonSIP.csv")
#255 family too many to plot on heatmap
#try above 0.5%
temp_filt_family_taxa_proportions_tab_nonSIP_z_score <- data.frame(family_taxa_proportions_tab_nonSIP[apply(family_taxa_proportions_tab_nonSIP, 1, max,na.rm=TRUE) > 0.5, ]) 
dim(temp_filt_family_taxa_proportions_tab_nonSIP_z_score) # now we have 54 rows, include all 13C incorporators except Piscirickettsiaceae
View(temp_filt_family_taxa_proportions_tab_nonSIP_z_score)
rownames(temp_filt_family_taxa_proportions_tab_nonSIP_z_score)
#remove unannotated
family_taxa_proportions_tab_nonSIP_z_score<-temp_filt_family_taxa_proportions_tab_nonSIP_z_score[-54,]
View(family_taxa_proportions_tab_nonSIP_z_score)

#can use scale function to calculate z score for each organism (minus mean and then divided by standard deviation)
#scale only works for matrix
temp_z_score<-matrix(ncol=20,nrow=53)
rownames(temp_z_score)<-rownames(family_taxa_proportions_tab_nonSIP_z_score)
colnames(temp_z_score)<-colnames(family_taxa_proportions_tab_nonSIP_z_score)
family_taxa_proportions_tab_nonSIP_z_score_matrix<-data.matrix(family_taxa_proportions_tab_nonSIP_z_score)
View(family_taxa_proportions_tab_nonSIP_z_score_matrix)
for (i in 1:53){
  temp_z_score[i,] <- scale(family_taxa_proportions_tab_nonSIP_z_score_matrix[i,],center=TRUE,scale=TRUE)  #not organized column order, center TRUE: minus mean, scale TRUE:divided by std
}
View(temp_z_score)
#convert matrix to data frame and then reorganize column into order
z_score_no_order<-data.frame(temp_z_score)
View(z_score_no_order)
#correct some oceanospirillales_NA to oligotrophs
z_score<-z_score_no_order[c(1,2,3,6,8,11,19,20,21,23,30,5,7,9,10,12:15,17,22,24:27,29,31:33,35:53,4,16,18,28,34),c(5,1,6,2,7,3,8,4,11,9,12,10,17,13,18,14,19,15,20,16)]#group copiotrophs, oligotrophs, and archaea
View(z_score)
colnames(z_score)<-c("July 2016 S Control T0","July 2016 S Control TF","July 2016 S TW lysate PPL T0","July 2016 S TW lysate PPL TF","July 2016 M Control T0","July 2016 M Control TF","July 2016 M TW lysate PPL T0","July 2016 M TW lysate PPL TF","July 2017 M Control T0","July 2017 M Control TF","July 2017 M Amino acid T0","July 2017 M Amino acid TF","Nov 2017 M Control T0","Nov 2017 M Control TF","Nov 2017 M Syn lysate T0","Nov 2017 M Syn lysate TF","Nov 2017 M Syn lysate PPL T0","Nov 2017 M Syn lysate PPL TF","Nov 2017 M Syn exudate PPL T0","Nov 2017 M Syn exudate PPL TF")
#if same name "M Control" for AE1712 and AE1726, won't plot the column, need unique column names
View(z_score)

#prepare for heatmap
z_score$Family_Taxa<-row.names(z_score)
##replace NA with empty in taxa name for plot
z_score$Family_Taxa<-gsub("_NA",'',z_score$Family_Taxa)
View(z_score)
z_score.long <- melt(z_score, id = "Family_Taxa")
View(z_score.long)

z_score.long$Family_Taxa <- factor(x = z_score.long$Family_Taxa,
                                     levels = z_score$Family_Taxa, 
                                     ordered = TRUE)  #keep original order

#heatmap
heatmap_z_score.plot <- ggplot(data = z_score.long, aes(x = variable, y = Family_Taxa)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="blue",mid="grey80",high="red")+  #default low is muted red, mid is white, high is muted blue
  theme(legend.position = "right") + theme(plot.subtitle = element_text(vjust = 1), 
                                           plot.caption = element_text(vjust = 1), 
                                           axis.text.x = element_text(angle=90,size=15,face="bold",color="black")) + 
  theme(legend.text = element_text(size = 15))  + theme(legend.title = element_text(size = 15))  + theme(legend.key.size = unit(3,"line")) + 
  theme(axis.title = element_text(size = 15)) +labs(x = NULL, y = NULL,fill="z score") + theme(axis.text.y = element_text(size = 20,face="bold"))
  
pdf("z_score_heatmap.pdf",width=20,height=20)
print(heatmap_z_score.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
dev.off()
#this may be the best way to present the data

##group S samples together, M samples together, group T0 together, TF together
heatmap_z_score_group.plot <- ggplot(data = z_score.long, aes(x = factor(variable,levels=c("July 2016 S Control T0","July 2016 S TW lysate PPL T0","July 2016 S Control TF","July 2016 S TW lysate PPL TF","July 2016 M Control T0","July 2016 M TW lysate PPL T0","July 2016 M Control TF","July 2016 M TW lysate PPL TF","July 2017 M Control T0","July 2017 M Amino acid T0","July 2017 M Control TF","July 2017 M Amino acid TF","Nov 2017 M Control T0","Nov 2017 M Syn lysate T0","Nov 2017 M Syn lysate PPL T0","Nov 2017 M Syn exudate PPL T0","Nov 2017 M Control TF","Nov 2017 M Syn lysate TF","Nov 2017 M Syn lysate PPL TF","Nov 2017 M Syn exudate PPL TF")), y = Family_Taxa)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="blue",mid="grey80",high="red")+  #default low is muted red, mid is white, high is muted blue
  theme(legend.position = "right") + theme(plot.subtitle = element_text(vjust = 1), 
                                           plot.caption = element_text(vjust = 1), 
                                           axis.text.x = element_text(angle=90,size=15,face="bold",color="black")) + 
  theme(legend.text = element_text(size = 15))  + theme(legend.title = element_text(size = 15))  + theme(legend.key.size = unit(3,"line")) + 
  theme(axis.title = element_text(size = 15)) +labs(x = NULL, y = NULL,fill="z-score") + theme(axis.text.y = element_text(size = 20,face="bold"))

pdf("z_score_heatmap_group.pdf",width=20,height=20)
print(heatmap_z_score_group.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
dev.off()

###if just compare z score between Control TF and amendment TF, calculate each taxa z score diff between these two
dim(z_score)
z_score_treatment_diff<-as.data.frame(matrix(nrow=53,ncol=6))
rownames(z_score_treatment_diff)<-z_score$Family_Taxa
colnames(z_score_treatment_diff)<-c("S TW lysate PPL","M TW lysate PPL","M Amino acid","M Syn lysate","M Syn lysate PPL","M Syn exudate PPL")
z_score_treatment_diff[,1]<-z_score[,4]-z_score[,2]
z_score_treatment_diff[,2]<-z_score[,8]-z_score[,6]
z_score_treatment_diff[,3]<-z_score[,12]-z_score[,10]
z_score_treatment_diff[,4]<-z_score[,16]-z_score[,14]
z_score_treatment_diff[,5]<-z_score[,18]-z_score[,14]
z_score_treatment_diff[,6]<-z_score[,20]-z_score[,14]
View(z_score_treatment_diff)
#prepare for heatmap
z_score_treatment_diff$Family_Taxa<-row.names(z_score_treatment_diff)
z_score_treatment_diff.long <- melt(z_score_treatment_diff, id = "Family_Taxa")
View(z_score_treatment_diff.long)

z_score_treatment_diff.long$Family_Taxa <- factor(x = z_score_treatment_diff.long$Family_Taxa,
                                   levels = z_score_treatment_diff$Family_Taxa, 
                                   ordered = TRUE)  #keep original order

#heatmap
heatmap_z_score_treatment_diff.plot <- ggplot(data = z_score_treatment_diff.long, aes(x = variable, y = Family_Taxa)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low="blue",mid="grey80",high="red")+  #default low is muted red, mid is white, high is muted blue
  theme(legend.position = "right") + theme(plot.subtitle = element_text(vjust = 1), 
                                           plot.caption = element_text(vjust = 1), 
                                           axis.text.x = element_text(angle=90,size=15,face="bold",color="black")) + 
  theme(legend.text = element_text(size = 15))  + theme(legend.title = element_text(size = 15))  + theme(legend.key.size = unit(3,"line")) + 
  theme(axis.title = element_text(size = 15)) +labs(x = NULL, y = NULL,fill="vs. control z score change") + theme(axis.text.y = element_text(size = 20,face="bold"))

pdf("z_score_treatment_diff_heatmap.pdf",width=20,height=20)
print(heatmap_z_score_treatment_diff.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
dev.off()

#subset z score vs.control TF change >1 or <-1
z_score_treatment_diff_above_1stdev<-z_score_treatment_diff[apply(z_score_treatment_diff[,-7],1, function(x) any(abs(x)>1)),]
View(z_score_treatment_diff_above_1stdev)
write.csv(z_score_treatment_diff_above_1stdev,"z_score_treatment_diff_above_1stdev.csv")



#also plot a scatterplot showing range of relative abundance for each organism, put next to the heatmap
family_taxa_proportions_tab_nonSIP_z_score$Family_Taxa<-row.names(family_taxa_proportions_tab_nonSIP_z_score)
family_taxa_proportions_tab_nonSIP_z_score$Family_Taxa<-gsub("_NA",'',family_taxa_proportions_tab_nonSIP_z_score$Family_Taxa)
family_taxa_proportions_tab_nonSIP_z_score.long <- melt(family_taxa_proportions_tab_nonSIP_z_score, id = "Family_Taxa")
View(family_taxa_proportions_tab_nonSIP_z_score.long)
family_taxa_proportions_tab_nonSIP_z_score.long$Family_Taxa <- factor(x = family_taxa_proportions_tab_nonSIP_z_score.long$Family_Taxa,
                                   levels = z_score$Family_Taxa, 
                                   ordered = TRUE)  #keep original order

#convert temp_z_score matrix to dataframe for melt
View(temp_z_score)
temp_z_score_dataframe<-data.frame(temp_z_score)
temp_z_score_dataframe$Family_Taxa<-row.names(temp_z_score_dataframe)
temp_z_score_dataframe.long <- melt(temp_z_score_dataframe, id = "Family_Taxa")
View(temp_z_score_dataframe.long)
family_taxa_proportions_tab_nonSIP_z_score.long$z_score<-temp_z_score_dataframe.long$value #temp_z_score family order same as taxa proportion
#add z score as color coding for raw scatterplot
View(family_taxa_proportions_tab_nonSIP_z_score.long)

pdf("raw_RA_scatterplot.pdf",width=5)
print(ggplot(data=family_taxa_proportions_tab_nonSIP_z_score.long,aes(x=value, y=Family_Taxa,color=z_score))+
        geom_point(size=2)+scale_color_gradient2(low="blue",mid="grey80",high="red")+
        theme(axis.text.x = element_text(size = 11)) +theme(axis.text.y = element_blank()) +
        labs(x = "Relative abundance (%)", y = NULL) + theme(plot.subtitle = element_text(vjust = 1), 
                                                             plot.caption = element_text(vjust = 1), 
                                                             axis.text.x = element_text(colour = "black"))  +#no row names for organisms
        theme(legend.position="none") +  #remove legend
        theme(panel.background = element_rect(fill = NA)) + theme(axis.line = element_line(linetype = "solid"))+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+  #add right top axis line
        scale_y_discrete(position="right")) #move tick to right y axis
        
       
dev.off()

pdf("raw_RA_scatterplot_log10scale.pdf",width=5)
print(ggplot(data=family_taxa_proportions_tab_nonSIP_z_score.long,aes(x=value, y=Family_Taxa,color=z_score))+
        geom_point(size=2)+scale_color_gradient2(low="blue",mid="grey80",high="red")+
        theme(axis.text.x = element_text(size = 11)) +theme(axis.text.y = element_blank()) +
        labs(x = expression("Percentage of relative abundance (log"[10]*" scale)"), y = NULL) + theme(plot.subtitle = element_text(vjust = 1), 
                                                             plot.caption = element_text(vjust = 1), 
                                                             axis.text.x = element_text(colour = "black"))  +#no row names for organisms
        theme(legend.position="none") +  #remove legend
        theme(panel.background = element_rect(fill = NA)) + theme(axis.line = element_line(linetype = "solid"))+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+  #add right top axis line
        scale_y_discrete(position="right") + #move tick to right y axis
        scale_x_log10(labels = scales::comma)+  #comma to use number instead of scientifc notation
        annotation_logticks(sides="b",short=unit(-0.1, "cm"), mid=unit(-0.2, "cm"), long=unit(-0.3,"cm")) +#negative to move tick outside of axis 
        # remove clipping
        coord_cartesian(clip="off")+ 
        # add space between ticks and labels
        theme(axis.text.x=element_text(margin = margin(t = 10))))#log scale ticks


dev.off()
write.csv(family_taxa_proportions_tab_nonSIP_z_score,"family_taxa_proportions_tab_nonSIP_z_score.csv")




####SIP fractionation data processing now
#If we try all family level with any SIP sample >0.1% to see coarse fraction OTU relative abundance vs density pattern, we can use above generated count tab as concatenated at family level and ggplot this time
family_taxa_proportions_tab_SIP <- family_taxa_proportions_tab[,grepl("SIP",colnames(family_taxa_proportions_tab))]
dim(family_taxa_proportions_tab_SIP)
temp_filt_family_taxa_proportions_tab_SIP <- data.frame(family_taxa_proportions_tab_SIP[apply(family_taxa_proportions_tab_SIP, 1, max) > 0.1, ])
dim(temp_filt_family_taxa_proportions_tab_SIP) #128 rows
View(temp_filt_family_taxa_proportions_tab_SIP)
  
filtered_family_proportions_SIP <- colSums(family_taxa_proportions_tab_SIP) - colSums(temp_filt_family_taxa_proportions_tab_SIP)
filt_family_taxa_proportions_tab_SIP <- rbind(temp_filt_family_taxa_proportions_tab_SIP, "Other(<=0.1%)"=filtered_family_proportions_SIP)  
View(filt_family_taxa_proportions_tab_SIP)
write.csv(filt_family_taxa_proportions_tab_SIP,"filt_family_taxa_proportions_tab_SIP.csv",row.names=T)
  
filt_family_taxa_proportions_tab_SIP_for_plot <- filt_family_taxa_proportions_tab_SIP
filt_family_taxa_proportions_tab_SIP_for_plot$Family_Taxa <- row.names(filt_family_taxa_proportions_tab_SIP_for_plot)
filt_family_taxa_proportions_tab_SIP_for_plot.g <- gather(filt_family_taxa_proportions_tab_SIP_for_plot, Sample, Proportion, -Family_Taxa) 
  
head(filt_family_taxa_proportions_tab_SIP_for_plot.g)
head(filt_family_taxa_proportions_tab_SIP_for_plot)
 
filt_family_taxa_proportions_tab_SIP_for_plot.g2 <- merge(filt_family_taxa_proportions_tab_SIP_for_plot.g, sample_info_for_merge[grepl("SIP",sample_info_for_merge$Sample),],sort=F) 
  
filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample <- factor(filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample, levels=unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample))
filt_family_taxa_proportions_tab_SIP_for_plot.g2$Family_Taxa <- factor(filt_family_taxa_proportions_tab_SIP_for_plot.g2$Family_Taxa, levels=unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2$Family_Taxa))
  
#we will only include F2-F9
View(filt_family_taxa_proportions_tab_SIP_for_plot.g2)
pdf("SIP1614_AB_Family_SIP.pdf") #we have 129 Family_Taxa, put 9 plots in one page, total 15 pages
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_AB <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1614_AB5_F2|SIP1614_AB5_F3|SIP1614_AB5_F4|SIP1614_AB5_F5|SIP1614_AB5_F6|SIP1614_AB5_F7|SIP1614_AB5_F8|SIP1614_AB5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_AB$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_AB[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_AB$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_AB$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
    geom_point() +
    geom_line() +
    facet_wrap(~Family_Taxa)+
    theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
    theme(strip.text.x = element_text(size = 5))+   #change facet title size to fit into one line
    theme(axis.text.x = element_text(angle = 15))  + 
    theme(axis.text.x = element_text(size = 4)))
  }
dev.off()
    
#repeat this for all others
    
pdf("SIP1614_CD_Family_SIP.pdf") #we have 130 Family_Taxa, put 9 plots in one page, total 15 pages
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_CD <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1614_CD5_F2|SIP1614_CD5_F3|SIP1614_CD5_F4|SIP1614_CD5_F5|SIP1614_CD5_F6|SIP1614_CD5_F7|SIP1614_CD5_F8|SIP1614_CD5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_AB$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_CD[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_CD$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_CD$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   #change facet title size to fit into one line
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

    
pdf("SIP1614_EF_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_EF <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1614_EF5_F2|SIP1614_EF5_F3|SIP1614_EF5_F4|SIP1614_EF5_F5|SIP1614_EF5_F6|SIP1614_EF5_F7|SIP1614_EF5_F8|SIP1614_EF5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_EF$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_EF[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_EF$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_EF$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   #change facet title size to fit into one line
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1614_GH_Family_SIP.pdf") ####note GH density fractions quite different from IJ and KL, use F1-F8 for GH to correct for this, F1 name use F1_,otherwise include F10 as well
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_GH <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1614_GH5_F1_|SIP1614_GH5_F2|SIP1614_GH5_F3|SIP1614_GH5_F4|SIP1614_GH5_F5|SIP1614_GH5_F6|SIP1614_GH5_F7|SIP1614_GH5_F8",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_GH$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_GH[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_GH$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_GH$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

  
pdf("SIP1614_IJ_Family_SIP.pdf")
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_IJ <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1614_IJ5_F2|SIP1614_IJ5_F3|SIP1614_IJ5_F4|SIP1614_IJ5_F5|SIP1614_IJ5_F6|SIP1614_IJ5_F7|SIP1614_IJ5_F8|SIP1614_IJ5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_IJ$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_IJ[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_IJ$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_IJ$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   #change facet title size to fit into one line
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1614_KL_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_KL <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1614_KL5_F2|SIP1614_KL5_F3|SIP1614_KL5_F4|SIP1614_KL5_F5|SIP1614_KL5_F6|SIP1614_KL5_F7|SIP1614_KL5_F8|SIP1614_KL5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_KL$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_KL[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_KL$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1614_KL$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

    
pdf("SIP1712_AB_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_AB <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1712_AB6_F2|SIP1712_AB6_F3|SIP1712_AB6_F4|SIP1712_AB6_F5|SIP1712_AB6_F6|SIP1712_AB6_F7|SIP1712_AB6_F8|SIP1712_AB6_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_AB$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_AB[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_AB$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_AB$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()
    
pdf("SIP1712_CD_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_CD <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1712_CD6_F2|SIP1712_CD6_F3|SIP1712_CD6_F4|SIP1712_CD6_F5|SIP1712_CD6_F6|SIP1712_CD6_F7|SIP1712_CD6_F8|SIP1712_CD6_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_CD$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_CD[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_CD$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_CD$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()
    
pdf("SIP1712_EF_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_EF <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1712_EF6_F2|SIP1712_EF6_F3|SIP1712_EF6_F4|SIP1712_EF6_F5|SIP1712_EF6_F6|SIP1712_EF6_F7|SIP1712_EF6_F8|SIP1712_EF6_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_EF$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_EF[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_EF$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1712_EF$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1726_AB_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_AB <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1726_AB5_F2|SIP1726_AB5_F3|SIP1726_AB5_F4|SIP1726_AB5_F5|SIP1726_AB5_F6|SIP1726_AB5_F7|SIP1726_AB5_F8|SIP1726_AB5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_AB$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_AB[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_AB$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_AB$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1726_CD_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_CD <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1726_CD5_F2|SIP1726_CD5_F3|SIP1726_CD5_F4|SIP1726_CD5_F5|SIP1726_CD5_F6|SIP1726_CD5_F7|SIP1726_CD5_F8|SIP1726_CD5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_CD$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_CD[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_CD$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_CD$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()
    
pdf("SIP1726_EF_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_EF <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1726_EF5_F2|SIP1726_EF5_F3|SIP1726_EF5_F4|SIP1726_EF5_F5|SIP1726_EF5_F6|SIP1726_EF5_F7|SIP1726_EF5_F8|SIP1726_EF5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_EF$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_EF[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_EF$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_EF$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()
    
pdf("SIP1726_GH_Family_SIP.pdf") 
filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_GH <- filt_family_taxa_proportions_tab_SIP_for_plot.g2[grepl("SIP1726_GH5_F2|SIP1726_GH5_F3|SIP1726_GH5_F4|SIP1726_GH5_F5|SIP1726_GH5_F6|SIP1726_GH5_F7|SIP1726_GH5_F8|SIP1726_GH5_F9",filt_family_taxa_proportions_tab_SIP_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_GH$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_GH[filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_GH$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_for_plot.g2_SIP1726_GH$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=Proportion))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="Proportion (%)")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()
    
 
####Since control also showed some GC smear and increase with density for certain OTUs, now I am going to do a log2(amend/control) on each OTU to compare this trend along density fraction from F2 to F9 (some from F1 to F8 as above).
filt_family_taxa_proportions_tab_SIP_logratio<- filt_family_taxa_proportions_tab_SIP
#create columns to put in the ratio numbers
#replace any 0 with 0.0001 (our minimal reads in SIP samples are 10604, threshold choose 2/10604=0.00019)
#Any possibility that a number is between 0 and 0.0001, if we assign 0 to 0.0001, then we make it larger than the number between 0 and 0.0001?
sum(filt_family_taxa_proportions_tab_SIP < 0.0001 & filt_family_taxa_proportions_tab_SIP != 0) 
#there is none, so we can just replace all 0 with 0.0001
filt_family_taxa_proportions_tab_SIP_logratio[filt_family_taxa_proportions_tab_SIP_logratio == 0] <- 0.0001
dim(filt_family_taxa_proportions_tab_SIP_logratio)
#then divide amended proportion of one OTU by the control proportion of that OTU and put in a new column, do this for F2 to F9 (except a few F1 to F8)
for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,218+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+11]/filt_family_taxa_proportions_tab_SIP_logratio[,i])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[218+i] <- paste("SIP1614_CD5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
#check ratio
filt_family_taxa_proportions_tab_SIP_logratio[,220]

#do this for all others as well
for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,226+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+22]/filt_family_taxa_proportions_tab_SIP_logratio[,i])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[226+i] <- paste("SIP1614_EF5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
filt_family_taxa_proportions_tab_SIP_logratio[,228]

#note GH fraction density quite different from IJ and KL, correct for this, use GH F1-F8
for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,234+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+44]/filt_family_taxa_proportions_tab_SIP_logratio[,i+33-1])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[234+i] <- paste("SIP1614_IJ5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
filt_family_taxa_proportions_tab_SIP_logratio[,236]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,242+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+55]/filt_family_taxa_proportions_tab_SIP_logratio[,i+33-1])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[242+i] <- paste("SIP1614_KL5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
filt_family_taxa_proportions_tab_SIP_logratio[,244]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,250+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+101]/filt_family_taxa_proportions_tab_SIP_logratio[,i+90])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[250+i] <- paste("SIP1712_CD6_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
filt_family_taxa_proportions_tab_SIP_logratio[,252]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,258+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+112]/filt_family_taxa_proportions_tab_SIP_logratio[,i+90])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[258+i] <- paste("SIP1712_EF6_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
filt_family_taxa_proportions_tab_SIP_logratio[,260]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,266+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+146]/filt_family_taxa_proportions_tab_SIP_logratio[,i+135])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[266+i] <- paste("SIP1726_CD5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
filt_family_taxa_proportions_tab_SIP_logratio[,268]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,274+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+156]/filt_family_taxa_proportions_tab_SIP_logratio[,i+135])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[274+i] <- paste("SIP1726_EF5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
filt_family_taxa_proportions_tab_SIP_logratio[,276]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_logratio[,282+i] = log2(filt_family_taxa_proportions_tab_SIP_logratio[,i+168]/filt_family_taxa_proportions_tab_SIP_logratio[,i+135])
  colnames(filt_family_taxa_proportions_tab_SIP_logratio)[282+i] <- paste("SIP1726_GH5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_logratio)
filt_family_taxa_proportions_tab_SIP_logratio[,284]

dim(filt_family_taxa_proportions_tab_SIP_logratio) 
filt_family_taxa_proportions_tab_SIP_logratio_for_plot <- filt_family_taxa_proportions_tab_SIP_logratio[,220:291]
colnames(filt_family_taxa_proportions_tab_SIP_logratio_for_plot)
filt_family_taxa_proportions_tab_SIP_logratio_for_plot$Family_Taxa <- row.names(filt_family_taxa_proportions_tab_SIP_logratio_for_plot)
View(filt_family_taxa_proportions_tab_SIP_logratio_for_plot)
write.csv(filt_family_taxa_proportions_tab_SIP_logratio_for_plot,"filt_family_taxa_proportions_tab_SIP_logratio_for_plot.csv",row.names=T)
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g <- gather(filt_family_taxa_proportions_tab_SIP_logratio_for_plot, Sample, logratio, -Family_Taxa) 
View(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g)

#create a new sample info with new sample names we created above
View(sample_info_tab_Sample)
sample_info_tab_SIP_logratio <- sample_info_tab_Sample[c(13:20,24:31,46:53,57:64,103:110,114:121,148:155,159:166,170:177),]
View(sample_info_tab_SIP_logratio)
rownames(sample_info_tab_SIP_logratio)<-sapply(strsplit(row.names(sample_info_tab_SIP_logratio),"_S"),`[`,1)
View(sample_info_tab_SIP_logratio)
sample_info_tab_SIP_logratio_for_merge<-data.frame("Sample"=row.names(sample_info_tab_SIP_logratio), "Experiment"=sample_info_tab_SIP_logratio$Experiment, "Treatment"=sample_info_tab_SIP_logratio$Treatment,"Density"=sample_info_tab_SIP_logratio$Density, stringsAsFactors=F)
  
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2 <- merge(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g, sample_info_tab_SIP_logratio_for_merge,sort=F) 
  
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample <- factor(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample, levels=unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample))
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Family_Taxa <- factor(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Family_Taxa, levels=unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Family_Taxa))
View(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2)
  
pdf("SIP1614_CD_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614CD <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1614_CD",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
View(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614CD)
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614CD$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614CD[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614CD$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614CD$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1614_EF_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614EF <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1614_EF",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
View(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614EF)
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614EF$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614EF[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614EF$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614EF$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1614_IJ_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614IJ <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1614_IJ",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614IJ$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614IJ[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614IJ$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614IJ$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1614_KL_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614KL <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1614_KL",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614KL$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614KL[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614KL$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1614KL$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1712_CD_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712CD <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1712_CD",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712CD$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712CD[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712CD$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712CD$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1712_EF_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712EF <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1712_EF",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712EF$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712EF[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712EF$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1712EF$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1726_CD_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726CD <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1726_CD",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726CD$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726CD[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726CD$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726CD$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1726_EF_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726EF <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1726_EF",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726EF$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726EF[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726EF$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726EF$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()

pdf("SIP1726_GH_Family_SIP_logratio.pdf") #only F2 to F9
filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726GH <- filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2[grepl("SIP1726_GH",filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2$Sample),]
for(i in seq(1,length(unique(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726GH$Family_Taxa)),9)){
  print(ggplot(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726GH[filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726GH$Family_Taxa %in% levels(filt_family_taxa_proportions_tab_SIP_logratio_for_plot.g2_1726GH$Family_Taxa)[i:(i+8)],], 
               aes(x=Density, y=logratio))  +
          geom_point() +
          geom_line() +
          facet_wrap(~Family_Taxa)+
          theme_bw()+ labs(x="Density (g/mL)", y="log2")+
          theme(strip.text.x = element_text(size = 5))+   
          theme(axis.text.x = element_text(angle = 15))  + 
          theme(axis.text.x = element_text(size = 4)))
}
dev.off()
  
#Now try mean density method for all 129 OTUs
#to calculate mean density for each OTU in every sample: 
#1)the relative abundance in a fraction was divided by the relative abundance in the same fraction of the unamended control (zero replace with 0.0001 or 2/min_reads): this gives a sense of "enrichment" of the bug to account for the fact that PCR amplifies stuff that isn't really there. Call this standardized abundance.
#2)The standardized abundance in each fraction is then divided by the total standardized abundance of all fractions...this number ranges from 0 to 1 and across all fractions totals 1. Conceptually this is a metric of the relative distribution across density.
#3)I then multiply this number times the density of each fraction...add this up and you get the median density of each bug. This approach gives you a quantitative density for each bug in each sample, which is better than a yes/no for uptake. The density is dictated by how well the population is labeled as a whole, so it works best on copiotrophs. 
  
#use above family SIP tab, total 129 OTUs
View(filt_family_taxa_proportions_tab_SIP)
filt_family_taxa_proportions_tab_SIP_mean_density<- filt_family_taxa_proportions_tab_SIP
  
filt_family_taxa_proportions_tab_SIP_mean_density[filt_family_taxa_proportions_tab_SIP_mean_density == 0] <- 0.0001
#then divide amended proportion of one OTU by the control proportion of that OTU and put in a new column, similar as above, do F2-F9 
for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,218+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+11]/filt_family_taxa_proportions_tab_SIP_mean_density[,i]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[218+i] <- paste("SIP1614_CD5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,220]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,226+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+22]/filt_family_taxa_proportions_tab_SIP_mean_density[,i]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[226+i] <- paste("SIP1614_EF5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,228]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,234+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+44]/filt_family_taxa_proportions_tab_SIP_mean_density[,i+33-1]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[234+i] <- paste("SIP1614_IJ5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,236]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,242+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+55]/filt_family_taxa_proportions_tab_SIP_mean_density[,i+33-1]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[242+i] <- paste("SIP1614_KL5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,244]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,250+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+101]/filt_family_taxa_proportions_tab_SIP_mean_density[,i+90]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[250+i] <- paste("SIP1712_CD6_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,252]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,258+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+112]/filt_family_taxa_proportions_tab_SIP_mean_density[,i+90]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[258+i] <- paste("SIP1712_EF6_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,260]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,266+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+146]/filt_family_taxa_proportions_tab_SIP_mean_density[,i+135]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[266+i] <- paste("SIP1726_CD5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,268]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,274+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+156]/filt_family_taxa_proportions_tab_SIP_mean_density[,i+135]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[274+i] <- paste("SIP1726_EF5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,276]

for (i in 2:9){
  filt_family_taxa_proportions_tab_SIP_mean_density[,282+i] = filt_family_taxa_proportions_tab_SIP_mean_density[,i+168]/filt_family_taxa_proportions_tab_SIP_mean_density[,i+135]
  colnames(filt_family_taxa_proportions_tab_SIP_mean_density)[282+i] <- paste("SIP1726_GH5_F",as.character(i),sep="") 
}
colnames(filt_family_taxa_proportions_tab_SIP_mean_density)
filt_family_taxa_proportions_tab_SIP_mean_density[,284]
  
#create a data frame to put in normalized proportion of above standardized relative abundance across fractions 2-9 (1-8) for each OTU
#subset standardized ralative abundance
sub_filt_family_taxa_proportions_tab_SIP_mean_density<-filt_family_taxa_proportions_tab_SIP_mean_density[,c(220:291)]
View(sub_filt_family_taxa_proportions_tab_SIP_mean_density)
dim(sub_filt_family_taxa_proportions_tab_SIP_mean_density)
normalized_standardized_RA <- data.frame(matrix(ncol=72,nrow=129))
#y for row, i for column, x to control move to next 8 fractions after every sample
for(y in 1:129){
    for(x in 0:8){
      for (i in (1+8*x):(8+8*x)){
      normalized_standardized_RA[y,i]=sub_filt_family_taxa_proportions_tab_SIP_mean_density[y,i]/sum(sub_filt_family_taxa_proportions_tab_SIP_mean_density[y,c((1+8*x):(8+8*x))])
      }
    }
  }
  
colnames(normalized_standardized_RA)<-colnames(sub_filt_family_taxa_proportions_tab_SIP_mean_density)
View(normalized_standardized_RA)
dim(normalized_standardized_RA)
sum(normalized_standardized_RA[1,1:8]) #check,should equal 1
sum(normalized_standardized_RA[2,9:16]) #check,should equal 1
  
#use density vector, density1 for 1614CD,density2 for 1614EF,density3 for 1614IJ,density4 for 1614KL,density5 for 1712CD,density6 for 1712EF,density7 for 1726CD,density8 for 1726EF,density9 for 1726GH
density1<-c(1.772,1.765,1.754,1.743,1.733,1.724,1.715,1.706,1.697,1.605)
density2<-c(1.776,1.766,1.756,1.746,1.737,1.726,1.717,1.707,1.697,1.619)
density3<-c(1.765,1.759,1.752,1.743,1.733,1.723,1.715,1.706,1.697,1.620)
density4<-c(1.766,1.761,1.752,1.743,1.734,1.724,1.716,1.706,1.697,1.644)
density5<-c(1.785,1.778,1.767,1.756,1.747,1.737,1.727,1.719,1.709,1.629)
density6<-c(1.784,1.777,1.767,1.756,1.747,1.737,1.727,1.718,1.709,1.646)
density7<-c(1.781,1.777,1.769,1.759,1.749,1.739,1.730,1.721,1.712,1.672)
density8<-c(1.778,1.773,1.766,1.757,1.748,1.738,1.728,1.720,1.712,1.639)
density9<-c(1.778,1.773,1.766,1.757,1.748,1.737,1.727,1.718,1.712,1.639)

###normalize density based on average whole density between 1614 and 1712_1726 batch to eliminate density difference due to operation and CsCl buffer difference
#whole density: 1.730,1.730,1.731,1.728,1.728,1.727,1.742,1.742,1.741,1.742,1.742,1.741,1.741
ave_whole_density<-mean(c(1.730,1.730,1.731,1.728,1.728,1.727,1.742,1.742,1.741,1.742,1.742,1.741,1.741))
density1_norm<-density1*ave_whole_density/1.730
density2_norm<-density2*ave_whole_density/1.731
density3_norm<-density3*ave_whole_density/1.728
density4_norm<-density4*ave_whole_density/1.727
density5_norm<-density5*ave_whole_density/1.742
density6_norm<-density6*ave_whole_density/1.741
density7_norm<-density7*ave_whole_density/1.742
density8_norm<-density8*ave_whole_density/1.741
density9_norm<-density9*ave_whole_density/1.741

#create a data frame to put in mean density
mean_density<-data.frame(matrix(ncol=9,nrow=129))
for(i in 1:129){
  mean_density[i,1]=normalized_standardized_RA[i,1]*density1_norm[2]+normalized_standardized_RA[i,2]*density1_norm[3]+normalized_standardized_RA[i,3]*density1_norm[4]+normalized_standardized_RA[i,4]*density1_norm[5]+normalized_standardized_RA[i,5]*density1_norm[6]+normalized_standardized_RA[i,6]*density1_norm[7]+normalized_standardized_RA[i,7]*density1_norm[8]+normalized_standardized_RA[i,8]*density1_norm[9]
  
  mean_density[i,2]=normalized_standardized_RA[i,9]*density2_norm[2]+normalized_standardized_RA[i,10]*density2_norm[3]+normalized_standardized_RA[i,11]*density2_norm[4]+normalized_standardized_RA[i,12]*density2_norm[5]+normalized_standardized_RA[i,13]*density2_norm[6]+normalized_standardized_RA[i,14]*density2_norm[7]+normalized_standardized_RA[i,15]*density2_norm[8]+normalized_standardized_RA[i,16]*density2_norm[9]
  
  mean_density[i,3]=normalized_standardized_RA[i,17]*density3_norm[2]+normalized_standardized_RA[i,18]*density3_norm[3]+normalized_standardized_RA[i,19]*density3_norm[4]+normalized_standardized_RA[i,20]*density3_norm[5]+normalized_standardized_RA[i,21]*density3_norm[6]+normalized_standardized_RA[i,22]*density3_norm[7]+normalized_standardized_RA[i,23]*density3_norm[8]+normalized_standardized_RA[i,24]*density3_norm[9]
  
  mean_density[i,4]=normalized_standardized_RA[i,25]*density4_norm[2]+normalized_standardized_RA[i,26]*density4_norm[3]+normalized_standardized_RA[i,27]*density4_norm[4]+normalized_standardized_RA[i,28]*density4_norm[5]+normalized_standardized_RA[i,29]*density4_norm[6]+normalized_standardized_RA[i,30]*density4_norm[7]+normalized_standardized_RA[i,31]*density4_norm[8]+normalized_standardized_RA[i,32]*density4_norm[9]
  
  mean_density[i,5]=normalized_standardized_RA[i,33]*density5_norm[2]+normalized_standardized_RA[i,34]*density5_norm[3]+normalized_standardized_RA[i,35]*density5_norm[4]+normalized_standardized_RA[i,36]*density5_norm[5]+normalized_standardized_RA[i,37]*density5_norm[6]+normalized_standardized_RA[i,38]*density5_norm[7]+normalized_standardized_RA[i,39]*density5_norm[8]+normalized_standardized_RA[i,40]*density5_norm[9]
  
  mean_density[i,6]=normalized_standardized_RA[i,41]*density6_norm[2]+normalized_standardized_RA[i,42]*density6_norm[3]+normalized_standardized_RA[i,43]*density6_norm[4]+normalized_standardized_RA[i,44]*density6_norm[5]+normalized_standardized_RA[i,45]*density6_norm[6]+normalized_standardized_RA[i,46]*density6_norm[7]+normalized_standardized_RA[i,47]*density6_norm[8]+normalized_standardized_RA[i,48]*density6_norm[9]
  
  mean_density[i,7]=normalized_standardized_RA[i,49]*density7_norm[2]+normalized_standardized_RA[i,50]*density7_norm[3]+normalized_standardized_RA[i,51]*density7_norm[4]+normalized_standardized_RA[i,52]*density7_norm[5]+normalized_standardized_RA[i,53]*density7_norm[6]+normalized_standardized_RA[i,54]*density7_norm[7]+normalized_standardized_RA[i,55]*density7_norm[8]+normalized_standardized_RA[i,56]*density7_norm[9]
  
  mean_density[i,8]=normalized_standardized_RA[i,57]*density8_norm[2]+normalized_standardized_RA[i,58]*density8_norm[3]+normalized_standardized_RA[i,59]*density8_norm[4]+normalized_standardized_RA[i,60]*density8_norm[5]+normalized_standardized_RA[i,61]*density8_norm[6]+normalized_standardized_RA[i,62]*density8_norm[7]+normalized_standardized_RA[i,63]*density8_norm[8]+normalized_standardized_RA[i,64]*density8_norm[9]
  
  mean_density[i,9]=normalized_standardized_RA[i,65]*density9_norm[2]+normalized_standardized_RA[i,66]*density9_norm[3]+normalized_standardized_RA[i,67]*density9_norm[4]+normalized_standardized_RA[i,68]*density9_norm[5]+normalized_standardized_RA[i,69]*density9_norm[6]+normalized_standardized_RA[i,70]*density9_norm[7]+normalized_standardized_RA[i,71]*density9_norm[8]+normalized_standardized_RA[i,72]*density9_norm[9]
  
  }

colnames(mean_density)<-c("12C_TW_lys_PPL_S","13C_TW_lys_PPL_S","12C_TW_lys_PPL_M","13C_TW_lys_PPL_M","Lignin_M","Amino_Acids_M","Syn_lys_M","Syn_lys_PPL_M","Syn_exu_PPL_M")
rownames(mean_density)<-rownames(sub_filt_family_taxa_proportions_tab_SIP_mean_density)
View(mean_density)


  
#To make our figure, we will build the two plots (the cluster diagram and the heatmap) separately, then use the grid framework to put them together. We start by making the dendrogram (or cluster).
# Run clustering
mean_density_cluster <- as.matrix(mean_density)
mean_density_cluster.dendro <- as.dendrogram(hclust(d = dist(x = mean_density_cluster)))
  
# Create dendro
dendro.plot <- ggdendrogram(data = mean_density_cluster.dendro, rotate = TRUE)
dendro.plot <- dendro.plot + theme(axis.text.y = element_text(size = 6))
# Preview the plot
print(dendro.plot)
  
#Now prepare for heatmap
#We need to start by transforming data to a “long” format, where each row only contains values for one measurement.
#We use the reshape2 package’s melt function to handle the data transformation:
#Frist make a variable to store taxa name
mean_density$Family_Taxa<-rownames(mean_density)
View(mean_density)
# Data wrangling
mean_density.long <- melt(mean_density, id = "Family_Taxa")
View(mean_density.long)

#To make the tips of the dendrogram in the same order as the y-axis of the heatmap, we’ll need to re-order the heatmap to match the structure of the clustering.
#It would be nice to have the dendrogram and heatmap closer together, without a legend separating them.
#The dendrogram should be vertically stretched so each tip lines up with a row in the heatmap plot.
  
#We’ll use the order of the tips in the dendrogram to re-order the rows in our heatmap. The order of those tips are stored in the mean_density_cluster.dendro object, and we can extract it with the order.dendro function:
mean_density.order <- order.dendrogram(mean_density_cluster.dendro)
  
#Using the order of the dendrogram tips (we stored this above), we can re-level the mean_density.long$Family_Taxa column to match the order of the tips in the trees.
# Order the levels according to their position in the cluster
mean_density.long$Family_Taxa <- factor(x = mean_density.long$Family_Taxa,
                                 levels = mean_density$Family_Taxa[mean_density.order], 
                                 ordered = TRUE)
  
# Heatmap
#In the ggplot package, we use the geom_tile layer for creating a heatmap. Given our prior experience with the y-axis labels being large, we will again use theme to make the accession numbers (the y-axis labels) a little smaller:
#We would like the legend to be somewere else besides in between the two graphs, so we will move it to the top of the heatmap. We do this with through the theeme layer when creating the heatmap, setting the value of legend.position to “top”:
#we can remove the Family_Taxa from the heatmap, as they are already displayed at the tips of the dendrogram (and we made sure the heatmap rows and dendrogram tips match, above). 
#We do so by setting the relevant theme elements to element_blank() in our call to ggplot:
heatmap.plot <- ggplot(data = mean_density.long, aes(x = variable, y = Family_Taxa)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colors = c("navy", "blue", "cyan", "lightcyan","yellow", "red", "red4")) +  #custom color for heat map
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top") + theme(plot.subtitle = element_text(vjust = 1), 
                                         plot.caption = element_text(vjust = 1), 
                                         axis.text.x = element_text(angle = 10))

#Finally, we will need to adjust the height and alignment of dendrogram to ensure the tips line up with the corresponding rows in the heatmap. Unfortunately, because the two figures are actually independent graphics objects, this part involves considerable trial and error. 
#In general, we will need to change two parameters in the viewport function when we call print on dendro.plot:
#y: This argument sets the vertical justification. A value of 0.5 centers the graphic in the middle of the viewport to which it is assigned. Values closer to zero move the graphic towards the bottom of the viewport and values closer to one move the graphic towards the top of the viewport.
#height: This is the proportion of the viewport’s vertical space the graphic should use. We initially set it to use all the vertical space (height = 1.0), but since the legend now occupies the real estate at the top of the graphic, we will need to make the dendrogram height smaller.
#Note: The actual values you use will depend on the graphics device you ultimately use
pdf("mean_density.pdf",width=20,height=20)
grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.485, width = 0.2, height = 1.032))#adjust y and height to match pdf size
dev.off()
  
#remove C12 and lignin
View(mean_density)
mean_density_13C <- mean_density[,-c(1,3,5,10)] #remember to remove the Family_Taxa column
mean_density_cluster_13C <- as.matrix(mean_density_13C)
mean_density_cluster_13C.dendro <- as.dendrogram(hclust(d = dist(x = mean_density_cluster_13C)))
dendro_13C.plot <- ggdendrogram(data = mean_density_cluster_13C.dendro, rotate = TRUE)
dendro_13C.plot <- dendro_13C.plot + theme(axis.text.y = element_text(size = 6))

mean_density_13C$Family_Taxa<-rownames(mean_density_13C)
mean_density_13C.long <- melt(mean_density_13C, id = "Family_Taxa")
View(mean_density_13C.long)
  
mean_density_13C.order <- order.dendrogram(mean_density_cluster_13C.dendro)
mean_density_13C.long$Family_Taxa <- factor(x = mean_density_13C.long$Family_Taxa,
                                          levels = mean_density_13C$Family_Taxa[mean_density_13C.order], 
                                          ordered = TRUE)
  
heatmap_13C.plot <- ggplot(data = mean_density_13C.long, aes(x = variable, y = Family_Taxa)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colors = c("navy", "blue", "cyan", "lightcyan","yellow", "red", "red4")) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top") + theme(plot.subtitle = element_text(vjust = 1), 
                                           plot.caption = element_text(vjust = 1), 
                                           axis.text.x = element_text(angle = 10))
  
pdf("mean_density_13C_nolig.pdf",width=20,height=20)
grid.newpage()
print(heatmap_13C.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro_13C.plot, vp = viewport(x = 0.90, y = 0.485, width = 0.2, height = 1.032))
dev.off()

  
###I am going to choose some major taxa for ploting the mean density. Two criteria:
#1) Taxa is widespread across amended treatments meaning they are non-zero in at least three fractions (out of F2-F9 or F1 to F8) in any amended treatment (if require this for every amended treatments instead of any one, rare in one treatment but abundance in another, then only get too few family (try & instead of | in below))
#2) Taxa relatively abundant >0.1% in either control "whole" or amend treatment "whole" sampels

#include whole samples
View(family_taxa_proportions_tab)
filt_family_taxa_proportions_tab_SIP_mean_density_major <- family_taxa_proportions_tab[,grepl("SIP1614_AB5_F|SIP1614_AB5_whole|SIP1614_CD5_F|SIP1614_CD5_whole|SIP1614_EF5_F|SIP1614_EF5_whole|SIP1614_GH5_F|SIP1614_GH5_whole|SIP1614_IJ5_F|SIP1614_IJ5_whole|SIP1614_KL5_F|SIP1614_KL5_whole|SIP1712_AB6_F|SIP1712_AB6_whole|SIP1712_CD6_F|SIP1712_CD6_whole|SIP1712_EF6_F|SIP1712_EF6_whole|SIP1726_AB5_F|SIP1726_AB5_whole|SIP1726_CD5_F|SIP1726_CD5_whole|SIP1726_EF5_F|SIP1726_EF5_whole|SIP1726_GH5_F|SIP1726_GH5_whole",colnames(family_taxa_proportions_tab))]
dim(filt_family_taxa_proportions_tab_SIP_mean_density_major)
write.csv(filt_family_taxa_proportions_tab_SIP_mean_density_major,file="filt_family_taxa_proportions_tab_SIP_mean_density_major.csv",row.names = T) #to check some columns cannot view in R here
View(filt_family_taxa_proportions_tab_SIP_mean_density_major)
#csv insert column number for selecting columns below

countSIP1614CD<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 13:20], 1, function(x) sum(x > 0))
countSIP1614EF<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 24:31], 1, function(x) sum(x > 0))
countSIP1614IJ<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 46:53], 1, function(x) sum(x > 0))
countSIP1614KL<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 57:64], 1, function(x) sum(x > 0))
countSIP1712CD<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 79:86], 1, function(x) sum(x > 0))
countSIP1712EF<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 90:97], 1, function(x) sum(x > 0))
countSIP1726CD<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 112:119], 1, function(x) sum(x > 0))
countSIP1726EF<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 123:130], 1, function(x) sum(x > 0))
countSIP1726GH<- apply(filt_family_taxa_proportions_tab_SIP_mean_density_major[, 134:141], 1, function(x) sum(x > 0))
filt_family_taxa_proportions_tab_SIP_mean_density_widespread <- filt_family_taxa_proportions_tab_SIP_mean_density_major[which(countSIP1614CD>=3 | countSIP1614EF>=3 | countSIP1614IJ>=3 | countSIP1614KL>=3 | countSIP1712CD>=3 | countSIP1712EF>=3 | countSIP1726CD>=3 | countSIP1726EF>=3 | countSIP1726GH>=3 ),]
View(filt_family_taxa_proportions_tab_SIP_mean_density_widespread) #85 rows
dim(filt_family_taxa_proportions_tab_SIP_mean_density_widespread)
  
colnames(filt_family_taxa_proportions_tab_SIP_mean_density_widespread)
filt_family_taxa_proportions_tab_SIP_mean_density_widespread_nonrare<-filt_family_taxa_proportions_tab_SIP_mean_density_widespread[apply(filt_family_taxa_proportions_tab_SIP_mean_density_widespread,1,function(x) x[11]>0.1|x[22]>0.1|x[33]>0.1|x[44]>0.1|x[55]>0.1|x[66]>0.1|x[77]>0.1|x[88]>0.1|x[99]>0.1|x[110]>0.1|x[121]>0.1|x[132]>0.1|x[143]>0.1),]
View(filt_family_taxa_proportions_tab_SIP_mean_density_widespread_nonrare) #62 rows
#get out mean density of these widespread nonrare taxa
mean_density_major <- mean_density[row.names(mean_density) %in% row.names(filt_family_taxa_proportions_tab_SIP_mean_density_widespread_nonrare),]
View(mean_density_major)
write.csv(mean_density_major,"mean_density_major.csv")

#remove 12C and lignin 
mean_density_major_13C <- mean_density_major[,-c(1,3,5,10)]
View(mean_density_major_13C)
mean_density_major_cluster_13C <- as.matrix(mean_density_major_13C)
mean_density_major_cluster_13C.dendro <- as.dendrogram(hclust(d = dist(x = mean_density_major_cluster_13C)))
dendro_major_13C.plot <- ggdendrogram(data = mean_density_major_cluster_13C.dendro, rotate = TRUE)
dendro_major_13C.plot <- dendro_major_13C.plot + theme(axis.text.y = element_text(size = 10))
  
mean_density_major_13C$Family_Taxa<-rownames(mean_density_major_13C)
mean_density_major_13C.long <- melt(mean_density_major_13C, id = "Family_Taxa")
View(mean_density_major_13C.long)
  
mean_density_major_13C.order <- order.dendrogram(mean_density_major_cluster_13C.dendro)
mean_density_major_13C.long$Family_Taxa <- factor(x = mean_density_major_13C.long$Family_Taxa,
                                              levels = mean_density_major_13C$Family_Taxa[mean_density_major_13C.order], 
                                              ordered = TRUE)
  
heatmap_major_13C.plot <- ggplot(data = mean_density_major_13C.long, aes(x = variable, y = Family_Taxa)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colors = c("navy", "blue", "cyan", "lightcyan","yellow", "red", "red4")) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top") + theme(plot.subtitle = element_text(vjust = 1), 
                                           plot.caption = element_text(vjust = 1), 
                                           axis.text.x = element_text(angle = 10))
  
pdf("mean_density_major_13C_nolig.pdf",width=30,height=20)
grid.newpage()
print(heatmap_major_13C.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro_major_13C.plot, vp = viewport(x = 0.90, y = 0.487, width = 0.2, height = 1.02))
dev.off()
 

####Determine 13C uptake taxa,increasing with density pattern (clear pattern, not >1 up and down) or V pattern (V down not in last four fractions,increasing in at least three heavy fractions) in SIP plot (at least three fraction non-zero, not include unannotated and other taxa), also clear increasing or V (not >1 up and down, V down not in last four fractions) pattern in logratio plot, underscore with red line
View(mean_density_major_13C)
write.csv(mean_density_major_13C, file="mean_density_major_13C.csv") #label uptake ones, save as mean_density_major_13C_uptake
#remove Family_taxa, otherwise cause NA in cluster analysis

mean_density_major_13C_uptake <- mean_density_major_13C[c(2,3,14,16:18,21:26,29,31,35,36,38,41,43,44,46,52,54),-7]
colnames(mean_density_major_13C_uptake) <- c("S_TW_lysate_SPE","M_TW_lysate_SPE","M_Amino_acid","M_Syn_lysate","M_Syn_lysate_SPE","M_Syn_exudate_SPE")   #edit name to real sample name
View(mean_density_major_13C_uptake)
write.csv(mean_density_major_13C_uptake,"mean_density_major_13C_uptake.csv")

mean_density_major_13C_uptake_cluster <- as.matrix(mean_density_major_13C_uptake)
mean_density_major_13C_uptake.dendro <- as.dendrogram(hclust(d = dist(x = mean_density_major_13C_uptake_cluster)))
dendro_major_13C_uptake.plot <- ggdendrogram(data = mean_density_major_13C_uptake.dendro, rotate = TRUE)
dendro_major_13C_uptake.plot <- dendro_major_13C_uptake.plot + theme(axis.text.y = element_text(size = 18))

#Add back Family_Taxa column
mean_density_major_13C_uptake$Family_Taxa<-row.names(mean_density_major_13C_uptake)
mean_density_major_13C_uptake.long <- melt(mean_density_major_13C_uptake, id = "Family_Taxa")
View(mean_density_major_13C_uptake.long)

mean_density_major_13C_uptake.order <- order.dendrogram(mean_density_major_13C_uptake.dendro)
mean_density_major_13C_uptake.long$Family_Taxa <- factor(x = mean_density_major_13C_uptake.long$Family_Taxa,
                                                                            levels = mean_density_major_13C_uptake$Family_Taxa[mean_density_major_13C_uptake.order], 
                                                                            ordered = TRUE)

#squeeze heatmap bar
heatmap_major_13C_uptake.plot <- ggplot(data = mean_density_major_13C_uptake.long, aes(x = variable, y = Family_Taxa)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn("mean density (g/mL)",colors = c("navy", "blue", "cyan", "lightcyan","yellow", "red", "red4"),guide=guide_colorbar(reverse=TRUE)) +  #colorbar to use continous color, reverse to put dense density at bottom of legend
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "left") + theme(plot.subtitle = element_text(vjust = 1), 
                                         plot.caption = element_text(vjust = 1), 
                                         axis.text.x = element_text(angle=90,size=25,face="bold",color="black")) + 
  theme(legend.text = element_text(size = 20))  + theme(legend.title = element_text(size = 20))  + theme(legend.key.size = unit(3.5,"line"))
#dendrogram text bold
dendro_major_13C_uptake.plot <- dendro_major_13C_uptake.plot + theme(axis.text.y = element_text(size = 18,face="bold",color="black"))

pdf("mean_density_major_13C_uptake.pdf",width=24,height=20)
grid.newpage()
print(dendro_major_13C_uptake.plot , vp = viewport(x = 0.597, y = 0.58, width = 0.5, height = 0.868))
print(heatmap_major_13C_uptake.plot, vp = viewport(x = 0.2, y = 0.5, width = 0.3, height = 1.0))
dev.off()

####draw histogram for each treatment mean density across taxa

pdf("histogram_S_TW.pdf")
ggplot(data=mean_density_major_13C_uptake, aes(x=S_TW_lysate_SPE)) + 
  geom_histogram(breaks=seq(1.69, 1.77, by=0.005), 
                 col="black",
                 fill="gray", 
                 alpha = .2) + 
  labs(title="Histogram for mean density in S TW lysate PPL", x="mean density", y="Count") + 
  scale_x_continuous(limits = c(1.69,1.77), expand = c(0.01, 0),breaks=pretty_breaks(n=7)) +   #no extra space between y axis and scale
  scale_y_continuous(limits = c(0,13), expand = c(0, 0),breaks=pretty_breaks(n=5)) +theme(axis.text=element_text(colour="black"),
                                                                                          panel.grid.major = element_line(linetype = "blank"), 
                                                                                          panel.grid.minor = element_line(linetype = "blank"), 
                                                                                          panel.background = element_rect(colour="black",size=0.9,fill=NA))
dev.off()
pdf("histogram_M_TW.pdf")
ggplot(data=mean_density_major_13C_uptake, aes(x=M_TW_lysate_SPE)) + 
  geom_histogram(breaks=seq(1.69, 1.77, by=0.005), 
                 col="black",
                 fill="gray", 
                 alpha = .2) + 
  labs(title="Histogram for mean density in M TW lysate PPL", x="mean density", y="Count") + 
  scale_x_continuous(limits = c(1.69,1.77), expand = c(0.01, 0),breaks=pretty_breaks(n=7)) +   #no extra space between y axis and scale
  scale_y_continuous(limits = c(0,13), expand = c(0, 0),breaks=pretty_breaks(n=5)) +theme(axis.text=element_text(colour="black"),
                                                                                          panel.grid.major = element_line(linetype = "blank"), 
                                                                                          panel.grid.minor = element_line(linetype = "blank"), 
                                                                                          panel.background = element_rect(colour="black",size=0.9,fill=NA))
dev.off()
pdf("histogram_AA.pdf")
ggplot(data=mean_density_major_13C_uptake, aes(x=M_Amino_acid)) + 
  geom_histogram(breaks=seq(1.69, 1.77, by=0.005), 
                 col="black",
                 fill="gray", 
                 alpha = .2) + 
  labs(title="Histogram for mean density in M Amino acid", x="mean density", y="Count") + 
  scale_x_continuous(limits = c(1.69,1.77), expand = c(0.01, 0),breaks=pretty_breaks(n=7)) +   #no extra space between y axis and scale
  scale_y_continuous(limits = c(0,13), expand = c(0, 0),breaks=pretty_breaks(n=5)) +theme(axis.text=element_text(colour="black"),
                      panel.grid.major = element_line(linetype = "blank"), 
                      panel.grid.minor = element_line(linetype = "blank"), 
                      panel.background = element_rect(colour="black",size=0.9,fill=NA))
dev.off()
pdf("histogram_Syn_lysate.pdf")
ggplot(data=mean_density_major_13C_uptake, aes(x=M_Syn_lysate)) + 
  geom_histogram(breaks=seq(1.69, 1.77, by=0.005), 
                 col="black",
                 fill="gray", 
                 alpha = .2) + 
  labs(title="Histogram for mean density in M Syn lysate", x="mean density", y="Count") + 
  scale_x_continuous(limits = c(1.69,1.77), expand = c(0.01, 0),breaks=pretty_breaks(n=7)) +   #no extra space between y axis and scale
  scale_y_continuous(limits = c(0,13), expand = c(0, 0),breaks=pretty_breaks(n=5)) +theme(axis.text=element_text(colour="black"),
                                                                                          panel.grid.major = element_line(linetype = "blank"), 
                                                                                          panel.grid.minor = element_line(linetype = "blank"), 
                                                                                          panel.background = element_rect(colour="black",size=0.9,fill=NA))
dev.off()
pdf("histogram_Syn_lysate_PPL.pdf")
ggplot(data=mean_density_major_13C_uptake, aes(x=M_Syn_lysate_SPE)) + 
  geom_histogram(breaks=seq(1.69, 1.77, by=0.005), 
                 col="black",
                 fill="gray", 
                 alpha = .2) + 
  labs(title="Histogram for mean density in M Syn lysate PPL", x="mean density", y="Count") + 
  scale_x_continuous(limits = c(1.69,1.77), expand = c(0.01, 0),breaks=pretty_breaks(n=7)) +   #no extra space between y axis and scale
  scale_y_continuous(limits = c(0,13), expand = c(0, 0),breaks=pretty_breaks(n=5)) +theme(axis.text=element_text(colour="black"),
                                                                                          panel.grid.major = element_line(linetype = "blank"), 
                                                                                          panel.grid.minor = element_line(linetype = "blank"), 
                                                                                          panel.background = element_rect(colour="black",size=0.9,fill=NA))
dev.off()
pdf("histogram_Syn_exudate_PPL.pdf")
ggplot(data=mean_density_major_13C_uptake, aes(x=M_Syn_exudate_SPE)) + 
  geom_histogram(breaks=seq(1.69, 1.77, by=0.005), 
                 col="black",
                 fill="gray", 
                 alpha = .2) + 
  labs(title="Histogram for mean density in M Syn exudate PPL", x="mean density", y="Count") + 
  scale_x_continuous(limits = c(1.69,1.77), expand = c(0.01, 0),breaks=pretty_breaks(n=7)) +   #no extra space between y axis and scale
  scale_y_continuous(limits = c(0,13), expand = c(0, 0),breaks=pretty_breaks(n=5)) +theme(axis.text=element_text(colour="black"),
                                                                                          panel.grid.major = element_line(linetype = "blank"), 
                                                                                          panel.grid.minor = element_line(linetype = "blank"), 
                                                                                          panel.background = element_rect(colour="black",size=0.9,fill=NA))
dev.off()



#To validate our approach of using unamended control without 12C, we need to show statistically unamended control is not different from 12C
#I will calculate a mean density for each treatment: 13C, 12C, unamended control (not standardized RA, just using relative abundance in each treatment normalized to sum and set scale to 0-1), then compare density shift between 13C vs.12C and 13C vs. unamended control
View(filt_family_taxa_proportions_tab_SIP)
#subset SIP1614
filt_family_taxa_proportions_tab_SIP1614<-filt_family_taxa_proportions_tab_SIP[,c(1:66)]
View(filt_family_taxa_proportions_tab_SIP1614)
#subset major taxa:
#1) Taxa is widespread across amended treatments meaning they are non-zero in at least three fractions (out of F2-F9 or F1 to F8) in any amended treatment (if require this for every amended treatments instead of any one, rare in one treatment but abundance in another, then only get too few family (try & instead of | in below))
#2) Taxa relatively abundant >0.1% in either control "whole" or amend treatment "whole" sampels

#include whole sample
SIP1614_countSIP1614CD<- apply(filt_family_taxa_proportions_tab_SIP1614[, 13:20], 1, function(x) sum(x > 0))
SIP1614_countSIP1614EF<- apply(filt_family_taxa_proportions_tab_SIP1614[, 24:31], 1, function(x) sum(x > 0))
SIP1614_countSIP1614IJ<- apply(filt_family_taxa_proportions_tab_SIP1614[, 46:53], 1, function(x) sum(x > 0))
SIP1614_countSIP1614KL<- apply(filt_family_taxa_proportions_tab_SIP1614[, 57:64], 1, function(x) sum(x > 0))

filt_family_taxa_proportions_tab_SIP1614_widespread <- filt_family_taxa_proportions_tab_SIP1614[which(SIP1614_countSIP1614CD>=3 | SIP1614_countSIP1614EF>=3 | SIP1614_countSIP1614IJ>=3 | SIP1614_countSIP1614KL>=3),]
View(filt_family_taxa_proportions_tab_SIP1614_widespread) 
dim(filt_family_taxa_proportions_tab_SIP1614_widespread) #63 rows

colnames(filt_family_taxa_proportions_tab_SIP1614_widespread)
filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare<-filt_family_taxa_proportions_tab_SIP1614_widespread[apply(filt_family_taxa_proportions_tab_SIP1614_widespread,1,function(x) x[11]>0.1|x[22]>0.1|x[33]>0.1|x[44]>0.1|x[55]>0.1|x[66]>0.1),]
View(filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare) #41 rows

#subset SIP1614 only fraction 2-8, except GH F1-F8
filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare_sub<-filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare[,c(2:9,13:20,24:31,34:41,46:53,57:64)]
View(filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare_sub)

SIP1614_normalized_RA <- data.frame(matrix(ncol=48,nrow=41))
#y for row, i for column, x to control move to next 8 fractions after every sample, note GH use F1-F8
for(y in 1:41){
  for(x in 0:5){
    for (i in (1+8*x):(8+8*x)){
      SIP1614_normalized_RA[y,i]=filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare_sub[y,i]/sum(filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare_sub[y,c((1+8*x):(8+8*x))])
    }
  }
}
rownames(SIP1614_normalized_RA)<-rownames(filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare_sub)
colnames(SIP1614_normalized_RA)<-colnames(filt_family_taxa_proportions_tab_SIP1614_widespread_nonrare_sub)
View(SIP1614_normalized_RA)
sum(SIP1614_normalized_RA[1,c(17:24)])
#remove taxa rows with NaN and unannotated
SIP1614_normalized_RA_noNA<-na.omit(SIP1614_normalized_RA)
SIP1614_normalized_RA_noNA<-SIP1614_normalized_RA_noNA[-21,]
View(SIP1614_normalized_RA_noNA) #20 rows

#add in density for AB and GH control
density1614_AB<-c(1.778,1.772,1.754,1.744,1.735,1.726,1.716,1.707,1.697,1.651)
density1614_GH<-c(1.769,1.761,1.748,1.736,1.725,1.716,1.710,1.698,1.627,1.234)
density1614_AB_norm<-density1614_AB*ave_whole_density/1.730
density1614_GH_norm<-density1614_GH*ave_whole_density/1.728
#create a data frame to put in SIP1614 mean density
mean_density_SIP1614<-data.frame(matrix(ncol=6,nrow=20))
for(i in 1:20){
  mean_density_SIP1614[i,1]=SIP1614_normalized_RA_noNA[i,1]*density1614_AB_norm[2]+SIP1614_normalized_RA_noNA[i,2]*density1614_AB_norm[3]+SIP1614_normalized_RA_noNA[i,3]*density1614_AB_norm[4]+SIP1614_normalized_RA_noNA[i,4]*density1614_AB_norm[5]+SIP1614_normalized_RA_noNA[i,5]*density1614_AB_norm[6]+SIP1614_normalized_RA_noNA[i,6]*density1614_AB_norm[7]+SIP1614_normalized_RA_noNA[i,7]*density1614_AB_norm[8]+SIP1614_normalized_RA_noNA[i,8]*density1614_AB_norm[9]
  
  mean_density_SIP1614[i,2]=SIP1614_normalized_RA_noNA[i,9]*density1_norm[2]+SIP1614_normalized_RA_noNA[i,10]*density1_norm[3]+SIP1614_normalized_RA_noNA[i,11]*density1_norm[4]+SIP1614_normalized_RA_noNA[i,12]*density1_norm[5]+SIP1614_normalized_RA_noNA[i,13]*density1_norm[6]+SIP1614_normalized_RA_noNA[i,14]*density1_norm[7]+SIP1614_normalized_RA_noNA[i,15]*density1_norm[8]+SIP1614_normalized_RA_noNA[i,16]*density1_norm[9]
  
  mean_density_SIP1614[i,3]=SIP1614_normalized_RA_noNA[i,17]*density2_norm[2]+SIP1614_normalized_RA_noNA[i,18]*density2_norm[3]+SIP1614_normalized_RA_noNA[i,19]*density2_norm[4]+SIP1614_normalized_RA_noNA[i,20]*density2_norm[5]+SIP1614_normalized_RA_noNA[i,21]*density2_norm[6]+SIP1614_normalized_RA_noNA[i,22]*density2_norm[7]+SIP1614_normalized_RA_noNA[i,23]*density2_norm[8]+SIP1614_normalized_RA_noNA[i,24]*density2_norm[9]
  
  mean_density_SIP1614[i,4]=SIP1614_normalized_RA_noNA[i,25]*density1614_GH_norm[1]+SIP1614_normalized_RA_noNA[i,26]*density1614_GH_norm[2]+SIP1614_normalized_RA_noNA[i,27]*density1614_GH_norm[3]+SIP1614_normalized_RA_noNA[i,28]*density1614_GH_norm[4]+SIP1614_normalized_RA_noNA[i,29]*density1614_GH_norm[5]+SIP1614_normalized_RA_noNA[i,30]*density1614_GH_norm[6]+SIP1614_normalized_RA_noNA[i,31]*density1614_GH_norm[7]+SIP1614_normalized_RA_noNA[i,32]*density1614_GH_norm[8]
  
  mean_density_SIP1614[i,5]=SIP1614_normalized_RA_noNA[i,33]*density3_norm[2]+SIP1614_normalized_RA_noNA[i,34]*density3_norm[3]+SIP1614_normalized_RA_noNA[i,35]*density3_norm[4]+SIP1614_normalized_RA_noNA[i,36]*density3_norm[5]+SIP1614_normalized_RA_noNA[i,37]*density3_norm[6]+SIP1614_normalized_RA_noNA[i,38]*density3_norm[7]+SIP1614_normalized_RA_noNA[i,39]*density3_norm[8]+SIP1614_normalized_RA_noNA[i,40]*density3_norm[9]
  
  mean_density_SIP1614[i,6]=SIP1614_normalized_RA_noNA[i,41]*density4_norm[2]+SIP1614_normalized_RA_noNA[i,42]*density4_norm[3]+SIP1614_normalized_RA_noNA[i,43]*density4_norm[4]+SIP1614_normalized_RA_noNA[i,44]*density4_norm[5]+SIP1614_normalized_RA_noNA[i,45]*density4_norm[6]+SIP1614_normalized_RA_noNA[i,46]*density4_norm[7]+SIP1614_normalized_RA_noNA[i,47]*density4_norm[8]+SIP1614_normalized_RA_noNA[i,48]*density4_norm[9]
  
}

colnames(mean_density_SIP1614)<-c("S_Control","S_12C_TW_lys_PPL","S_13C_TW_lys_PPL","M_Control","M_12C_TW_lys_PPL","M_13C_TW_lys_PPL")
rownames(mean_density_SIP1614)<-rownames(SIP1614_normalized_RA_noNA)
View(mean_density_SIP1614)

#calculate mean density shift 13C vs. 12C, and 13C vs. control (shift to heavy in well labeled ones, but not well in partially labeled, so still mean density not a standard to identify, but can still calculate difference, show not different from 12C)
mean_density_SIP1614_shift<-data.frame(matrix(ncol=4,nrow=20))
mean_density_SIP1614_shift[,1]<-mean_density_SIP1614[,3]-mean_density_SIP1614[,1]
mean_density_SIP1614_shift[,2]<-mean_density_SIP1614[,3]-mean_density_SIP1614[,2]
mean_density_SIP1614_shift[,3]<-mean_density_SIP1614[,6]-mean_density_SIP1614[,4]
mean_density_SIP1614_shift[,4]<-mean_density_SIP1614[,6]-mean_density_SIP1614[,5]
rownames(mean_density_SIP1614_shift)<-rownames(mean_density_SIP1614)
colnames(mean_density_SIP1614_shift)<-c("S_13C_vs_control","S_13C_vs_12C","M_13C_vs_control","M_13C_vs_12C")
View(mean_density_SIP1614_shift)

#paired t test for identified incorporators
mean_density_SIP1614_shift_incorporator_S<-mean_density_SIP1614_shift[c(2,8:10,12),c(1,2)]
colnames(mean_density_SIP1614_shift_incorporator_S)<-c("13C_vs_control","13C_vs_12C")
mean_density_SIP1614_shift_incorporator_M<-mean_density_SIP1614_shift[c(10,19),c(3,4)]
colnames(mean_density_SIP1614_shift_incorporator_M)<-c("13C_vs_control","13C_vs_12C")
mean_density_SIP1614_shift_incorporator<-rbind(mean_density_SIP1614_shift_incorporator_S,mean_density_SIP1614_shift_incorporator_M)
View(mean_density_SIP1614_shift_incorporator)

t.test(mean_density_SIP1614_shift_incorporator[,1],mean_density_SIP1614_shift_incorporator[,2],paired=TRUE) #p 0.03585, not significantly different between 12C and unamended control


########Check neg control
colSums(count_tab_order[,grepl("Blank|neg|CsCl|frac|GB",colnames(count_tab_order))])
#take out blank and neg control samples to generate a separate phyloseq object,blank_whole 0 reads, take out, otherwise NAs will affect later analysis
count_tab_blank <- count_tab_order[,grepl("Blank_F|neg|CsCl|frac|GB",colnames(count_tab_order))]
write.csv(count_tab_blank,"count_tab_blank.csv")
sample_info_tab_blank <- sample_info_tab[grepl("Blank_F|neg|CsCl|frac|GB",rownames(sample_info_tab)),,drop=FALSE] #drop false to keep row name
OTU_blank = otu_table(count_tab_blank, taxa_are_rows = TRUE)
SAM_blank = sample_data(sample_info_tab_blank)
physeqhet_blank = phyloseq(OTU_blank,TAX,SAM_blank)
physeqhet_blank
  
family_counts_tab_blank <- otu_table(tax_glom(physeqhet_blank, taxrank="multi_Taxonomy")) 
dim(family_counts_tab_blank)
  
family_tax_vec_blank <- data.frame(tax_table(tax_glom(physeqhet_blank, taxrank="multi_Taxonomy")))$multi_Taxonomy # making a vector of family names to set as row names
family_tax_vec_blank 
  
rownames(family_counts_tab_blank) <- as.vector(family_tax_vec_blank)
head(family_counts_tab_blank)

#multiple rows as NA since I did not change NA to empty well in txt when reading in the tax tab
annotated_family_counts_tab_blank <- subset(family_counts_tab_blank,rownames(family_counts_tab_blank) != "NA_NA_NA")
unannotated_family_tax_counts_blank <- colSums(count_tab_blank) - colSums(annotated_family_counts_tab_blank)
unannotated_family_tax_counts_blank   #0

View(annotated_family_counts_tab_blank)
  
family_taxa_proportions_tab_blank <- apply(annotated_family_counts_tab_blank, 2, function(x) x/sum(x)*100)
colSums(family_taxa_proportions_tab_blank) 
  
dim(family_taxa_proportions_tab_blank)
# too many rows
View(family_taxa_proportions_tab_blank)

temp_filt_family_taxa_proportions_tab_blank <- data.frame(family_taxa_proportions_tab_blank[apply(family_taxa_proportions_tab_blank, 1, max) > 0.1, ]) 
# checking how many we have that were above this threshold
dim(temp_filt_family_taxa_proportions_tab_blank) # 48 rows
View(temp_filt_family_taxa_proportions_tab_blank)
  
filtered_proportions_blank <- colSums(family_taxa_proportions_tab_blank) - colSums(temp_filt_family_taxa_proportions_tab_blank)
filt_family_taxa_proportions_tab_blank <- rbind(temp_filt_family_taxa_proportions_tab_blank, "Other(<=0.1%)"=filtered_proportions_blank)  
View(filt_family_taxa_proportions_tab_blank)
 
filt_family_taxa_proportions_tab_blank_for_plot <- filt_family_taxa_proportions_tab_blank
 
filt_family_taxa_proportions_tab_blank_for_plot$Family_Taxa <- row.names(filt_family_taxa_proportions_tab_blank_for_plot)
  
filt_family_taxa_proportions_tab_blank_for_plot.g <- gather(filt_family_taxa_proportions_tab_blank_for_plot, Sample, Proportion, -Family_Taxa) #for ggplot
  
head(filt_family_taxa_proportions_tab_blank_for_plot.g)
head(filt_family_taxa_proportions_tab_blank_for_plot)

sample_info_blank_for_merge<-data.frame("Sample"=row.names(sample_info_tab_blank), "Experiment"=sample_info_tab_blank$Experiment, "Treatment"=sample_info_tab_blank$Treatment, stringsAsFactors=F)

filt_family_taxa_proportions_tab_blank_for_plot.g2 <- merge(filt_family_taxa_proportions_tab_blank_for_plot.g, sample_info_blank_for_merge,sort=F) #sort=F to keep original row order
  

pdf("blank_family_barplot.pdf",height=10,width=30)
ggplot(filt_family_taxa_proportions_tab_blank_for_plot.g2, aes(x=Sample, y=Proportion, fill=Family_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial family diversity via 16S rRNA sequences relative abundance in blanks")+
    theme(axis.text.x = element_text(angle = 15)) + theme(axis.text.x = element_text(size = 8))
dev.off()
  
  
######Check mock community,first check Emma's mock community, need to go to genus level
count_tab_amplicon_mock <- count_tab_order[,grepl("mock_amplicon",colnames(count_tab_order))]
write.csv(count_tab_amplicon_mock,"count_tab_amplicon_mock.csv")
sample_info_tab_amplicon_mock <- sample_info_tab[grepl("mock_amplicon",rownames(sample_info_tab)),,drop=FALSE] #drop false to keep row name
#create tax tab with class_order_family_genus names
tax_tab_mock<-data.frame(tax_tab) #convert to data frame first, otherwise next step won't work
tax_tab_mock$multi_Genus <- paste(tax_tab_mock$Phylum,tax_tab_mock$Class,tax_tab_mock$Order,tax_tab_mock$Family,tax_tab_mock$Genus, sep ="_")
#convert back to matrix for phyloseq object
tax_tab_mock<-as.matrix(tax_tab_mock)

OTU_amplicon_mock = otu_table(count_tab_amplicon_mock, taxa_are_rows = TRUE)
TAX_mock = tax_table(tax_tab_mock)
SAM_amplicon_mock = sample_data(sample_info_tab_amplicon_mock)
physeqhet_amplicon_mock = phyloseq(OTU_amplicon_mock,TAX_mock,SAM_amplicon_mock)
physeqhet_amplicon_mock
  
counts_tab_amplicon_mock <- otu_table(tax_glom(physeqhet_amplicon_mock, taxrank="multi_Genus")) 
dim(counts_tab_amplicon_mock)
  
tax_vec_amplicon_mock <- data.frame(tax_table(tax_glom(physeqhet_amplicon_mock, taxrank="multi_Genus")))$multi_Genus
tax_vec_amplicon_mock 
  
rownames(counts_tab_amplicon_mock) <- as.vector(tax_vec_amplicon_mock)
head(counts_tab_amplicon_mock)

annotated_counts_tab_amplicon_mock <- subset(counts_tab_amplicon_mock,rownames(counts_tab_amplicon_mock) != "NA_NA_NA_NA_NA")
unannotated_tax_counts_amplicon_mock <- colSums(count_tab_amplicon_mock) - colSums(annotated_counts_tab_amplicon_mock)
unannotated_tax_counts_amplicon_mock #none

genus_taxa_proportions_tab_amplicon_mock <- apply(counts_tab_amplicon_mock, 2, function(x) x/sum(x)*100)
colSums(genus_taxa_proportions_tab_amplicon_mock) 
  
dim(genus_taxa_proportions_tab_amplicon_mock)
View(genus_taxa_proportions_tab_amplicon_mock)
# too many rows, 472 rows
write.csv(genus_taxa_proportions_tab_amplicon_mock,"genus_taxa_proportions_tab_amplicon_mock.csv")
  
temp_filt_genus_taxa_proportions_tab_amplicon_mock <- data.frame(genus_taxa_proportions_tab_amplicon_mock[apply(genus_taxa_proportions_tab_amplicon_mock, 1, max) >  0.1, ])
# checking how many we have that were above this threshold
dim(temp_filt_genus_taxa_proportions_tab_amplicon_mock) # 22 rows
View(temp_filt_genus_taxa_proportions_tab_amplicon_mock)
  
filtered_proportions_amplicon_mock <- colSums(genus_taxa_proportions_tab_amplicon_mock) - colSums(temp_filt_genus_taxa_proportions_tab_amplicon_mock)
filt_genus_taxa_proportions_tab_amplicon_mock <- rbind(temp_filt_genus_taxa_proportions_tab_amplicon_mock, "Other(<=0.1%)"=filtered_proportions_amplicon_mock)  
View(filt_genus_taxa_proportions_tab_amplicon_mock) #23 rows
  
filt_genus_taxa_proportions_tab_amplicon_mock_for_plot <- filt_genus_taxa_proportions_tab_amplicon_mock
  
filt_genus_taxa_proportions_tab_amplicon_mock_for_plot$Genus_Taxa <- row.names(filt_genus_taxa_proportions_tab_amplicon_mock_for_plot)
  
filt_genus_taxa_proportions_tab_amplicon_mock_for_plot.g <- gather(filt_genus_taxa_proportions_tab_amplicon_mock_for_plot, Sample, Proportion, -Genus_Taxa) #for ggplot
  
head(filt_genus_taxa_proportions_tab_amplicon_mock_for_plot.g)
head(filt_genus_taxa_proportions_tab_amplicon_mock_for_plot)
  
sample_info_amplicon_mock_for_merge<-data.frame("Sample"=row.names(sample_info_tab_amplicon_mock), "Experiment"=sample_info_tab_amplicon_mock$Experiment, "Treatment"=sample_info_tab_amplicon_mock$Treatment, stringsAsFactors=F)
  
filt_genus_taxa_proportions_tab_amplicon_mock_for_plot.g2 <- merge(filt_genus_taxa_proportions_tab_amplicon_mock_for_plot.g, sample_info_amplicon_mock_for_merge,sort=F) #sort=F to keep original row order
  
Rhgcols_genus_amplicon_mock<-c("darkolivegreen3","indianred4","dark blue","darkkhaki","red","darkseagreen","darkmagenta","deepskyblue","darkorange1","blue","forestgreen","brown","aquamarine3","salmon1","green","olivedrab1","plum","darkorchid1","rosybrown1","darkslategray4","cyan","lightskyblue1","grey")

pdf("amplicon_mock_genus_barplot.pdf",width=20,height=10)
ggplot(filt_genus_taxa_proportions_tab_amplicon_mock_for_plot.g2, aes(x=Sample, y=Proportion, fill=Genus_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_genus_amplicon_mock)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial genus diversity via 16S rRNA sequences relative abundance in amplicon mock")+
    theme(axis.text.x = element_text(angle = 15)) + theme(axis.text.x = element_text(size = 8))
dev.off()


  
######Now Check earth microbiome mock community, go to genus level
count_tab_EM_mock <- count_tab_order[,grepl("mock78",colnames(count_tab_order))]
sample_info_tab_EM_mock <- sample_info_tab[grepl("mock78",rownames(sample_info_tab)),,drop=FALSE] #drop false to keep row name
OTU_EM_mock = otu_table(count_tab_EM_mock, taxa_are_rows = TRUE)
SAM_EM_mock = sample_data(sample_info_tab_EM_mock)
physeqhet_EM_mock = phyloseq(OTU_EM_mock,TAX_mock,SAM_EM_mock)
physeqhet_EM_mock

counts_tab_EM_mock <- otu_table(tax_glom(physeqhet_EM_mock, taxrank="multi_Genus")) 
dim(counts_tab_EM_mock) 

tax_vec_EM_mock <- data.frame(tax_table(tax_glom(physeqhet_EM_mock, taxrank="multi_Genus")))$multi_Genus
tax_vec_EM_mock 
  
rownames(counts_tab_EM_mock) <- as.vector(tax_vec_EM_mock)
head(counts_tab_EM_mock)

annotated_counts_tab_EM_mock <- subset(counts_tab_EM_mock,rownames(counts_tab_EM_mock) != "NA_NA_NA_NA_NA")
unannotated_tax_counts_EM_mock <- colSums(count_tab_EM_mock) - colSums(annotated_counts_tab_EM_mock)
unannotated_tax_counts_EM_mock #none 

genus_taxa_proportions_tab_EM_mock <- apply(counts_tab_EM_mock, 2, function(x) x/sum(x)*100)
colSums(genus_taxa_proportions_tab_EM_mock) 
  
dim(genus_taxa_proportions_tab_EM_mock)
# too many rows,472 rows
  
temp_filt_genus_taxa_proportions_tab_EM_mock <- data.frame(genus_taxa_proportions_tab_EM_mock[apply(genus_taxa_proportions_tab_EM_mock, 1, max) >  0.1, ])
# checking how many we have that were above this threshold
dim(temp_filt_genus_taxa_proportions_tab_EM_mock) # 19 rows
View(temp_filt_genus_taxa_proportions_tab_EM_mock)
  
filtered_proportions_EM_mock <- colSums(genus_taxa_proportions_tab_EM_mock) - colSums(temp_filt_genus_taxa_proportions_tab_EM_mock)
filt_genus_taxa_proportions_tab_EM_mock <- rbind(temp_filt_genus_taxa_proportions_tab_EM_mock, "Other(<=0.1%)"=filtered_proportions_EM_mock)  
View(filt_genus_taxa_proportions_tab_EM_mock) #20 rows
  
filt_genus_taxa_proportions_tab_EM_mock_for_plot <- filt_genus_taxa_proportions_tab_EM_mock
  
filt_genus_taxa_proportions_tab_EM_mock_for_plot$Genus_Taxa <- row.names(filt_genus_taxa_proportions_tab_EM_mock_for_plot)
  
filt_genus_taxa_proportions_tab_EM_mock_for_plot.g <- gather(filt_genus_taxa_proportions_tab_EM_mock_for_plot, Sample, Proportion, -Genus_Taxa) #for ggplot
  
head(filt_genus_taxa_proportions_tab_EM_mock_for_plot.g)
head(filt_genus_taxa_proportions_tab_EM_mock_for_plot)
  
sample_info_EM_mock_for_merge<-data.frame("Sample"=row.names(sample_info_tab_EM_mock), "Experiment"=sample_info_tab_EM_mock$Experiment, "Treatment"=sample_info_tab_EM_mock$Treatment, stringsAsFactors=F)
  
filt_genus_taxa_proportions_tab_EM_mock_for_plot.g2 <- merge(filt_genus_taxa_proportions_tab_EM_mock_for_plot.g, sample_info_EM_mock_for_merge,sort=F) #sort=F to keep original row order
  
Rhgcols_genus_EM_mock<-c("darkolivegreen3","indianred4","dark blue","darkkhaki","red","darkseagreen","darkmagenta","deepskyblue","darkorange1","blue","forestgreen","brown","aquamarine3","salmon1","green","olivedrab1","plum","darkorchid1","lightskyblue1","grey")
pdf("EM_mock_genus_barplot.pdf",width=15,height=10)
ggplot(filt_genus_taxa_proportions_tab_EM_mock_for_plot.g2, aes(x=Sample, y=Proportion, fill=Genus_Taxa))  +
    geom_bar(width=0.6, stat="identity") +
    scale_fill_manual(values =Rhgcols_genus_EM_mock)+
    theme_bw()+ labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="Bacterial genus diversity via 16S rRNA sequences relative abundance in EM mock")+
    theme(axis.text.x = element_text(angle = 15)) + theme(axis.text.x = element_text(size = 8))
dev.off()

