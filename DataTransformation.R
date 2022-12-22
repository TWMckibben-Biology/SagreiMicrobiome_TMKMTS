###Sagrei Gut Microbiome Project###
##Taylor McKibben, Kaitlyn Murphy, and Tonia Schwartz##
###note to self: add in t.tests for phenotype data###

#If using NCBI database, get an api key from ncbi.gov and set it
#edit_r_environ(scope = "user")
#In the new pane enter ENTREZ_KEY=PutYourAPIKeyHere
#Save that code and restart R


#Set working directory, install packages, and load them.
##Due to version differences, Phyloseq has a different installation method than usual

rm()
setwd("C:/Users/taylo/Box/TS_Lab_Microbiome/2021_Taylor_Anole_Development_Microbiome/Data/DataAnalyses/R")
getwd()


#install these packages
install.packages('lme4')
install.packages('nlme')
install.packages('emmeans')
install.packages('taxize')
install.packages('broom')
install.packages('knitr')
install.packages('breakaway')
install.packages('eulerr')
devtools::install_github('microsud/microbiomeutilities')
devtools::install_github(repo = "malucalle/selbal")
devtools::install_github("benjjneb/dada2", ref="v1.16")
.cran_packages <- c("tidyverse", "cowplot", "picante", "vegan", "HMP", "dendextend", "rms", "devtools")
.bioc_packages <- c("phyloseq", "DESeq2", "microbiome", "metagenomeSeq", "ALDEx2", "PERFect")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(.bioc_packages, version = "3.9")

install.packages('phyloseq')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("PERFect")

#call on these libraries
library(lme4)
library(nlme)
library(vegan)
library(devtools)
library(tidyverse)
library(emmeans)
library(taxize)
library(broom)
library(phyloseq)
library(DESeq2)
library(microbiome)
library(picante)
library(ALDEx2)
library(metagenomeSeq)
library(HMP)
library(dendextend)
library(rms)
library(cowplot)
library(PERFect)
library(dada2)
library(breakaway)
library(selbal)
library(eulerr)
library(microbiomeutilities)

#####Data Manipulation######

#Step 1: Load data
###I'm loading in my sequencing run data, MB_040722_epi2me
MB_040722_epi2me <- read.csv("C:/MSOfficeFiles/Excel/MicrobiomeProject/MB_040722_epi2me.csv")

#Step 2: Reduce dataset to just barcodes and species
MB_BarcodesandSpecies = data.frame(MB_040722_epi2me$barcode, MB_040722_epi2me$species)

#Step 3: Rename columns and sort data by ascending order in barcodes
colnames(MB_BarcodesandSpecies) = c("barcode", "species")
MB_BAS_Sorted = data.frame(MB_BarcodesandSpecies[order(MB_BarcodesandSpecies$barcode),])
rm(MB_BarcodesandSpecies)

#Step 4: Gain a frequency by species by barcode count and rearrange the array
MB_SpeciesBarcode = data.frame(table(MB_BAS_Sorted$species, MB_BAS_Sorted$barcode)) 
rm(MB_BAS_Sorted)
MB_Whatnow = data.frame(MB_SpeciesBarcode %>%
  gather(key, value, Freq) %>%
  spread(Var2, value)) 
rm(MB_SpeciesBarcode)
MB_SpeciesPerBarcode <- data.frame(subset (MB_Whatnow, select = -key)) 
rm(MB_Whatnow)
rownames(MB_SpeciesPerBarcode) = MB_SpeciesPerBarcode$Var1


#I would highly recommend saving this csv as it is finally small
#enough to open in excel without melting your computer
write.csv(MB_SpeciesPerBarcode, "C:/MSOfficeFiles/Excel/Microbiome_Clean/SpeciesPerBarcode.csv")

######Data Analysis#####
#For metrics such as Shannon's, I am creating a file that I can add
#to an independent "phenotypic data" file that contains other data about my
#samples such as mass and SVL. The file is organized such that each row is
#an individual and each column is a phenotype.

#Filtering data

##Test code for PERFect
###Format data table to fit package
MB_SPB_PERF = MB_SpeciesPerBarcode[-1,-1]
specname = unique(row.names(MB_SPB_PERF))
MB_SPB_PERF = MB_SPB_PERF[,-24]
MB_SPB_PERF = data.frame(t(MB_SPB_PERF))
colnames(MB_SPB_PERF) = specname

#PERFect_sim
#Returns a list of significant taxa
MB_PERF_FILT = PERFect_sim(MB_SPB_PERF)
MB_SigTax = colnames(MB_PERF_FILT$filtX)

#Filter data using PERFect results
MB_Filt = data.frame(t(MB_SpeciesPerBarcode))
MB_Filtrow = grep(paste(MB_SigTax, collapse = "|"), MB_SpeciesPerBarcode$Var1)
MB_Filt = data.frame(MB_SpeciesPerBarcode[MB_Filtrow, -c(1,26)])
write.csv(MB_Filt, "C:/MSOfficeFiles/Excel/Microbiome_Clean/MB_Filt.csv")
rm(MB_Filtrow)
rm(MB_SigTax)

#Phyloseq
#import if you havent already("C:/MSOfficeFiles/Excel/MicrobiomeProject/MB_040722_epi2me")
#Microbiome Statistical Analysis
#taken from https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
#Modified to use ONT minKNOW/EPI2ME output, data used in the tutorial was a QIIME2 rds file
#Many of the following steps prior to data analysis are transformation to fit my data into the pipeline

  
#Create a taxonomy df using taxize
#create a vector of species name using our filtered dataset
MB_Uspec = row.names(MB_Filt)
#Get taxonomy hierarchy using unique taxid for NCBI-Taxize pkg
MB_Hier = classification(MB_Uspec, db = 'ncbi', return_id = TRUE)
rm(MB_040722_epi2me)

#create a tax hierarchy name vector
hierarchy = c("no rank", "superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species", "strain")
#create a df using our hierarchy vector and species vector
MB_Taxonomy = data.frame(setNames(data.frame(matrix(ncol = length(hierarchy), nrow = length(MB_Uspec))), hierarchy))
MB_Hier = rbind(MB_Hier)
MB_Taxonomy = cbind(MB_Taxonomy, MB_Uspec)
names(MB_Taxonomy)[names(MB_Taxonomy) == 'MB_Uspec'] = 'query'
MB_Taxonomy = data.frame(MB_Taxonomy[, 11])
row.names(MB_Taxonomy) = MB_Uspec
colnames(MB_Taxonomy) = c("query")

#subset out each taxonomy rank and prep it to merge
norank = subset(MB_Hier, rank == "no rank")
norank = data.frame(norank[, -2:-3])
colnames(norank) = c("norank", "query")

superkingdom = subset(MB_Hier, rank == "superkingdom")
superkingdom = data.frame(superkingdom[, -2:-3])
colnames(superkingdom) = c("superkingdom", "query")

clade = subset(MB_Hier, rank == "clade")
clade = data.frame(clade[, -2:-3])
colnames(clade) = c("clade", "query")

phylum = subset(MB_Hier, rank == "phylum")
phylum = data.frame(phylum[, -2:-3])
colnames(phylum) = c("phylum", "query")

class = subset(MB_Hier, rank == "class")
class = data.frame(class[, -2:-3])
colnames(class) = c("class", "query")

order = subset(MB_Hier, rank == "order")
order = data.frame(order[, -2:-3])
colnames(order) = c("order", "query")

family = subset(MB_Hier, rank == "family")
family = data.frame(family[, -2:-3])
colnames(family) = c("family", "query")

genus = subset(MB_Hier, rank == "genus")
genus = data.frame(genus[, -2:-3])
colnames(genus) = c("genus", "query")

species = subset(MB_Hier, rank == "species")
species = data.frame(species[, -2:-3])
colnames(species) = c("species", "query")

strain = subset(MB_Hier, rank == "strain")
strain = data.frame(strain[, -2:-3])
colnames(strain) = c("strain", "query")

#Merge all of your rank df with the taxonomy df

MB_Taxonomy = merge(MB_Taxonomy, norank, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, superkingdom, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, clade, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, phylum, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, class, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, order, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, family, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, genus, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, species, by = "query", all.x = TRUE)
MB_Taxonomy = merge(MB_Taxonomy, strain, by = "query", all.x = TRUE)

#Remove species with no data and duplicated species before cleaning up
##For phylum firmicutes, clades CFB and Bacteroidetes are inconsistently used. Check your own data.
###I will convert them to the same name or phyloseq wont be able to glom
MB_Taxonomy$clade = gsub("Bacteroidetes/Chlorobi group", "FCB group", MB_Taxonomy$clade)
#If you want analysis on any level lower than phylum, you will need to repeat this down to your preferred level
#Special note: Phyloseq gloms by checking that all of your samples at a given level, "phylum" in this case,
#have identical upstream taxonomy. If one phylum has differing names upstream, the function will not work.
#For example, proteobacteria have no clade classification, but as long as we replace the NA values with the same text, ie "clade, the function will still work
#If your data isn't cleaned up, row.names wont work and you cant proceed.
MB_Taxonomy = data.frame(MB_Taxonomy[!duplicated(MB_Taxonomy$query),])
row.names(MB_Taxonomy) = MB_Taxonomy$query

#Now we have to match our otu table to the taxonomy table
#And we need to make our sample_data object

#Finally make ps objects out of our shiny new dataframes

#Tax table
MB_Taxonomy$norank <- MB_Taxonomy$norank %>% replace_na('norank')
MB_Taxonomy$superkingdom <- MB_Taxonomy$superkingdom %>% replace_na('superkingdom')
MB_Taxonomy$clade <- MB_Taxonomy$clade %>% replace_na('clade')
MB_Taxonomy$phylum <- MB_Taxonomy$phylum %>% replace_na('phylum')
MB_Taxonomy$class <- MB_Taxonomy$class %>% replace_na('class')
MB_Taxonomy$order <- MB_Taxonomy$order %>% replace_na('order')
MB_Taxonomy$family <- MB_Taxonomy$family %>% replace_na('family')
MB_Taxonomy$genus <- MB_Taxonomy$genus %>% replace_na('genus')
MB_Taxonomy$species <- MB_Taxonomy$species %>% replace_na('species')
MB_Taxonomy$strain <- MB_Taxonomy$strain %>% replace_na('strain')
MB_Taxonomy = data.frame(MB_Taxonomy[,-1:-2])
#OTU table-Needs row names before Tax becomes a matrix
MB_Filt = subset(MB_Filt, row.names(MB_Filt)%in%row.names(MB_Taxonomy))
OTUmat = as.matrix(MB_Filt)
OTU = otu_table(OTUmat, taxa_are_rows = TRUE)
#Tax Table
Taxmat = as.matrix(MB_Taxonomy)
TAX = tax_table(Taxmat)

#Sample data table
#I am loading my MB_Pheno_Data_WMaternal file
mydat = read.csv(file.choose())
mydat = sample_data(mydat)
mydat$BarcodeID = formatC(mydat$BarcodeID, width = 2, format = "d", flag = "0")
barcodes = as.character(mydat$BarcodeID)
barcodes = paste("barcode", barcodes, sep = "")
row.names(mydat) = barcodes
mydat = mydat[,-3]
SAM = sample_data(mydat)

#Now we can start phyloseq analysis
MBphy = phyloseq(OTU,TAX,SAM)

sort(phyloseq::sample_sums(MBphy))

MBphyfw = subset_samples(MBphy, Treatment != "Maternal")

MBphy %>% 
  sample_data %>%
  dplyr::count(Treatment)

MBphyfw %>% 
  sample_data %>%
  dplyr::count(Treatment)

(MBphy <- phyloseq::subset_samples(MBphy, phyloseq::sample_sums(MBphy) > 5000)) 
(MBphy <- phyloseq::prune_taxa(phyloseq::taxa_sums(MBphy) > 0, MBphy))

(MBphyfw <- phyloseq::subset_samples(MBphyfw, phyloseq::sample_sums(MBphyfw) > 5000)) 
(MBphyfw <- phyloseq::prune_taxa(phyloseq::taxa_sums(MBphyfw) > 0, MBphyfw)) 

#Relative Abundance
ps_rel_abund = phyloseq::transform_sample_counts(MBphy, function(x){x / sum(x)})
phyloseq::otu_table(MBphy)
phyloseq::otu_table(ps_rel_abund)

ps_rel_abundfw = phyloseq::transform_sample_counts(MBphyfw, function(x){x / sum(x)})
phyloseq::otu_table(MBphyfw)
phyloseq::otu_table(ps_rel_abundfw)

#Relative Abundance Stacked Bar Test
png(filename = "RelAbundAll.png", width = 140, height = 180, units = "mm", res = 1000, pointsize = 12)

phyloseq::plot_bar(ps_rel_abund, fill = "phylum") +
  geom_bar(aes(color = phylum, fill = phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

png(filename = "RelAbundFW.png", width = 140, height = 180, units = "mm", res = 1000, pointsize = 12)

phyloseq::plot_bar(ps_rel_abundfw, fill = "phylum") +
  geom_bar(aes(color = phylum, fill = phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Treatment, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

rank_names(physeq = MBphy)
is.na(MBphy@tax_table[,5])

#Faceted Box Plots
#aggs data into phylum level
ps_phylum = phyloseq::tax_glom(MBphy, "phylum")
(ps_phylum <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_phylum) > 0, ps_phylum)) 
(ps_phylum <- phyloseq::prune_taxa(phyloseq::taxa_names(ps_phylum) != "phylum", ps_phylum)) 

ps_phylumfw = phyloseq::tax_glom(MBphyfw, "phylum")
(ps_phylumfw <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_phylumfw) > 0, ps_phylumfw)) 
(ps_phylumfw <- phyloseq::prune_taxa(phyloseq::taxa_names(ps_phylumfw) != "phylum", ps_phylumfw)) 

#should rename rows based on phyla
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "phylum"]
phyloseq::taxa_names(ps_phylumfw) <- phyloseq::tax_table(ps_phylumfw)[, "phylum"]

#Sanity check that those steps worked
phyloseq::otu_table(ps_phylum)
(ps_phylum <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_phylum) > 100, ps_phylum)) 
taxa_sums(ps_phylum@otu_table)

phyloseq::otu_table(ps_phylumfw)
(ps_phylumfw <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_phylumfw) > 100, ps_phylumfw)) 
taxa_sums(ps_phylumfw@otu_table)

png(filename = "facboxplot.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Log 10 Abundance\n", ) +
  facet_wrap(~ OTU, scales = "free") +
  scale_y_log10()
  theme(strip.text = element_text(size = 8))

dev.off()

png(filename = "facboxplotfw.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

phyloseq::psmelt(ps_phylumfw) %>%
  ggplot(data = ., aes(x = Treatment, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Log 10 Abundance\n", ) +
  facet_wrap(~ OTU, scales = "free") +
  scale_y_log10()
theme(strip.text = element_text(size = 8))

dev.off()

#Dirichlet-Multinomial Distribution
control <- phyloseq::subset_samples(ps_phylum, Treatment == "Control")
diverse <- phyloseq::subset_samples(ps_phylum, Treatment == "Diverse")
#Output OTU tables
control_otu <- data.frame(phyloseq::otu_table(control))
div_otu <- data.frame(phyloseq::otu_table(diverse))
#Group rare phyla
control_otu <- control_otu %>%
  t(.) %>%
  as.data.frame(.)
  #mutate(Other = Cyanobacteria + Euryarchaeota + Tenericutes + Verrucomicrobia + Fusobacteria) %>%
  #dplyr::select(-Cyanobacteria, -Euryarchaeota, -Tenericutes, -Verrucomicrobia, -Fusobacteria)
div_otu <- div_otu %>%
  t(.) %>%
  as.data.frame(.)
  #mutate(Other = Cyanobacteria + Euryarchaeota + Tenericutes + Verrucomicrobia + Fusobacteria) %>%
  #dplyr::select(-Cyanobacteria, -Euryarchaeota, -Tenericutes, -Verrucomicrobia, -Fusobacteria)
#HMP test
group_data <- list(control_otu, div_otu)
(xdc <- HMP::Xdc.sevsample(group_data))         

#Bray-curtis dissimilarity-Hiararchical Clustering
#Val 0 = samples share all taxa, Val 1 = samples share no taxa
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abundfw))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
as.matrix(bc_dist)[1:5, 1:5]

ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(ps_rel_abundfw))
colorCode <- c(Control = "red", `Diverse` = "blue")
labels_colors(ward) <- colorCode[meta$Treatment][order.dendrogram(ward)]
#Plot
png(filename = "wardfw.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

plot(ward)

dev.off()

###Breakaway method###

#This is an ALTERNATE measure of Alpha-diversity
#Do not include both Shannon and Breakaway
#Or do, and have a comparison discussion
ba_adiv <- breakaway::breakaway(MBphyfw)
ba_adiv[1]
#Plot estimates
png(filename = "breakaway.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

plot(ba_adiv, MBphyfw, color = "Treatment")  

dev.off()

#Examine models
summary(ba_adiv) %>%
  add_column("SampleNames" = MBphyfw %>% otu_table %>% sample_names)  
#Test for group differnce
bt <- breakaway::betta(summary(ba_adiv)$estimate,
                       summary(ba_adiv)$error,
                       make_design_matrix(MBphyfw, "Treatment"))
bt$table    

###Shannons anyways for good measure
png(filename = "Shannon.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(MBphyfw),
                         "observed" = phyloseq::estimate_richness(MBphyfw, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")

dev.off()

#Subsample for rare reads/group differences
(ps_rare <- phyloseq::rarefy_even_depth(ps_phylumfw, replace = FALSE)) 
head(phyloseq::sample_sums(ps_rare))
adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(ps_rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps_rare, measures = "Shannon"),
  #"PD" = picante::pd(samp = data.frame(t(data.frame(phyloseq::otu_table(ps_rare)))), tree = phyloseq::phy_tree(ps_rare))[, 1],
  "Group" = phyloseq::sample_data(ps_rare)$Treatment)
head(adiv)
#Plot
png(filename = "adiv.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon"))) %>%
  ggplot(aes(x = Group, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Group), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free") +
  theme(legend.position="none")

dev.off()

#Wilcoxon Diversity-Location
adiv %>%
  group_by(Group) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon))
#Shannon test
wilcox.test(Observed ~ Shannon, data = adiv, exact = FALSE, conf.int = TRUE)

#Beta-diversity
#Aitchison Distance
(ps_clr <- microbiome::transform(ps_phylumfw, "clr")) 
phyloseq::otu_table(ps_phylumfw)
phyloseq::otu_table(ps_clr)

#PCA plot-data generation
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")

#Plot scree plot
png(filename = "ScreePCA.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
head(ord_clr$CA$eig)  
sapply(ord_clr$CA$eig, function(x) x / sum(ord_clr$CA$eig))

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)

#PCA plotting
phyloseq::plot_ordination(ps_phylumfw, ord_clr, type="samples", color="Treatment") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = Treatment), linetype = 2)

dev.off()

#PCA PERMANOVA
#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$Treatment)
#Dispersion
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$Treatment)
dispr
#Dispersion Plot
png(filename = "Dispersion.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr)

dev.off()

#UniFrac Method
#Generate distances
png(filename = "UniFrac.png", width = 560, height = 560, units = "mm", res = 1000, pointsize = 12)

ord_unifrac <- ordinate(ps_rare, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(ps_rare, method = "PCoA", distance = "unifrac")   

dev.off()

#Plot ordinations
a <- plot_ordination(ps_rare, ord_unifrac, color = "Treatment") + geom_point(size = 2)
b <- plot_ordination(ps_rare, ord_unifrac_un, color = "Treatment") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))

#Finding potential causal bacteria
#Generate data.frame with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$Treatment <- phyloseq::sample_data(ps_clr)$Treatment
#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ Treatment, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ Treatment, data = df)$p.value
}
#Create nested data frames by OTU and loop over each using map 
wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -Treatment) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       
#Show results
head(wilcox_results)
head(wilcox_results$data[[1]])
wilcox_results$wilcox_test[[1]]
wilcox_results$p_value[[1]]
#Unnesting
wilcox_results <- wilcox_results %>%
  dplyr::select(OTU, p_value) %>%
  unnest()
head(wilcox_results)  
#Adding taxonomic labels
taxa_info <- data.frame(tax_table(ps_clr))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")
#Computing FDR corrected p-values
wilcox_results <- wilcox_results %>%
  full_join(taxa_info) %>%
  arrange(p_value) %>%
  mutate(BH_FDR = p.adjust(p_value, "BH")) %>%
  filter(BH_FDR < 0.05) %>%
  dplyr::select(OTU, p_value, BH_FDR, everything())
print.data.frame(wilcox_results)  

#Run ALDEx2- ANOVA-like Differential Expression (ALDEx2 pkg)
aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(MBphy)), phyloseq::sample_data(MBphy)$Treatment, test="t", effect = TRUE, denom="iqlr")
#ALDE plot
#Plot effect sizes
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#Clean presentation
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
sig_aldex2

#Predictive modeling
#Generate data.frame
clr_pcs <- data.frame(
  "pc1" = ord_clr$CA$u[,1],
  "pc2" = ord_clr$CA$u[,2],
  "pc3" = ord_clr$CA$u[,3],
  "Group" = phyloseq::sample_data(ps_clr)$Group
)
clr_pcs$Group_num <- ifelse(clr_pcs$Group == "Control", 0, 1)
head(clr_pcs)
#Specify a datadist object (for rms)
dd <- datadist(clr_pcs)
options(datadist = "dd")
#Plot the unconditional associations
a <- ggplot(clr_pcs, aes(x = pc1, y = Group_num)) +
  Hmisc::histSpikeg(Group_num ~ pc1, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC1", y = "Pr(Diverse)\n")
b <- ggplot(clr_pcs, aes(x = pc2, y = Group_num)) +
  Hmisc::histSpikeg(Group_num ~ pc2, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC2", y = "Pr(Diverse)\n")
c <- ggplot(clr_pcs, aes(x = pc3, y = Group_num)) +
  Hmisc::histSpikeg(Group_num ~ pc3, lowess = TRUE, data = clr_pcs) +
  labs(x = "\nPC3", y = "Pr(Diverse)\n")
cowplot::plot_grid(a, b, c, nrow = 2, ncol = 2, scale = .9, labels = "AUTO")
#Fit full model with splines (5 knots each)
m1 <- rms::lrm(Group_num ~ rcs(pc1, 5) + rcs(pc2, 5) + rcs(pc3, 5), data = clr_pcs, x = TRUE, y = TRUE)

#Grid search for penalties
pentrace(m1, list(simple = c(0, 1, 2), nonlinear = c(0, 100, 200)))
pen_m1 <- update(m1, penalty = list(simple = 1, nonlinear = 200))
pen_m1

#Plot log odds
ggplot(Predict(pen_m1))

#Obtain optimism corrected estimates
(val <- rms::validate(pen_m1))
#Compute corrected c-statistic
(c_opt_corr <- 0.5 * (val[1, 5] + 1))
#Plot calibration
cal <- rms::calibrate(pen_m1, B = 200)
plot(cal)
#Output pred. probs
head(predict(pen_m1, type ="fitted"))

#Heatmap
png(file = "phylumheatmap.png", width = 380, height = 190, units = "mm", res = 1000, pointsize = 12)

c = plot_heatmap(MBphy, taxa.label = "phylum")

plot_grid(c, labels = "Phylum", align = "hv", scale = 0.9)

dev.off()

#Shared Core Microbiota -microbiome, microbiomeutilities, and eulerr pkgs
#This should use phyloseq objects, I think

table(meta(ps_phylum)$Treatment, useNA = "always")
pseq.rel <- microbiome::transform(ps_phylum, "compositional")
groups = unique(as.character(meta(pseq.rel)$Treatment))
print(groups)

list_core <- c() # an empty object to store information

for (n in groups){
  
  ps.sub <- subset_samples(pseq.rel, groups == n) 
  
  core_m <- core_members(ps.sub, 
                         detection = 0,  
                         prevalence = 0)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) 
  list_core[[n]] <- core_m 
  print(list_core)
}
mycols <- c(Diverse="#d6e2e9", Control="#cbf3f0", Maternal="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)
