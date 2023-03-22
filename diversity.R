# Importing data 
library(phyloseq)
library(ggplot2)
library(ape)

#import OTU table (from .csv file)
Install.packages(“here”)
Library(here)
in_path <- here("phyloseq", "otu_table3.txt")
otu <- read.table(file = in_path, header = TRUE, check.names=FALSE)

rownames(otu)<-otu$OTUID
otu_table<-otu_table[-c(1)]
dim(otu_table)
otu <- as.matrix(otu)

taxonomy <- read.csv("taxonomy.csv", sep = ",", row.names = 1)
taxonomy <- as.matrix(taxonomy)

#import as csv – changed rownames to sample IDs 
metadata4<-read.csv(file.choose(Fungal_metadata_endohisto2.tsv), header=TRUE)
rownames(metadata4)<-metadata4$sample.id
metadata4<-metadata4[-c(1)]
META<-sample_data(metadata4)
phy_tree <- read_tree("tree.nwk")


taxmat = matrix(nrow = nrow(taxonomy), ncol = 7)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

library(stringr)
split<-str_split_fixed(taxonomy$Taxon, ";", 7)
rownames(split)<-rownames(taxonomy)
colnames(split)<-colnames(taxmat)
taxonomy2<-as.matrix(split)

OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy2)
META <- sample_data(metadata4)

# Sanity checks for consistent OTU names
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

# Same sample names
sample_names(OTU)
sample_names(META)

# Finally merge!
ps <- phyloseq(OTU, TAX, META, phy_tree)
ps



#ALPHA DIVERSITY 
#removed zero counts 
pruned = prune_taxa(taxa_sums(ps) > 0, ps)
#subsetted to Endoscopic Activity vs Remission 
lod<-subset_samples(ps.genus, Remission.State!= "NA")
sample_data(lod)$Remission.State<-as.character(sample_data(lod)$Remission.State)
#plot alpha diversity
plot_richness(lod, x="Remission.State", measures=c("Observed", "Shannon"), color="Remission.State")+geom_point(size=5, alpha=0.7) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_discrete(c("0","1")) + scale_color_manual(values = c("navyblue", "firebrick")) + geom_boxplot(width=0.5)+ geom_jitter(height = 0, width = 0.25, size=3)

results = estimate_richness(lod, measures = c(“Shannon”, “Observed”)
d = sample_data(lod)
remission = results[d[,'Remission.State'] == '0',]
activity = results[d[,'Remission.State'] == '1',]
wilcox.test(remission,activity)
summary(remission)
summary(activity)

#Shannon effective numbers
H2 <- exp(estimate_richness(lod, measures="Shannon")
d = sample_data(lod)
remission = H2[d[,'Remission.State'] == '0',]
activity = H2[d[,'Remission.State'] == '1',]
wilcox.test(remission,activity)
summary(remission)
summary(activity)




#BETA DIVERSITY 
pruned = prune_samples(sample_sums(ps) > 0, ps)
sample_sums(pruned)
lod<-subset_samples(pruned, Remission.State!= "NA")
sample_data(lod)$Remission.State<-as.character(sample_data(lod)$Remission.State)

#For weighted unifrac + NMDS
wunifrac_dist = phyloseq::distance(lod, method="unifrac", weighted=T)
NMDS.lod <- ordinate(lod, "NMDS", wunifrac_dist, weakties=FALSE)
plot_ordination(lod, NMDS.lod, color="Remission.State") + theme_bw() +theme(aspect.ratio=1)+ theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#For unweighted unifrac + PCOA
#Unifrac
wunifrac_dist = phyloseq::distance(lod, method="unifrac", weighted=F)
ordination = ordinate(lod, method="PCoA", distance=wunifrac_dist)
plot_ordination(lod, ordination, color="Remission.State") + theme_bw() +theme(aspect.ratio=1)+ theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#Plotted Beta Diversity
plot_ordination(lod, ordination, color="Remission.State") + theme_bw() +theme(aspect.ratio=1)+ theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

results2<-adonis2(wunifrac_dist ~ sample_data(lod.rarefied)$Remission.State)
summary(results2)
write_tsv(results2, "betadivadonis.tsv")

#FRACTIONAL PREVALENCE OF PHYLUM AND GENUS
ps.filter<-phyloseq_filter_prevalence(ps, prev.trh = 0.05, abund.trh = NULL,                                       threshold_condition = "OR", abund.type = "total")
ps.genus = tax_glom(ps.filter, taxrank="Genus", NArm=FALSE)
prevdf2 <- apply(X = otu_table(ps.genus), MARGIN = ifelse(taxa_are_rows(ps.genus), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
prevdf2 <- data.frame(Prevalence = prevdf2, TotalAbundance = taxa_sums(ps.genus), tax_table(ps.genus))
plyr::ddply(prevdf2, "Genus", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
ggplot(prevdf2, aes(TotalAbundance, Prevalence / nsamples(ps.genus),color=Genus)) +geom_point(size = 7, alpha = 0.7) +geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +scale_x_log10() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]")

#DIFFERENTIAL ABUNDANCE ANALYSIS

#adding 1 pseudocount to all OTU counts and remaking physeq object 
Otu2<-otu+1
OTU <- otu_table(realotu, taxa_are_rows = TRUE)
ps2<-phyloseq(OTU,TAX,META,phy_tree)

#removed rare organisms
ps3<-phyloseq_filter_taxa_tot_fraction(ps2, frac = 0.01)
ps.genus = tax_glom(ps3, taxrank="Genus", NArm=FALSE)

head(sample_data(ps)$Remission.State, 25)
ps = subset_samples(ps, Remission.State != "NA")
head(sample_data(ps)$Remission.State, 25)
sample_data(ps)$Remission.State<-as.character(sample_data(ps)$Remission.State)

diagdds = phyloseq_to_deseq2(ps, ~ Remission.State)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

write_tsv(sigtab, “deseqresults.tsv”)



#Adjusting for confounders in relative abundance

ps.na = subset_samples(ps, Remission.State != "NA")
ps.combo = subset_samples(ps.na, biologics!=”NA”)
ps.combo2 = tax_glom(ps.combo, taxrank="Genus", NArm=FALSE)
diagdds = phyloseq_to_deseq2(ps.combo2, design=~biologics+Genderscore+Agescore+Remission.State) 

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.combo2)[rownames(sigtab), ], "matrix"))
head(sigtab)

x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#RAW ABUNDANCE BOX PLOTS
ps <- phyloseq(OTU, TAX, META, phy_tree)
ps4 = tax_glom(ps, taxrank="Genus", NArm=FALSE)
ps.subset<-subset_samples(ps4, Remission.State != "NA")
sample_data(ps.subset)$Remission.State<-as.character(sample_data(ps.subset)$Remission.State)
GP.candida <- subset_taxa(ps.subset, Genus=="g__Candida")
p<-phyloseq::otu_table(GP.candida)[1:1, 1:44]
p<-phyloseq::psmelt(GP.candida)
library(scales)
ggplot(p, aes(x = Remission.State, y = Abundance)) +scale_y_continuous(trans=pseudo_log_trans(base=10))+geom_boxplot(outlier.shape  = NA, width=0.25, color="blue1") +geom_jitter(color="blue1", height = 0, width = 0.2, size=4) +labs(x = "", y = "Abundance\n") +facet_wrap(~ OTU, scales = "free")+theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()


#HEAT TREES
Install.packages(“metacoder”)
Library(metacoder) 
Library(tidverse)


#set up for metacoder

ps.filter<- phyloseq_filter_taxa_tot_fraction(ps, frac = 0.01)
ps.genus = tax_glom(ps.filter, taxrank="Genus", NArm=FALSE)
otu4<-otu_table(ps.genus)
class(otu4)
otu4<-as.matrix(otu4)
otu4<-as.data.frame(otu4)
otu4<-add_column(otu4, lineage = NA, .before = "928540")
tax4<-tax_table(ps.genus)
tax4<-as.matrix(tax4)
tax4<-as.data.frame(tax4)
tax4$taxonomy <- paste(tax4$Domain,tax4$Phylum, tax4$Class, tax4$Order, tax4$Family, tax4$Genus, sep=";")
otu4$lineage<-tax4$taxonomy
otu4<-add_column(otu4, OTU_ID = NA, .before = "lineage")
otu4$OTU_ID<-rownames(otu4)
rownames(otu4)<-1:nrow(otu4)
meta4<-sample_data(ps.genus)
meta4<-as.matrix(meta4)
meta4<-as.data.frame(meta4)
meta4<-add_column(meta4, sample.id = NA, .before = "DEIDENTIFIED_MASTER_PATIENT_ID")
meta4$sample.id<-rownames(meta4)



#compare mean abundances in remission state

obj <- parse_tax_data(otu4,
+                       class_cols = "lineage", # the column that contains taxonomic information
+                       class_sep = ";", # The character used to separate taxa in the classification
+                       class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
+                       class_key = c(tax_rank = "info", # A key describing each regex capture group
+                                     tax_name = "taxon_name"))


obj$data$tax_data <- zero_low_counts(obj, dataset = "tax_data", min_count = 5)
no_reads <- rowSums(obj$data$tax_data[, rownames(meta4)]) == 0
obj <- filter_obs(obj, target = "tax_data", ! no_reads, drop_taxa = TRUE)
meta4$sample.id<-as.character(meta4$sample.id)
obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data",
+                                        cols = meta4$sample.id)

obj$data$otu_props <- calc_obs_props(obj, "tax_abund", other_cols = TRUE)

testmeta<-meta4[!is.na(meta4$Remission.State),]
testmeta$Remission.State<-as.character(testmeta$Remission.State)
obj$data$type_abund <- calc_group_mean(obj, "tax_abund",
+                                        cols = testmeta$sample.id,
+                                        groups = testmeta$Remission.State)
print(obj$data$type_abund)


set.seed(1)
heat_tree(obj, 
          node_label = obj$taxon_names(),
          node_size = obj$n_obs(),
          node_color = obj$data$type_abund $"allgroups", 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations


# compare mean differences heat tree

testmeta<-meta4[!is.na(meta4$Remission.State),]

trial = function(abund_1, abund_2) {log_ratio <- log(mean(abund_1) / mean(abund_2))
+ if(is.nan(log_ratio)) {log_ratio <- 0}
+ list(log_mean_ratio = log_ratio, median_diff = median(abund_1) - median(abund_2), mean_diff = mean(abund_1) - mean(abund_2), wilcox_p_value = wilcox.test(abund_1, abund_2)$p.value)}

obj$data$diff_table <- compare_groups(obj,dataset = "tax_abund",cols = testmeta$sample.id, func = trial, groups = testmeta$Remission.State)

print(obj$data$diff_table)
pcr_success_color_scale = c(viridis::viridis(50)[25:45])

OR
makepalette<-colorRampPalette(brewer.pal(9,"Blues"))(100)
pcr_success_color_scale = c(makepalette[30:70])

heat_tree(obj, node_label = taxon_names, node_size = n_obs, node_color = -log_mean_ratio, node_color_interval = c(5, -5), node_color_range = pcr_success_color_scale, node_size_axis_label = "ASV count", node_color_axis_label = "Log Mean Ratio of Relative Abundance")


#Calculate proportions of specific taxa reads compared to total 

ps.subset = subset_samples(ps, Remission.State != "NA")
ps.phylum = tax_glom(ps.subset, taxrank="Phylum", NArm=FALSE)
get_taxa_unique(ps.phylum, "Phylum")
sum(sample_sums(ps.phylum))
GP.chl <- subset_taxa(ps.phylum, Phylum=="p__unidentified")
sum(sample_sums(GP.chl) / sum(sample_sums(ps.phylum))

