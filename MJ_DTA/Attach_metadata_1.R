# Load Packages -----------------------------------------------------------
require(dplyr)
require(seqinr)
require(ape)
require(phangorn)
require(ggtree)
require(ggplot2)
require(stringr)

latestdate <- "20220902"

# Load data ---------------------------------------------------------------
meta <- read.csv(paste0("./In/latest_curated_metadata_AFGPAKMWIMOZ_", latestdate, ".csv"))
seqs <- read.FASTA(paste0("./In/latest_curated_seq_AFGPAKMWIMOZ_", latestdate,".fasta"), type = "DNA")


# Match and add metadata to labels ----------------------------------------
## check all are present
table(names(seqs) %in% meta$seqname)
names(seqs)[which(!(names(seqs) %in% meta$seqname))] 
## drop? 
seqs <- seqs[-which(!(names(seqs) %in% meta$seqname))]

## New column with additional info
meta$seqlabs <- paste(meta$seqname, meta$sample_type, meta$adm0_name, meta$sample_date_num, meta$cluster, gsub("_", "-",meta$region),    sep = "_")

## order metadata to match sequence data
order_meta <- left_join(data.frame(names(seqs)), meta, by = c("names.seqs."="seqname"))

# Export new data ---------------------------------------------------------
write.csv(order_meta, paste0("./Meta_matched/longlabs_metadata_AFGPAKMWIMOZ_", latestdate,".csv"))

names(seqs) <- order_meta$seqlabs

write.FASTA(seqs, paste0("./Meta_matched/longlabs_seq_AFGPAKMWIMOZ_",latestdate,".fasta"))
seqs <- read.FASTA(paste0("./Meta_matched/longlabs_seq_AFGPAKMWIMOZ_",latestdate,"_Sabin.fasta")) ## add in Sabin to see these 2012s
# Check root-to-tip -------------------------------------------------------
require(ggtree)
require(stringr)
require(ggplot2)

## convert to phydat
seqs_phydat <- phyDat(seqs, type = "DNA")

## run a model test - v.slow
# mt <- modelTest(seqs_phydat)

# HKY + G is fine

dna_dist <- dist.ml(seqs_phydat, model="JC69")

NJ_tree <- NJ(dna_dist)


## Run IQtree from the path
fn = paste0("./Meta_matched/longlabs_seq_AFGPAKMWIMOZ_",latestdate,".fasta")
system( paste0( 'iqtree -nt AUTO -redo -m HKY -s ', fn ), intern=FALSE)


# load in the generated tree and root it
IQtr <- read.tree(paste0(fn, ".treefile"))
dates <- sapply( strsplit( IQtr$tip.label, '\\_' ), function(x){ as.numeric( tail(x,3)[1] )})

# root to give a consistent molecular clock
IQ_root <- rtt(IQtr, dates)

# devtools::install_github("lucymli/EpiGenR")

# Returns a data frame with time and divergence for each sequence
mc.df <- EpiGenR::root2tip.divergence(tr=IQ_root,tip.dates = dates)

# Linear regression
lm.res = lm(divergence~time,data=mc.df)

ggplot(mc.df)+
  geom_point(aes(x=time, y=divergence))+
  geom_abline(intercept = lm.res$coefficients[1], slope = lm.res$coefficients[2], color = "red", size = 1.5)+
  theme_classic()+
  annotate("text", label = paste0("No. subst/site/year = ",round(lm.res$coefficients[2],3)), x= 2014.5, y = 0.16)+
  annotate("text", label = paste0("R2 = ",round(summary(lm.res)$r.squared,2)), x= 2014.5, y = 0.156)+
  labs(x = "Year", y = "Root-to-tip distance")

ggsave(paste0("./rtt",latestdate,".png"), width = 5, height = 5)

tipdata <- do.call(rbind.data.frame, str_split(IQ_root$tip.label, "_"))
colnames(tipdata) <- c("seqname", "sample_type", "adm0_name", "sample_date_num", "cluster", "region")
tipdata$node <- as.numeric(rownames(tipdata))

ggtree(IQ_root)  %<+% tipdata +
  geom_tippoint(aes(fill = region), shape = 21, size =2, color = "black") +
  scale_fill_manual(values = getPalette, name = "Location")+
  guides(fill = guide_legend(override.aes = list(size = 5, shape = 22), ncol=2))+
  theme(legend.position = c(.25, .8),
        legend.box.background = element_rect(colour = "black"))

ggsave(paste0("./ML_root_location_",latestdate,".png"), width = 10, height = 7)
ggsave(paste0("./ML_root_location_",latestdate,".svg"), width = 10, height = 7)
