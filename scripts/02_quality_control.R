# 02_qc.R

# Quality control and exploration of TCGA glioma RNA-seq data

# Input: "data/processed/tcga_glioma_se.rds"
# Output: 

library(SummarizedExperiment)

# 1. Loading and confirming the presence of data in the processed file

se = readRDS("data/processed/tcga_glioma_se.rds")
dim(se)                               # Should be equal to 60660 888

# 2. Understanding the samples

names(colData(se))
idh = se$paper_IDH.status
grade = se$paper_Grade
table(idh, grade, useNA = "ifany")
se = se[,!duplicated(se$patient)]     # One sample per patient (keeping the first), samples originating from the same patient are not independent
dim(se)                               # Should now be equal to 60660 800

# 3. Defining the groups of interest

se = se[,!is.na(se$paper_IDH.status)] # Removing samples without IDH status, these cannot be assigned to a group
dim(se)                               # Should now be equal to 60660 780
idh = se$paper_IDH.status             # Re-extracting the idh status after filtering
se$IDH = factor(idh, levels = c("WT","Mutant"))   # Creating the IDH factor, with "WT" as the reference, used downstream for analysis
levels(se$IDH)                        # Verifying the levels: "WT" and "Mutant" 
table(se$IDH, useNA = "ifany")        # Number of cases for each level

# 4. Extracting raw counts and computing library sizes

assayNames(se)                        # Show the matrices present in the se object
raw_counts = assay(se, "unstranded")  # Extracting the raw counts
class(raw_counts)                     # Verifying that the object "raw_counts" is a matrix
dim(raw_counts)                       # Verifying the dimensions of the "raw_counts", should be equal to 60660 780
lib_size = colSums(raw_counts) / 1e6             # Sum of all columns + Adjustment to millions of reads
summary(lib_size)                     # Summary of the library sizes
hist(lib_size)                        # Histogram (x = library size, in millions; y = how many samples fall into each range)
 
# 5. Filtering genes with low counts

n_samples = rowSums(raw_counts >= 10) # For each gene: how many samples >= 10
n = min(table(se$IDH))                # Determining the size of the smallest group
keep = n_samples >= n                 # Deciding which genes to keep
length(keep)                          # Confirmation that the number of genes remains the same
sum(keep)                             # Quantity of genes that are kept after filtering
se = se[keep,]                        # Keeping only the non-low count genes
nrow(se)                              # Verifying that the number of genes has changed, should be equal to 21923

# 6. DESeq 2 object vs vst

library(DESeq2)                       # Installation of the DESeq2 package
dds = DESeqDataSet(se, design = ~ IDH)           # DESeq2 object: counts + metadata from se, comparing IDH groups
