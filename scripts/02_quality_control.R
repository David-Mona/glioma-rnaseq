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
