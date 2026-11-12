
### These are the inital steps that must be taken in order to download the TCGA counts

# 01_download_tcga.R
# Download TCGA-LGG and TCGA-GBM RNA-seq counts (STAR - Counts) from the GDC

library(TCGAbiolinks)
library(SummarizedExperiment)

# 1. Query the GDC

query <- GDCquery(project       = c("TCGA-GBM", "TCGA-LGG"),
                  data.category = "Transcriptome Profiling",
                  data.type     = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

# 2. Checking what will be downloaded

res_query <- getResults(query)
table(res_query$project, res_query$sample_type)

# 2a. Saving the exact list of files used

dir.create("metadata", showWarnings = FALSE)
write.csv(res_query[, c("id", "file_name", "md5sum", "project", "sample_type", "cases")],
          "metadata/gdc_manifest.csv", row.names = FALSE)

# 3. Download (several GB)

# GDCdownload() failed to move files into the correct folder (project folder) on Windows,
# so files were downloaded directly from the GDC API into a short path
# outside the project. Change gdc_dir on another computer (e.g. "~/gdc").

gdc_dir <- "C:/gdc"

for (i in seq_len(nrow(res_query))) {
  dest <- file.path(gdc_dir, res_query$project[i], "Transcriptome_Profiling",
                    "Gene_Expression_Quantification", res_query$id[i])
  file <- file.path(dest, res_query$file_name[i])
  if (!file.exists(file)) {                         # skips files already downloaded
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    download.file(paste0("https://api.gdc.cancer.gov/data/", res_query$id[i]),
                  file, mode = "wb", quiet = TRUE)
  }
}

# 3a. Checking that all files are complete (md5 checksums)

files <- file.path(gdc_dir, res_query$project, "Transcriptome_Profiling",
                   "Gene_Expression_Quantification", res_query$id, res_query$file_name)
all(file.exists(files))                              # should be TRUE
all(tools::md5sum(files) == res_query$md5sum)        # should be TRUE

# 4. Build the SummarizedExperiment object

gdc_dir <- "C:/gdc"
se <- GDCprepare(query, directory = gdc_dir)
table(se$project_id, se$sample_type)

# GDCprepare() leaves intermediate files in the project folder -> Remove them

unlink(c("df.rds", "results.rds"))

# 5. Keeps the primary tumours only (avoiding bias from post-treatment tumours)

se <- se[, se$sample_type == "Primary Tumor"]

# 6. Look for the IDH status column

grep("idh|grade", names(colData(se)), value = TRUE, ignore.case = TRUE)

# 7. Save

saveRDS(se, "data/processed/tcga_glioma_se.rds")
sessionInfo()
