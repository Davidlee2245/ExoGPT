# TCGA Data for Pancreatic cancer (PAAD)

## How to Download Full TCGA Data

1. **Using GDC Data Portal**:
   - Visit: https://portal.gdc.cancer.gov/
   - Select project: PAAD
   - Download mRNA expression data (HTSeq-FPKM or HTSeq-Counts)

2. **Using cBioPortal**:
   - Visit: https://www.cbioportal.org/
   - Select study: TCGA PAAD
   - Download mRNA expression data

3. **Using R/Bioconductor**:
   ```r
   library(TCGAbiolinks)
   query <- GDCquery(
     project = "TCGA-PAAD",
     data.category = "Transcriptome Profiling",
     data.type = "Gene Expression Quantification",
     workflow.type = "HTSeq - FPKM"
   )
   GDCdownload(query)
   data <- GDCprepare(query)
   ```

## Expected CSV Format

The CSV file should have columns:
- gene_symbol: Gene symbol (e.g., "TP53")
- ensembl_id: Ensembl gene ID (optional)
- sample_type: "Primary Tumor" or "Solid Tissue Normal"
- log2_fpkm: Log2-transformed FPKM values
- project_code: TCGA project code (e.g., "PAAD")

## Current Status

This is a placeholder file. Replace it with actual TCGA data for full functionality.
