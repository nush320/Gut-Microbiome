# Gut Microbiome Meta-Analysis: Insights into African American Representation  

## Overview  
This repository hosts the workflows and scripts used in my gut microbiome meta-analysis project, centered on African American representation in microbiome studies.  

The project compiles and integrates publicly available gut microbiome datasets, applying standardized pipelines for quality control, taxonomic profiling, diversity analysis, and statistical testing.  

The overarching goal is to:  
- Highlight disparities in microbiome research participation.  
- Provide reproducible workflows for conducting large-scale microbiome meta-analyses.  
- Generate insights into microbial diversity and composition across ethnicities.  

---

## Project Workflow  

### 1. Data Collection  
- Systematic literature review.  
- Selection of datasets from public repositories, i.e., NCBI SRA.  
- Inclusion of fecal microbiome datasets with available ethnicity metadata.  

### 2. Data Processing  
- Raw FASTQ file download.
- Data cleaning and metadata arrangement.
- DADA2 Analysis: filter, denoise, and QC  
- Creation of ASV/OTU table, Taxonomy Table (alignment with SILVA 138.1), and Phyloseq Object   

### 3. Meta-Analysis  
- Diversity analyses (alpha & beta diversity).  
- Differential abundance testing.  
- Subsetting samples for **African American populations** and comparison across disease versus healthy cohorts.  
- Machine learning application  

---

## Repository Structure  
```bash
├── data/               # Metadata files and processed datasets (no raw FASTQ)
├── scripts/            # Data cleaning, preprocessing, and analysis scripts
│   ├── preprocessing/  # QC, trimming, filtering steps
│   ├── analysis/       # Diversity, statistical, and meta-analysis code
│   └── visualization/  # Plotting and figure generation
├── results/            # Outputs, figures, and statistical summaries
├── docs/               # Documentation, notes, and references
└── README.md           # Project description (this file)
