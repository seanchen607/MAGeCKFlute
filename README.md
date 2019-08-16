# MAGeCKFlute

## 1. sgRNA count by MAGeCK v2.0.5
    mageck count --trim-5 1 -l CRISPR_sgRNA_library.csv -n Cal27_NSD1KO --sample-label Cal27_day0,Cal27_day14,Cal27_NSD1KO_day0,Cal27_NSD1KO_day14  --fastq Cal27_day0.fastq Cal27_day14.fastq Cal27_NSD1KO_day0.fastq Cal27_NSD1KO_day14.fastq
----------------------------------------
## 2. Enrichment analysis using MAGeCK mle module
    mageck mle -k Cal27_NSD1KO.count_normalized.txt -d designmat.txt -n Cal27_NSD1KO_mle_normalized  
----------------------------------------
## 3. Integrative analysis pipeline for pooled CRISPR functional genetic screens - MAGeCKFlute
    ################### Run following codes in R #################
    # source("http://www.bioconductor.org/biocLite.R")
    # biocLite("MAGeCKFlute")
    library(MAGeCKFlute)
    ## Load gene summary data in MAGeCK MLE results
    mle.gene_summary = read.delim("Cal27_NSD1KO_mle_normalized.gene_summary.txt", header=T, sep="\t")
    ## Run the downstream analysis pipeline for MAGeCK MLE
    FluteMLE(mle.gene_summary, ctrlname="Cal27_day14", treatname="Cal27_NSD1KO_day14", prefix="MLE", organism="hsa")
    ###############################################################
----------------------------------------
## 4. Quality control
    ################### Run following codes in R #################
    library(MAGeCKFlute)
    ## Load gene summary data in MAGeCK MLE results
    countsummary = read.delim("Cal27_NSD1KO.countsummary.txt", header=T, sep="\t")

    MapRatesView(countsummary)

    IdentBarView(countsummary, x = "Label", y = "GiniIndex", ylab = "Gini index", main = "Evenness of sgRNA reads")

    countsummary$Missed = log10(countsummary$Zerocounts)
    IdentBarView(countsummary, x = "Label", y = "Missed", fill = "#394E80", ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")
    ###############################################################
----------------------------------------
