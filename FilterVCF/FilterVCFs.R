##Load library
library(mutationalProcesses)



#
##
###Functions
##
#

## Load VCF file
loadVcfFile <- function (vcfFilePath, bsgGenomeName, chromosomeSet) 
{
  require(VariantAnnotation)
  require(logging)
  CHROMOSOMES_TO_LOAD <- loadPredefinedChromosomeSet(bsgGenomeName, 
                                                     chromosomeSet)
  referenceBsGenome <- loadReferenceBsGenome(bsgGenomeName)
  ref_style <- seqlevelsStyle(referenceBsGenome)
  ref_organism <- GenomeInfoDb::organism(referenceBsGenome)
  vcf <- VariantAnnotation::readVcf(vcfFilePath, bsgGenomeName)
  loginfo(paste0("loaded ", dim(vcf)[[1]], " variants from file: ", 
                 vcfFilePath))
  loginfo(paste0("changing style from: ", paste(seqlevelsStyle(vcf), 
                                                collapse = ", "), " to: ", ref_style, " for vcf-file: ", 
                 vcfFilePath))
  seqlevelsStyle(vcf) <- ref_style
  loginfo(paste0("keep seq levels ", paste(CHROMOSOMES_TO_LOAD, 
                                           collapse = ", "), " for vcf-file: ", vcfFilePath))
  chromosomes_present_and_to_keep <- CHROMOSOMES_TO_LOAD[CHROMOSOMES_TO_LOAD %in% 
                                                           levels(unique(seqnames(rowRanges(vcf))))]
  vcf <- keepSeqlevels(x = vcf, value = chromosomes_present_and_to_keep, 
                       pruning.mode = "coarse")
  return(vcf)
}


## Load chromosome set
loadPredefinedChromosomeSet <- function (referenceGenomeBuild = c("BSgenome.Hsapiens.UCSC.hg19", 
                                                                  "BSgenome.Hsapiens.NCBI.GRCh38", "BSgenome.Hsapiens.UCSC.hg38"), 
                                         chromosomeSet = c("chromosomes", "autosomes", "autosomes_and_x")) 
{
  require(logging)
  referenceGenomeBuild <- match.arg(referenceGenomeBuild)
  chromosomeSet <- match.arg(chromosomeSet)
  HG_ALL_CHROMOSOMES <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                          "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
                          "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
  HG_AUTOSOMES <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                    "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                    "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
                    "chr19", "chr20", "chr21", "chr22")
  HG_AUTOSOMES_AND_X <- c("chr1", "chr2", "chr3", "chr4", "chr5", 
                          "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
                          "chr19", "chr20", "chr21", "chr22", "chrX")
  B_ALL_CHROMOSOMES <- c("1", "2", "3", "4", "5", "6", "7", 
                         "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
                         "18", "19", "20", "21", "22", "X", "Y")
  B_AUTOSOMES <- c("1", "2", "3", "4", "5", "6", "7", "8", 
                   "9", "10", "11", "12", "13", "14", "15", "16", "17", 
                   "18", "19", "20", "21", "22")
  B_AUTOSOMES_AND_X <- c("1", "2", "3", "4", "5", "6", "7", 
                         "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
                         "18", "19", "20", "21", "22", "X")
  predefinedChromosomeSet <- NULL
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.NCBI.GRCh38") && 
      (chromosomeSet == "chromosomes")) {
    predefinedChromosomeSet <- B_ALL_CHROMOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.NCBI.GRCh38") && 
      (chromosomeSet == "autosomes")) {
    predefinedChromosomeSet <- B_AUTOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.NCBI.GRCh38") && 
      (chromosomeSet == "autosomes_and_x")) {
    predefinedChromosomeSet <- B_AUTOSOMES_AND_X
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg19") && 
      (chromosomeSet == "chromosomes")) {
    predefinedChromosomeSet <- HG_ALL_CHROMOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg19") && 
      (chromosomeSet == "autosomes")) {
    predefinedChromosomeSet <- HG_AUTOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg19") && 
      (chromosomeSet == "autosomes_and_x")) {
    predefinedChromosomeSet <- HG_AUTOSOMES_AND_X
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg38") && 
      (chromosomeSet == "chromosomes")) {
    predefinedChromosomeSet <- HG_ALL_CHROMOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg38") && 
      (chromosomeSet == "autosomes")) {
    predefinedChromosomeSet <- HG_AUTOSOMES
  }
  if ((referenceGenomeBuild == "BSgenome.Hsapiens.UCSC.hg38") && 
      (chromosomeSet == "autosomes_and_x")) {
    predefinedChromosomeSet <- HG_AUTOSOMES_AND_X
  }
  loginfo(paste0("loaded: ", chromosomeSet, " for genome build: ", 
                 referenceGenomeBuild))
  return(predefinedChromosomeSet)
}


## Load reference genome
loadReferenceBsGenome <- function (bsGenomeReferenceGenomeName) 
{
  require(BSgenome)
  require(logging)
  library(bsGenomeReferenceGenomeName, character.only = TRUE)
  ref_genome <- base::get(bsGenomeReferenceGenomeName)
  ref_organism <- GenomeInfoDb::organism(ref_genome)
  ref_style <- seqlevelsStyle(ref_genome)
  genome_name <- genome(ref_genome)[[1]]
  loginfo(paste0("Loading genome: ", bsGenomeReferenceGenomeName, 
                 " (organism: ", ref_organism, " with reference style: ", 
                 ref_style, ")"))
  return(ref_genome)
}



#
##
###Run code
##
#

##FILTERING SNVs/MNPs

##Command needed later in the filtering process
'%nin%' = Negate('%in%')
##Load VCF
<INDIVIDUAL_ID>_Mutect2 <- loadVcfFile("<SNV_VCF_FILE>",
                               "BSgenome.Hsapiens.UCSC.hg38",
                               "chromosomes")
##Load associated TXT file (filtering is performed on the TXT file)
<INDIVIDUAL_ID>_Mutect2_dfNP <- read.delim("<SNV_TXT_FILE>", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
##Correct applied phasing
x=<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.GT
y=<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.REFCOUNT
z=<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.ALTCOUNT
tot_lines=nrow(<INDIVIDUAL_ID>_Mutect2_dfNP)

for (i in 1:tot_lines){
  ifelse(x[i]=="1|0", {y[i]=y[i]+z[i];z[i]=y[i]-z[i];y[i]=y[i]-z[i];x[i]=paste("0|1i")}, next)
}
<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.GT = x
<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.REFCOUNT = y
<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.ALTCOUNT = z
<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.VAF= <INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.ALTCOUNT/(<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.ALTCOUNT+<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.REFCOUNT)
##PASS filter data
<INDIVIDUAL_ID>_Mutect2_df_PASS = subset(<INDIVIDUAL_ID>_Mutect2_dfNP, <INDIVIDUAL_ID>_Mutect2_dfNP$FILTER=="PASS")
##Create rownames
KR_rownames = paste(<INDIVIDUAL_ID>_Mutect2_df_PASS$CHROM, ":", <INDIVIDUAL_ID>_Mutect2_df_PASS$POS, "_", <INDIVIDUAL_ID>_Mutect2_df_PASS$REF, "/", <INDIVIDUAL_ID>_Mutect2_df_PASS$ALT, sep = "")
row.names(<INDIVIDUAL_ID>_Mutect2_df_PASS) = KR_rownames
##Filter to keep only variants with 5 or more alternative reads
<INDIVIDUAL_ID>_Mutect2_df_ALTR_5 <- subset(<INDIVIDUAL_ID>_Mutect2_df_PASS, <INDIVIDUAL_ID>_Mutect2_df_PASS$<TUMOR_SAMPLE_ID>.ALTCOUNT>=5)
##Filter to keep only variants with an alternative allele frequence of 0.1 or higher
<INDIVIDUAL_ID>_Mutect2_df_AF0.1 <- subset(<INDIVIDUAL_ID>_Mutect2_df_ALTR_5, <INDIVIDUAL_ID>_Mutect2_df_ALTR_5$<TUMOR_SAMPLE_ID>.VAF>=0.1)
##Filter to keep only variants with 0 reads in the normal sample
<INDIVIDUAL_ID>_Mutect2_df_GL0 <- subset(<INDIVIDUAL_ID>_Mutect2_df_AF0.1, <INDIVIDUAL_ID>_Mutect2_df_AF0.1$<NORMAL_SAMPLE_ID>.VAF==0)
##Seperate multinucleotide variants (mnp) from single nucleotide variants (snv)
<INDIVIDUAL_ID>_Mutect2_df_mnp = data.frame(matrix(ncol = <number_of_columns_in_df, nrow = 0))
colnames(<INDIVIDUAL_ID>_Mutect2_df_mnp) = colnames(<INDIVIDUAL_ID>_Mutect2_df_GL0)
variant = NULL
position = NULL
next_position = NULL
chromosome = NULL
next_chromosome = NULL
for (variant in 1:(nrow(<INDIVIDUAL_ID>_Mutect2_df_GL0)-1)){
  position <- <INDIVIDUAL_ID>_Mutect2_df_GL0$POS[variant]
  next_position <- <INDIVIDUAL_ID>_Mutect2_df_GL0$POS[variant+1]
  chromosome <- <INDIVIDUAL_ID>_Mutect2_df_GL0$CHROM[variant]
  next_chromosome <- <INDIVIDUAL_ID>_Mutect2_df_GL0$CHROM[variant+1]
  if (chromosome == next_chromosome & (next_position > (position -2)) & (next_position < (position + 2))){
    <INDIVIDUAL_ID>_Mutect2_df_mnp <- rbind(<INDIVIDUAL_ID>_Mutect2_df_mnp, <INDIVIDUAL_ID>_Mutect2_df_GL0[variant:(variant+1),])
  }
}

<INDIVIDUAL_ID>_Mutect2_df_snv = <INDIVIDUAL_ID>_Mutect2_df_GL0[rownames(<INDIVIDUAL_ID>_Mutect2_df_GL0) %nin% rownames(<INDIVIDUAL_ID>_Mutect2_df_mnp),]

##Convert GnomAD and GoNL (population frequency) values to numeric
<INDIVIDUAL_ID>_Mutect2_df_snv$gnomADg_AF <- as.numeric(<INDIVIDUAL_ID>_Mutect2_df_snv$gnomADg_AF)
<INDIVIDUAL_ID>_Mutect2_df_snv$goNL_AF <- as.numeric(<INDIVIDUAL_ID>_Mutect2_df_snv$goNL_AF)

##Filter to keep only variants without GnomAD and GoNL frequencies, or frequencies of 0.01 or lower
<INDIVIDUAL_ID>_Mutect2_df_GnomAD_NA <- subset(<INDIVIDUAL_ID>_Mutect2_df_snv, is.na(<INDIVIDUAL_ID>_Mutect2_df_snv$gnomADg_AF))
<INDIVIDUAL_ID>_Mutect2_df_GnomAD_0.01 <- subset(<INDIVIDUAL_ID>_Mutect2_df_snv, <INDIVIDUAL_ID>_Mutect2_df_snv$gnomADg_AF<= 0.01)
<INDIVIDUAL_ID>_Mutect2_df_GnomAD <- rbind(<INDIVIDUAL_ID>_Mutect2_df_GnomAD_NA, <INDIVIDUAL_ID>_Mutect2_df_GnomAD_0.01)

<INDIVIDUAL_ID>_Mutect2_df_GoNL_NA <- subset(<INDIVIDUAL_ID>_Mutect2_df_GnomAD, is.na(<INDIVIDUAL_ID>_Mutect2_df_GnomAD$goNL_AF))
<INDIVIDUAL_ID>_Mutect2_df_GoNL_0.01 <- subset(<INDIVIDUAL_ID>_Mutect2_df_GnomAD, <INDIVIDUAL_ID>_Mutect2_df_GnomAD$goNL_AF<=0.01)
<INDIVIDUAL_ID>_Mutect2_df_GoNL <- rbind(<INDIVIDUAL_ID>_Mutect2_df_GoNL_NA, <INDIVIDUAL_ID>_Mutect2_df_GoNL_0.01)

##Output filtered VCF with snvs
<INDIVIDUAL_ID>_Mutect2_filtered <- <INDIVIDUAL_ID>_Mutect2[rownames(<INDIVIDUAL_ID>_Mutect2) %in% 
                                              rownames(<INDIVIDUAL_ID>_Mutect2_df_GoNL), ]
writeVcf(<INDIVIDUAL_ID>_Mutect2_filtered,"<OUTPUT_SNV_VCF_FILE>")
##Output filtered VCF with mnps
<INDIVIDUAL_ID>_Mutect2_mnp_filtered <- <INDIVIDUAL_ID>_Mutect2[rownames(<INDIVIDUAL_ID>_Mutect2) %in% 
                                                  rownames(<INDIVIDUAL_ID>_Mutect2_df_mnp), ]
writeVcf(<INDIVIDUAL_ID>_Mutect2_mnp_filtered,"<OUTPUT_MNP_VCF_FILE>")

##FILTERING INDELS

##Command needed later in the filtering process
'%nin%' = Negate('%in%')
##Load VCF
<INDIVIDUAL_ID>_Mutect2 <- loadVcfFile("<INDEL_VCF_FILE>",
                               "BSgenome.Hsapiens.UCSC.hg38",
                               "chromosomes")
##Load associated TXT file (filtering is performed on the TXT file)
<INDIVIDUAL_ID>_Mutect2_dfNP <- read.delim("<INDEL_TXT_FILE>", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
##Correct applied phasing
x=<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.GT
y=<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.REFCOUNT
z=<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.ALTCOUNT
tot_lines=nrow(<INDIVIDUAL_ID>_Mutect2_dfNP)

for (i in 1:tot_lines){
  ifelse(x[i]=="1|0", {y[i]=y[i]+z[i];z[i]=y[i]-z[i];y[i]=y[i]-z[i];x[i]=paste("0|1i")}, next)
}
<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.GT = x
<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.REFCOUNT = y
<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.ALTCOUNT = z
<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.VAF= <INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.ALTCOUNT/(<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.ALTCOUNT+<INDIVIDUAL_ID>_Mutect2_dfNP$<TUMOR_SAMPLE_ID>.REFCOUNT)
##PASS filter data
<INDIVIDUAL_ID>_Mutect2_df_PASS = subset(<INDIVIDUAL_ID>_Mutect2_dfNP, <INDIVIDUAL_ID>_Mutect2_dfNP$FILTER=="PASS")
##Create rownames
KR_rownames = paste(<INDIVIDUAL_ID>_Mutect2_df_PASS$CHROM, ":", <INDIVIDUAL_ID>_Mutect2_df_PASS$POS, "_", <INDIVIDUAL_ID>_Mutect2_df_PASS$REF, "/", <INDIVIDUAL_ID>_Mutect2_df_PASS$ALT, sep = "")
row.names(<INDIVIDUAL_ID>_Mutect2_df_PASS) = KR_rownames
##Filter to keep only variants with 5 or more alternative reads
<INDIVIDUAL_ID>_Mutect2_df_ALTR_5 <- subset(<INDIVIDUAL_ID>_Mutect2_df_PASS, <INDIVIDUAL_ID>_Mutect2_df_PASS$<TUMOR_SAMPLE_ID>.ALTCOUNT>=5)
##Filter to keep only variants with an alternative allele frequence of 0.1 or higher
<INDIVIDUAL_ID>_Mutect2_df_AF0.1 <- subset(<INDIVIDUAL_ID>_Mutect2_df_ALTR_5, <INDIVIDUAL_ID>_Mutect2_df_ALTR_5$<TUMOR_SAMPLE_ID>.VAF>=0.1)
##Filter to keep only variants with 0 reads in the normal sample
<INDIVIDUAL_ID>_Mutect2_df_GL0 <- subset(<INDIVIDUAL_ID>_Mutect2_df_AF0.1, <INDIVIDUAL_ID>_Mutect2_df_AF0.1$<NORMAL_SAMPLE_ID>.VAF==0)
##Convert GnomAD and GoNL (population frequency) values to numeric
<INDIVIDUAL_ID>_Mutect2_df_GL0$gnomADg_AF <- as.numeric(<INDIVIDUAL_ID>_Mutect2_df_GL0$gnomADg_AF)
<INDIVIDUAL_ID>_Mutect2_df_GL0$goNL_AF <- as.numeric(<INDIVIDUAL_ID>_Mutect2_df_GL0$goNL_AF)

##Filter to keep only variants without GnomAD and GoNL frequencies, or frequencies of 0.01 or lower
<INDIVIDUAL_ID>_Mutect2_df_GnomAD_NA <- subset(<INDIVIDUAL_ID>_Mutect2_df_GL0, is.na(<INDIVIDUAL_ID>_Mutect2_df_GL0$gnomADg_AF))
<INDIVIDUAL_ID>_Mutect2_df_GnomAD_0.01 <- subset(<INDIVIDUAL_ID>_Mutect2_df_GL0, <INDIVIDUAL_ID>_Mutect2_df_GL0$gnomADg_AF<= 0.01)
<INDIVIDUAL_ID>_Mutect2_df_GnomAD <- rbind(<INDIVIDUAL_ID>_Mutect2_df_GnomAD_NA, <INDIVIDUAL_ID>_Mutect2_df_GnomAD_0.01)

<INDIVIDUAL_ID>_Mutect2_df_GoNL_NA <- subset(<INDIVIDUAL_ID>_Mutect2_df_GnomAD, is.na(<INDIVIDUAL_ID>_Mutect2_df_GnomAD$goNL_AF))
<INDIVIDUAL_ID>_Mutect2_df_GoNL_0.01 <- subset(<INDIVIDUAL_ID>_Mutect2_df_GnomAD, <INDIVIDUAL_ID>_Mutect2_df_GnomAD$goNL_AF<=0.01)
<INDIVIDUAL_ID>_Mutect2_df_GoNL <- rbind(<INDIVIDUAL_ID>_Mutect2_df_GoNL_NA, <INDIVIDUAL_ID>_Mutect2_df_GoNL_0.01)

##Output filtered VCF with snvs
<INDIVIDUAL_ID>_Mutect2_filtered <- <INDIVIDUAL_ID>_Mutect2[rownames(<INDIVIDUAL_ID>_Mutect2) %in% 
                                              rownames(<INDIVIDUAL_ID>_Mutect2_df_GoNL), ]
writeVcf(<INDIVIDUAL_ID>_Mutect2_filtered,"<OUTPUT_INDEL_VCF_FILE>")

