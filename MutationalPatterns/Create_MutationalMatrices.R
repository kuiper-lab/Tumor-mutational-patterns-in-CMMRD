##CREATE MATRICES

library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

##SNV 
##Load VCF files
vcf_files <- list.files(path = "<PATH_TO_VCF_FILES", pattern = ".vcf", full.names = TRUE)
sample_names = c("SAMPLE_NAMES_VCF_FILES")
CMMRD_vcfs <- read_vcfs_as_granges(vcf_files,
                                   sample_names,
                                   ref_genome,
                                   type="snv")
##Create mutation matrix
mut_mat <- mut_matrix(vcf_list = CMMRD_vcfs, ref_genome = ref_genome)
mut_mat_control = read.delim("<CONTROL_COHORT1", row.names = 1)
mut_mat_all <- as.matrix(cbind(mut_mat, mut_mat_control))
save(mut_mat_all, file = "<FILE_PATH>/<FILE>.RData")

##INDEL 
##Load VCF files
vcf_files <- list.files(path = "<PATH_TO_VCF_FILES", pattern = ".vcf", full.names = TRUE)
sample_names = c("SAMPLE_NAMES_VCF_FILES")
CMMRD_vcfs <- read_vcfs_as_granges(vcf_files,
                                   sample_names,
                                   ref_genome,
                                   type="indel")
##Create mutation matrix
CMMRD_indel_vcfs <- get_indel_context(CMMRD_vcfs, ref_genome)
CMMRD_indel_counts <- count_indel_contexts(CMMRD_indel_vcfs)
mut_mat_control = read.delim("<CONTROL_COHORT1", row.names = 1)
mut_mat_all <- as.matrix(cbind(CMMRD_indel_counts, mut_mat_control))
save(mut_mat_all, file = "<FILE_PATH>/<FILE>.RData")





