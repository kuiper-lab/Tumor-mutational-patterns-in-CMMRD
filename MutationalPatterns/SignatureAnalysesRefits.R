##SIGNATURE ANALYSIS AND REFITS

library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

#_______________________________________________________________________________________________________________________________________________________________________________________________________________________#
##Load Rdata SNV signatures extracted, then proceed below
nmf_res = nmf_res_sbs_mutation_matrix_combined
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", "Signature H", "Signature I", "Signature J", "Signature K", "Signature L")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", "Signature H", "Signature I", "Signature J", "Signature K", "Signature L")
##Load COSMIC signatures
signatures = get_known_signatures()
##Adapt signature names to COSMIC signatures based on cosine similarity
nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
##For signatures with a cosine similarity to COSMIC signatures lower than 0.85
signatures_snv_denovo = nmf_res$signatures[,c("<SIGNATURES>")]
for (col in 1:ncol(signatures_indel_denovo)){
  signatures_denovo_scaled[, col] <- signatures_snv_denovo[, col] / sum(signatures_snv_denovo[, col])
}
strict_refit <- fit_to_signatures_strict(signatures_denovo_scaled[], signatures, max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res
##If the extracted signature can be explained by 3 COSMIC signatures or less, check if the cosine similarity between the original and recontructed profile is >0.8
cos_sim(signatures_snv_denovo[], fit_res_strict$reconstructed[])

##ONCE THE DEFINITIVE SET OF COSMIC AND DE NOVO SIGNATURES IS DETERMINED

##Load the CMMRD mutational matrix
vcf_files <- list.files(path = "<PATH_TO_VCF_FILES", pattern = ".vcf", full.names = TRUE)
sample_names = c("SAMPLE_NAMES_VCF_FILES")
CMMRD_vcfs <- read_vcfs_as_granges(vcf_files,
                                   sample_names,
                                   ref_genome,
                                   type="snv")
mut_mat <- mut_matrix(vcf_list = CMMRD_vcfs, ref_genome = ref_genome)
##Perform strict refit
signatures_snv_12_denovo = signatures[,c( "SBS1", "SBS15", "SBS14", "SBS7a", "SBS20","SBS26", "SBS5", "SBS11", "SBS18", "SBS17b", "SBS2", "SBS13", "SBS31", "SBS8")]
strict_refit <- fit_to_signatures_strict(mut_mat, signatures_snv_12_denovo, max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res
##Calculate cosine similarity between orginal and reconstructed profiles
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat[], fit_res_strict$reconstructed)

#_______________________________________________________________________________________________________________________________________________________________________________________________________________________#
##Load Rdata INDEL signatures extracted, then proceed below
nmf_res = nmf_res_indel_mutation_matrix_combined
colnames(nmf_res$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
##Load COSMIC signatures
signatures_indel = get_known_signatures(muttype = "indel")
##Adapt signature names to COSMIC signatures based on cosine similarity
nmf_res <- rename_nmf_signatures(nmf_res, signatures_indel, cutoff = 0.85)
##For signatures with a cosine similarity to COSMIC signatures lower than 0.85
signatures_indel_denovo = nmf_res$signatures[,c("<SIGNATURES>")]
for (col in 1:ncol(signatures_indel_denovo)){
  signatures_denovo_scaled[, col] <- signatures_indel_denovo[, col] / sum(signatures_indel_denovo[, col])
}
strict_refit <- fit_to_signatures_strict(signatures_denovo_scaled[], signatures_indel, max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res
##If the extracted signature can be explained by 3 COSMIC signatures or less, check if the cosine similarity between the original and recontructed profile is >0.8
cos_sim(signatures_indel_denovo[], fit_res_strict$reconstructed[])

##ONCE THE DEFINITIVE SET OF COSMIC AND DE NOVO SIGNATURES IS DETERMINED

##Load the CMMRD mutational matrix
vcf_files <- list.files(path = "<PATH_TO_VCF_FILES", pattern = ".vcf", full.names = TRUE)
sample_names = c("SAMPLE_NAMES_VCF_FILES")
CMMRD_vcfs <- read_vcfs_as_granges(vcf_files,
                                   sample_names,
                                   ref_genome,
                                   type="indel")
CMMRD_indel_vcfs <- get_indel_context(CMMRD_vcfs, ref_genome)
CMMRD_indel_counts <- count_indel_contexts(CMMRD_indel_vcfs)
##Perform strict refit
signatures_indel_6_denovo = signatures[,c("ID1", "ID2", "ID12", "IDA", "IDB", "IDC")]
strict_refit <- fit_to_signatures_strict(CMMRD_indel_counts, signatures_indel_6_denovo, max_delta = 0.033)
fit_res_strict <- strict_refit$fit_res
##Calculate cosine similarity between orginal and reconstructed profiles
cos_sim_samples_signatures <- cos_sim_matrix(CMMRD_indel_counts[], fit_res_strict$reconstructed)

#_________________________________________________________________________________________________________________#
##Additional analyses SBS11
##Load CMMRD samples
vcf_files <- list.files(path = "<PATH_TO_VCF_FILES", pattern = ".vcf", full.names = TRUE)
sample_names = c("SAMPLE_NAMES_VCF_FILES")
CMMRD_vcfs <- read_vcfs_as_granges(vcf_files,
                                   sample_names,
                                   ref_genome,
                                   type="snv")
mut_mat <- mut_matrix(vcf_list = CMMRD_vcfs, ref_genome = ref_genome)
##Select tumors with SBS11
mut_mat_SBS11 = mut_mat[,c("<tumors with SBS11>")]
##Load COSMIC signatures
signatures = get_known_signatures()
signatures_extraction = signatures[,c( "SBS1", "SBS15", "SBS14", "SBS7a", "SBS20","SBS26", "SBS5", "SBS11", "SBS18", "SBS17b", "SBS2", "SBS13", "SBS31", "SBS8")]
signatures_extraction_MMR = signatures[,c( "SBS1", "SBS15", "SBS14", "SBS7a", "SBS20","SBS26", "SBS5", "SBS11", "SBS18", "SBS17b", "SBS2", "SBS13", "SBS31", "SBS8","SBS6", "SBS21","SBS44")]
##Perform bootstrapped refit
contri_boots_tumSig1 <- fit_to_signatures_bootstrapped(mut_mat_SBS11[],
                                                       signatures,
                                                       n_boots = 200,
                                                       max_delta = 0.033,
                                                       method = "strict")
contri_boots_tumSig2 <- fit_to_signatures_bootstrapped(mut_mat_SBS11[],
                                                       signatures_extraction,
                                                       n_boots = 200,
                                                       max_delta = 0.033,
                                                       method = "strict")
contri_boots_tumSig3 <- fit_to_signatures_bootstrapped(mut_mat_SBS11[],
                                                       signatures_extraction_MMR,
                                                       n_boots = 200,
                                                       max_delta = 0.033,
                                                       method = "strict")
##Strict refit between tumors with SBS11 and extracted signatures with and without SBS11
signatures_extraction_noSBS11 = signatures[,c( "SBS1", "SBS15", "SBS14", "SBS7a", "SBS20","SBS26", "SBS5", "SBS11", "SBS18", "SBS17b", "SBS2", "SBS13", "SBS31", "SBS8")]
strict_refit1 <- fit_to_signatures_strict(mut_mat_SBS11[], signatures_extraction, max_delta = 0.033)
fit_res_strict1 <- strict_refit1$fit_res
strict_refit2 <- fit_to_signatures_strict(mut_mat_SBS11[], signatures_extraction_noSBS11, max_delta = 0.033)
fit_res_strict2 <- strict_refit2$fit_res

cos_sim_samples_signatures <- cos_sim_matrix(mut_mat_SBS11[], fit_res_strict1$reconstructed)
cos_sim_samples_signatures <- cos_sim_matrix(mut_mat_SBS11[], fit_res_strict2$reconstructed)





