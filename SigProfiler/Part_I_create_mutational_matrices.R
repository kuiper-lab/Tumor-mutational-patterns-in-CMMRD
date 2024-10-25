#Create mutational matrices
library(dplyr)
library(stringr)

#extract SBS and indels separately for each control file
files <- list.files("location/of/input/files", pattern="data_mutations_extended_control_")


ref_SNV <-c ("A", "T", "C", "G") 

for (file in files) {

  setwd("location/of/input/files")
 
  dataset <- read.delim(file, sep="\t", header = T)
  n <- str_replace(file, "data_mutations_extended_control_", "")
  name <- str_replace(n, ".txt", "")
  print(name)
  
  #SNVs
  SNV <- dataset %>%
    filter(Variant_Type == "SNP") %>%
    filter(Reference_Allele %in% ref_SNV) %>%
    select(Tumor_Sample_Barcode, Variant_Type, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2)
  
  
  SNV$Project <- name
  SNV$ID <- "1"
  SNV$Genome <- "GRCh37"
  SNV$Type <- "SOMATIC"
  
  SNV <- SNV[c(8,1,9,10,2,3,4,5,6,7,11)]
  
    
  SNV$Variant_Type <- str_replace(SNV$Variant_Type, "SNP", "SNV")
  
  colnames(SNV) <- c("Project",	"Sample",	"ID", "Genome",	"mut_type", "Chrom", "pos_start", "pos_end", "ref",	"alt",	"Type")
  
  
  #INDELS
  Indel <- dataset %>%
    filter(Variant_Type == "INS" | Variant_Type == "DEL") %>%
    select(Tumor_Sample_Barcode, Variant_Type, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2)
    
  
  Indel$Project <- name
  Indel$ID <- "1"
  Indel$Genome <- "GRCh37"
  Indel$Type <- "SOMATIC"
  
  Indel <- Indel[c(8,1,9,10,2,3,4,5,6,7,11)]
  
  
  
  Indel$Variant_Type <- str_replace(Indel$Variant_Type, "DEL", "INDEL")
  Indel$Variant_Type <- str_replace(Indel$Variant_Type, "INS", "INDEL")
  
  
  
  colnames(Indel) <- c("Project",	"Sample",	"ID", "Genome",	"mut_type", "Chrom", "pos_start", "pos_end", "ref",	"alt",	"Type")
  
  
  #save
  setwd("location/output/files")
  
  write.table(SNV, paste0(name, "_SNV_for_MMG_controls.txt"), row.names = F, col.names = T, sep="\t")
  

  setwd("location/output/files")
  
  write.table(Indel, paste0(name, "_Indel_for_MMG_controls.txt"), row.names = F, col.names = T, sep="\t")
}



################################
#extract SBS and indels separately for each test file
files <- list.files("location/of/input/files", pattern="data_mutations_extended_test_")

ref_SNV <-c ("A", "T", "C", "G") 

for (file in files) {
 
  setwd("location/of/input/files)
 
  
  dataset <- read.delim(file, sep="\t", header = T)
  n <- str_replace(file, "data_mutations_extended_test_", "")
  name <- str_replace(n, ".txt", "")
  print(name)
  
  #SNVs
  SNV <- dataset %>%
    filter(Variant_Type == "SNP") %>%
    filter(Reference %in% ref_SNV) %>%
    select(Tumor_ID, Variant_Type, Chromosome, Start_Position, End_Position, Reference, Variant)
  
  
  SNV$Project <- name
  SNV$ID <- "1"
  SNV$Genome <- "GRCh37"
  SNV$Type <- "SOMATIC"
  
  SNV <- SNV[c(8,1,9,10,2,3,4,5,6,7,11)]
  
    
  SNV$Variant_Type <- str_replace(SNV$Variant_Type, "SNP", "SNV")
  
  colnames(SNV) <- c("Project",	"Sample",	"ID", "Genome",	"mut_type", "Chrom", "pos_start", "pos_end", "ref",	"alt",	"Type")
  
  
  #INDELS
  Indel <- dataset %>%
    filter(Variant_Type == "INS" | Variant_Type == "DEL") %>%
    select(Tumor_ID, Variant_Type, Chromosome, Start_Position, End_Position, Reference, Variant)
    
  
  Indel$Project <- name
  Indel$ID <- "1"
  Indel$Genome <- "GRCh37"
  Indel$Type <- "SOMATIC"
  
  Indel <- Indel[c(8,1,9,10,2,3,4,5,6,7,11)]
  
  
  
  Indel$Variant_Type <- str_replace(Indel$Variant_Type, "DEL", "INDEL")
  Indel$Variant_Type <- str_replace(Indel$Variant_Type, "INS", "INDEL")
  
  
  
  colnames(Indel) <- c("Project",	"Sample",	"ID", "Genome",	"mut_type", "Chrom", "pos_start", "pos_end", "ref",	"alt",	"Type")
  
  
  #save
  setwd("location/output/files")
  
  write.table(SNV, paste0(name, "_SNV_for_MMG_test.txt"), row.names = F, col.names = T, sep="\t")
  

  setwd("location/output/files")
  
  write.table(Indel, paste0(name, "_Indel_for_MMG_test.txt"), row.names = F, col.names = T, sep="\t")
}



#################
#merge test and controls samples files

#SNVs

setwd("location/output/files")


files <- list.files("location/output/files", pattern="_SNV_for_MMG_test.txt")

SNV_for_MMG <- NULL

for (file in files) {
  
  data <- read.delim(file, sep="\t", header = T)
  name <- str_replace(file, "_SNV_for_MMG_test.txt", "")
  data$name <- name
  
  SNV_for_MMG <- rbind(data, SNV_for_MMG)


}

all_SNVs <- SNV_for_MMG

write.table(all_SNVs, "all_SNVs_for_MMG.txt", row.names = F, col.names = T, sep="\t")



###

#indels

setwd("location/output/files")

files <- list.files("location/output/files", pattern="_Indel_for_MMG_test.txt")


Indel_for_MMG <- NULL

for (file in files) {
  
  data <- read.delim(file, sep="\t", header = T)
  name <- str_replace(file, "_Indel_for_MMG_test.txt", "")
  data$name <- name
  
  Indel_for_MMG <- rbind(data, Indel_for_MMG)


}

all_Indels <- Indel_for_MMG

write.table(all_Indels, "all_Indels_for_MMG.txt", row.names = F, col.names = T, sep="\t")






