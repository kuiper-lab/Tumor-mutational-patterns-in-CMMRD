library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(plyr)
library(dplyr)
library(stringr)
require(ggplot2)
require(ggseqlogo)


setwd("locatio/of/the/input/files")


files <- list.files(pattern="_indels_for_context.txt")

InsC5 <- NULL

for (file in files){
  
  
  setwd("location/of/the/file")
  
  f <- read.delim(file, sep="\t", header = T)
  
  f_name <- file
  f_name2 <- str_replace(f_name, "_indels_for_context.txt", "")
  
  f1 <- f %>%
    filter(Variant.type == "Insertion") %>%
    filter(Variant == "C" | Variant == "G")
  
  
  
  f1$context <- as.character(getSeq(Hsapiens,
                                    f1$Chromosome,
                                    start = f1$Start.position - 5,
                                    end = f1$Start.position + 5,
                                    strand = "+"
  ))
  
  f1$context_long <- as.character(getSeq(Hsapiens,
                                    f1$Chromosome,
                                    start = f1$Start.position - 10,
                                    end = f1$Start.position + 10,
                                    strand = "+"
  ))
  
  
  f1$fivePrime <- substring(f1$context, 1, 5)
  f1$threePrime <- substring(f1$context, 7, 11)
  
  f2 <- f1 %>%
    filter(fivePrime == "CCCCC" | fivePrime == "GGGGG" | threePrime == "CCCCC" | threePrime == "GGGGG")
  
  
  InsC5 <- rbind(InsC5, f2)
  
  
}

setwd("location/for/the/output")

df.indels <- InsC5 


ggplot() + 
  geom_logo(df.indels$context_long) +
  theme_logo()


p1 = ggseqlogo(df.indels$context_long, method = 'bits')

