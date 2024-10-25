#plot mutational contributions CMMRD 
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(stringr)
library(forcats)


setwd("location/of/mutational/signatures")

###########################
#CMMRD tumor samples - SBS
SBS <- read_excel("mutational_signatures.xlsx", sheet=1)

SBS <- SBS %>%
  arrange(Order)

SBS$ID <- SBS$Samples

order_3 <- as.vector(SBS$ID)

SBS <- SBS %>%
  mutate(Samples2 = fct_relevel(ID, order_3)) 

decomposed_SBS_long <- gather(SBS, Signature, Contribution, `SBS1`:`SBS44`, factor_key=TRUE)


SBS_cols <- c(#"#E69F00", #SBS1-orange
              "#e0bd77", # - SBS1 new orange
              #"#56B4E9", #SBS6-light blue
              #"#999999", #SBS5-grey
              "#91c9e2", # - SBS5 new light blue
              "#c3ebdf", #SBS10a-light green
              "#ae82f5", #SBS10b-violet
              #"#000000", #SBS11-black
              "#757574", # - SBS11 new dark grey
              #"#009E73", #SBS14-green
              "#72b29e", # - new SBS14 green
              #"#D55E00", #SBS15-red
              "#d38d67", # - new SBS15 red
              "#fca59a", #SBS17a-light red
              "#c4e2ff", #SBS19-baby blue
              #"#F0E442", #SBS20-yellow
              "#f2eaa0", # - new SBS20 yellow
              #"#0072B2", #SBS26-dark blue
              "#7199aa", # - new SBS26 darker blue
              #"#CC79A7", #SBS40-pink
              "#78f5da" #SBS44-azul
)


ggplot() + geom_bar(aes(y = Contribution, x = Samples2, fill = Signature), data = decomposed_SBS_long,
                    stat="identity", color="black",  size=0.3) + 
  scale_fill_manual(values=SBS_cols) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), panel.background = element_blank()) + 
  labs(title="Relative contributions of known SBS COSMIC signatures CMMRD", 
       x =" ", y = "Relative contribution (%)")





###########################
#CMMRD tumor samples - indels
I <- read_excel("mutational_signatures.xlsx", sheet=1)


#arrange samples order
Indel <- I %>%
  arrange(Order)

Indel$ID <- Indel$Samples

order_3 <- as.vector(Indel$ID)

Indel <- Indel %>%
  mutate(Samples2 = fct_relevel(ID, order_3)) 

decomposed_Indel_long <- gather(Indel, Signature, Contribution, `ID1`:`ID83D`, factor_key=TRUE)


Indel_cols <- c("#D55E00", # - ID1 red
                "#f7cf6d", # - ID2 yellow
                "#009E73",# - ID7 green
                "#5f547a" # - ID83D purple
              
)


ggplot() + geom_bar(aes(y = Contribution, x = Samples2, fill = Signature), data = decomposed_Indel_long,
                    stat="identity", color="black",  size=0.3) + 
  scale_fill_manual(values=Indel_cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), panel.background = element_blank()) + 
  labs(title="Relative contributions of known Indel COSMIC signatures CMMRD", 
       x =" ", y = "Relative contribution (%)")



