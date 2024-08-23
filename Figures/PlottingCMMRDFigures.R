##PLOTS 

##______________________________________Plots constructed with MutationalPatterns________________________________________
library(MutationalPatterns)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

##Single mutational profile
plot_96_profile(mut_mat[,"<NUMBER>", drop = FALSE], ymax = "<VALUE>")

##Single indel profile
plot_indel_contexts(CMMRD_indel_counts[,"<NUMBER>"])

##Plot SNV types
type_occurrences <- mut_type_occurrences(CMMRD_vcfs[], ref_genome)
plot_spectrum(type_occurrences, CT = TRUE, 
              indv_points = TRUE, legend = FALSE)

##Plot indel types
plot_main_indel_contexts(CMMRD_indel_counts[])

##Plot cosine similarity heatmap
plot_cosine_heatmap(cos_sim_samples_signatures, cluster_rows = TRUE, cluster_cols = TRUE)

##Barplot cosine similarity samples
plot_original_vs_reconstructed("<ORIGINAL_PROFILE>"[,c("<SAMPLES>")], "<RECONSTRUCTED_PROFILE"[,c("<SAMPLES>")], 
                               y_intercept = 0.9)

##Plot refit results
plot_contribution(fit_res_strict$contribution[c("<SIGNATURES>"),c("<SAMPLES>")],
                  coord_flip = TRUE,
                  mode = "relative",
                  palette = c("<COLORS>"))

##Plot bootstrapping refit results
plot_bootstrapped_contribution(contri_boots_tumSig1[], 
                               mode = "relative", 
                               plot_type = "jitter")

##___________________________________________Plots constructed with ggplot_______________________________________________

library(ggplot2)

##Boxplots
box = ggplot("<DATA>", aes(x="<TUMORTYPE> or <GENE>", y="<SIGNATURE>")) +
  geom_boxplot()
Box = box+
  scale_fill_manual(name="<VALUE>", values = c("<COLORS>"))+
  geom_dotplot(binaxis = 'y', stackdir = 'center', stackratio = 0.5, dotsize = 0.5, fill = 'black') +
  ylab("Relative contribution") + 
  xlab("")+ theme_bw() + ggtitle("<SIGNATURE>")+
  theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(color = 'black'), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text = element_text(size=12), axis.title  = element_text(size=14),axis.text.x = element_text(angle=45, hjust=1))+
  scale_y_continuous(limits = c(0.00,1.20), breaks = c(0.00,0.25,0.50,0.75,1.00))

##Mutational load plots
MutBurden = read.delim("<MUT_LOAD_FILE>", header= TRUE, stringsAsFactors = TRUE)
MutBurden$total_burden = as.numeric(MutBurden$total_burden)
MutBurden_total = MutBurden[order(-MutBurden$total_burden),]

ggplot(data=MutBurden_total, aes(x=sample, y=total_burden)) +
  geom_point(stat="identity", color="orangered3", size=3)+
  ylab("Mutations/Mb")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "grey"))+
  scale_y_log10(limits=c(1,1000),labels = scales::comma)

##Plot SNV types cell line
ggplot(Cell_Line_SNV, aes(fill="<MUTATION_TYPE>", y="<VALUES>", x="<SAMPLE_NAME>"))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values = c("#30baed", "#000000", "#e98c7b", "#df2118", "#d4d1d2", "#adcc54", "#f1d1ce"))+
  ylab("SNVs per day")+theme_light()+theme(panel.grid = element_blank(), text = element_text(size=12))

##Plot indel types cell line
ggplot(Cell_Line_Indel, aes(fill="<MUTATION_TYPE>", y="<VALUES>", x="<SAMPLE_NAME>"))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values = c("#fabd6f", "#b0d28a", "#bdbdbd", "#f07e1a", "#35a137"))+
  ylab("Indels per day")+theme_light()+theme(panel.grid = element_blank(), text = element_text(size=12))

##VAF plots
ggplot(<DATAFRAME>, aes(x=<COLUMN>)) +
  xlim(0,1.0)+
  geom_density()+
  labs(x="VAF somatic variants", title="")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5), axis.title.x = element_text(size=10))

