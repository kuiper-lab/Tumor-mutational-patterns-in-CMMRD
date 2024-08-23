SCHEDULER="SLURM"
BATCHDIR="<NAME>"
PROJECTDIR="<PATH>"
INPUTDIR="$PROJECTDIR/<PATH>/$BATCHDIR/"
JOBDIR="$PROJECTDIR/<PATH>/$BATCHDIR/"
OUTPUTDIR="$PROJECTDIR/<PATH>/$BATCHDIR/"

mkdir -p "$JOBDIR"
mkdir -p "$OUTPUTDIR"

# set global varaiables
REF_GENOME="<PATH>/Homo_sapiens_assembly38.fasta"
NUM_CORES=2
NUM_THREADS=$(expr $NUM_CORES \* 2)


# VEP paths
REFERENCE_DIR="<PATH>"
DIR_CACHE="/<PATH>/vep92/cache/"
FASTA_FILE="<PATH>/Homo_sapiens_assembly38.fasta"
PLUGIN_PATH="/<PATH>/vep92/plugin/"
NUM_FORKS=$NUM_THREADS # recommended by the VEP team to use 4 forks
#NUM_CORES=4

# Java VM values
XMX="32g" #Gigabytes

#for VCF in $(ls $INPUTDIR/*chr10.vcf)
for VCF in $(ls $INPUTDIR/*chr10.vcf)
do

SAMPLE=$( basename $VCF _chr10.vcf)

echo "Generating job for sample: $SAMPLE";

declare -a chromosomes=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

for CHR in "${chromosomes[@]}"
do

# merge stats files
## ADAPTED FOR TWIST
stats_files_input_line=$(for CHR in ${chromosomes[@]}; do printf "%s %s %s\n" "-stats" $INPUTDIR/"$SAMPLE"_chr$CHR.vcf.stats "\\"; done)

# MERGE F1R2 stats
## ADAPTED FOR TWIST
f1r2_files_input_line=$(for CHR in ${chromosomes[@]}; do printf "%s %s %s\n" "-I" $INPUTDIR/"$SAMPLE"_chr$CHR\-f1r2.tar.gz "\\"; done)

# Merge Vcfs
## ADAPTED FOR TWIST
vcf_files_input_line=$(for CHR in ${chromosomes[@]}; do printf "%s %s %s\n" "-I" $INPUTDIR/"$SAMPLE"_chr$CHR.vcf "\\"; done)

done



JOBNAME="filterMutect2.$SAMPLE.sh"
OUTPUTLOG="$JOBDIR/filterMutect2.$SAMPLE.sh.out"
ERRORLOG="$JOBDIR/filterMutect2.$SAMPLE.sh.err"
WALLTIME="23:59:00"
NUMTASKS=1
NUMCPUS=$NUM_CORES
MEM="40G"
TMPSPACE="300G"
JOBFILE="$JOBDIR/filterMutect2.$SAMPLE.sh"

if [[ $SCHEDULER == "SLURM" ]]
then
    cat <<- EOF > $JOBFILE
#!/bin/bash
#SBATCH --job-name=$JOBNAME
#SBATCH --output=$OUTPUTLOG
#SBATCH --error=$ERRORLOG
#SBATCH --partition=cpu
#SBATCH --time=$WALLTIME
#SBATCH --ntasks=$NUMTASKS
#SBATCH --cpus-per-task $NUMCPUS
#SBATCH --mem=$MEM
#SBATCH --gres=tmpspace:$TMPSPACE
#SBATCH --nodes=1
#SBATCH --open-mode=append

EOF

elif [[ $SCHEDULER == "SGE" ]]
then

    cat <<- EOF > $JOBFILE
#$ -S /bin/bash
#$ -cwd
#$ -o $OUTPUTLOG
#$ -e $ERRORLOG
#$ -N $JOBNAME
#$ -l h_rt=$WALLTIME
#$ -l h_vmem=$MEM
#$ -l tmpspace=$TMPSPACE
#$ -pe threaded $NUMCPUS

EOF

else
    echo "Type of scheduler not known: $SCHEDULER"
    exit
fi




echo -e """

set -e # exit if any subcommand or pipeline returns a non-zero status
set -u # exit if any uninitialised variable is used


startTime=\$(date +%s)
echo \"startTime: \$startTime\"


set -e


# load the required modules
module load gatk/4.2.0.0
module load Java/1.8.0_60
module load ensemblvep/92
module load samtools/1.9


# merge stats
java -d64 -Xmx$XMX -Djava.io.tmpdir=\$TMPDIR -jar \$GATK4 \\
MergeMutectStats \\
$stats_files_input_line
-O $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.vcf.stats

# learn read orientation model
java -d64 -Xmx$XMX -Djava.io.tmpdir=\$TMPDIR -jar \$GATK4 \\
LearnReadOrientationModel \\
$f1r2_files_input_line
-O $OUTPUTDIR/$SAMPLE-read-orientation-model.tar.gz

# merge the vcf files
java -d64 -Xmx$XMX -Djava.io.tmpdir=\$TMPDIR -jar \$GATK4 \\
GatherVcfs \\
$vcf_files_input_line
-O $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.vcf \\
--REFERENCE_SEQUENCE $REF_GENOME

# add the filter column to the merged vcf
java -Djava.io.tmpdir=\$TMPDIR -Xmx$XMX -jar \$GATK4 \\
FilterMutectCalls \\
-R $REF_GENOME \\
-V $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.vcf \\
-O $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.vcf \\
--orientation-bias-artifact-priors $OUTPUTDIR/$SAMPLE-read-orientation-model.tar.gz \\

java -Xmx$XMX -Djava.io.tmpdir=\$TMPDIR -jar \$GATK4 \\
SelectVariants \\
-R $REF_GENOME \\
-V $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.vcf  \\
-O $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vcf  \\
--select-type-to-include SNP \\

# extract the indels from the PASS-filtered variants (--restrictAllelesTo BIALLELIC enforces that mulitallelic sites are excluded)

java -Xmx$XMX -Djava.io.tmpdir=\$TMPDIR -jar \$GATK4 \\
SelectVariants \\
-R $REF_GENOME \\
-V $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.vcf  \\
-O $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vcf  \\
--select-type-to-include INDEL \\



# annotate vcf-file
flag_fields=\"--fields Location,Allele,VARIANT_CLASS,Gene,Feature,SYMBOL,CCDS,STRAND,IMPACT,Consequence,SIFT,PolyPhen,CADD_phred,EXON,DISTANCE,COSMIC_ID,COSMIC_CNT,AC,AN,CADD_phred,gnomAD_exomes_AF,gnomAD_exomes_NFE_AF,ExAC_nonTCGA_AF,ExAC_nonTCGA_NFE_AF,gnomAD_genomes_AF,gnomAD_genomes_NFE_AF,phyloP100way_vertebrate,clinvar_rs,clinvar_clnsig,clinvar_trait,clinvar_golden_stars,Diag_Germline_Gene\"

#vep \\
#--offline \\
#--cache \\
#--dir_plugins $PLUGIN_PATH \\
#--fork $NUM_FORKS \\
#--dir_cache $DIR_CACHE \\
#--fasta $REF_GENOME \\
#--assembly GRCh38 \\
#--input_file $OUTPUTDIR/$SAMPLE.filterAnnotated.PASSonly.chr$CHR.vcf \\
#--output_file $OUTPUTDIR/$SAMPLE.filterAnnotated.PASSonly.vepAnnotated.chr$CHR.vcf \\
#--vcf $flag_fields\\
#--plugin ExAC,/<PATH>/grch38/ExAC.r0.3.1.sites.vep.vcf.gz,AC,AN \\
#--plugin dbNSFP,/<PATH>/grch38/dbNSFP3.5a.txt.gz,CADD_phred,gnomAD_exomes_AF,gnomAD_exomes_NFE_AF,ExAC_nonTCGA_AF,ExAC_nonTCGA_NFE_AF,gnomAD_genomes_AF,gnomAD_genomes_NFE_AF,phyloP100way_vertebrate,clinvar_rs,clinvar_clnsig,clinvar_trait,clinvar_golden_stars \\
#--custom /<PATH>/genome_targets/hg38_diagnostic_germline_1.0.bed.gz,Diag_Germline_Gene,bed,overlap,0,\\
#--canonical \\
#--exclude_predicted \\
#--variant_class --ccds --hgvs --symbol \\
#--pick --tsl --af_gnomad --numbers --polyphen b --protein --sift b

# annotate substitutions vcf-file
vep \\
--offline \\
--cache \\
--fork $NUM_FORKS \\
--dir_cache $DIR_CACHE \\
--dir_plugins $PLUGIN_PATH \\
--assembly GRCh38 \\
--input_file $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vcf \\
--output_file $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.vcf \\
--force_overwrite \\
--fasta $REF_GENOME \\
--vcf \\
--allele_number \\
--minimal \\
--pick \\
--pick_order tsl,appris,rank \\
--sift s \\
--polyphen s \\
--custom $REFERENCE_DIR/goNL/hg38/multisample.parents_only.info_only.vcf.gz,goNL,vcf,exact,0,AC,AN,AF \\
--custom $REFERENCE_DIR/gnomAD/gnomadGenomesHg38/genomewide/gnomad.genomes.r3.0.sites.vcf.gz,gnomADg,vcf,exact,0,AC,AN,AF,AF_male,AF_female,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_ami,AF_sas,AF_oth \\
--custom $REFERENCE_DIR/Cosmic/v91/CosmicVariants.normal.merged.vcf.gz,cosmic_v91,vcf,exact,GENOMIC_ID,COSMIC_ID,LEGACY_ID,CNT,SNP,HGVSC,HGVSG,HGVSP \\
--custom $REFERENCE_DIR/clinvar/clinvar_20200210.vcf.gz,ClinVar,vcf,exact,0,CLNALLELEID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG,AF_TGP \\
--custom $REFERENCE_DIR/phyloP/hg38/hg38.phyloP100way.bed.gz,phyloP100way,bed \\
--custom $REFERENCE_DIR/phyloP/hg38/hg38.phyloP20way.bed.gz,phyloP20way,bed \\
--custom $REFERENCE_DIR/dbSNP/All_20180418.vcf.gz,dbSNP151,vcf,exact,0,dbSNPBuildID \\
--custom $REFERENCE_DIR/annotation/cytoband/hg38.cytoband.bed.gz,cytoband,bed \\
--custom $REFERENCE_DIR/vep92/pluginReferences/gnomad.v2.1.1.hg38coordinates.pLIscores.bed.gz,gnomAD_pLI,bed \\
--custom $REFERENCE_DIR/GWAScatalog/gwas_catalog_v1.0.2-associations_e98_r2020-02-08.vcf.gz,GWAScatalogue,vcf,overlap,0,PUBMEDID,LINK,DISEASE/TRAIT,STRONGEST_SNP-RISK_ALLELE \\
--custom $REFERENCE_DIR/OMIM/generated_on_2020-03-02/genemap2.vcf.gz,OMIM,vcf,overlap,0,Gene_Name,Gene_Symbols,MIM_Number,Comments,Phenotypes \\
--plugin CADD,$REFERENCE_DIR/cadd/hg38/whole_genome_SNVs.tsv.gz,/hpc/pmc_kuiper/References/cadd/hg38/InDels.tsv.gz \\
--plugin MaxEntScan,/<PATH/maxEntScan/maxEntScan/,SWA,NCSS \\
--plugin SpliceAI,snv=/<PATH>/spliceAI/spliceai_scores.masked.snv.hg38.vcf.gz,indel=/hpc/pmc_kuiper/References/spliceAI/spliceai_scores.masked.indel.hg38.vcf.gz \\
--plugin GeneSplicer,/<PATH>/geneSplicer/geneSplicer/bin/linux/genesplicer,/hpc/pmc_kuiper/References/geneSplicer/geneSplicer/human/,context=200,cache_size=50 \\
--plugin LoF,loftee_path:$REFERENCE_DIR/vep92/plugin/,human_ancestor_fa:$REFERENCE_DIR/vep92/pluginReferences/loftee/hg38/human_ancestor.fa.gz,exonic_denovo_only:0 \\
--plugin ExACpLI,$REFERENCE_DIR/vep92/pluginReferences/ExACpLI_values.txt \\
--plugin dbscSNV,$REFERENCE_DIR/dbscSNV/hg38/dbscSNV1.1_GRCh38.txt.gz \\
--plugin dbNSFP,/<PATH>//dbNSFP/dbNSFP4.0a.gz,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,Interpro_domain

# annotate indel vcf-file
vep \\
--offline \\
--cache \\
--fork $NUM_FORKS \\
--dir_cache $DIR_CACHE \\
--dir_plugins $PLUGIN_PATH \\
--assembly GRCh38 \\
--input_file $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vcf  \\
--output_file $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.vcf \\
--force_overwrite \\
--fasta $REF_GENOME \\
--vcf \\
--allele_number \\
--minimal \\
--pick \\
--pick_order tsl,appris,rank \\
--sift s \\
--polyphen s \\
--custom $REFERENCE_DIR/goNL/hg38/multisample.parents_only.info_only.vcf.gz,goNL,vcf,exact,0,AC,AN,AF \\
--custom $REFERENCE_DIR/gnomAD/gnomadGenomesHg38/genomewide/gnomad.genomes.r3.0.sites.vcf.gz,gnomADg,vcf,exact,0,AC,AN,AF,AF_male,AF_female,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_ami,AF_sas,AF_oth \\
--custom $REFERENCE_DIR/Cosmic/v91/CosmicVariants.normal.merged.vcf.gz,cosmic_v91,vcf,exact,GENOMIC_ID,COSMIC_ID,LEGACY_ID,CNT,SNP,HGVSC,HGVSG,HGVSP \\
--custom $REFERENCE_DIR/clinvar/clinvar_20200210.vcf.gz,ClinVar,vcf,exact,0,CLNALLELEID,CLNDN,CLNDISDB,CLNREVSTAT,CLNSIG,AF_TGP \\
--custom $REFERENCE_DIR/phyloP/hg38/hg38.phyloP100way.bed.gz,phyloP100way,bed \\
--custom $REFERENCE_DIR/phyloP/hg38/hg38.phyloP20way.bed.gz,phyloP20way,bed \\
--custom $REFERENCE_DIR/dbSNP/All_20180418.vcf.gz,dbSNP151,vcf,exact,0,dbSNPBuildID \\
--custom $REFERENCE_DIR/annotation/cytoband/hg38.cytoband.bed.gz,cytoband,bed \\
--custom $REFERENCE_DIR/vep92/pluginReferences/gnomad.v2.1.1.hg38coordinates.pLIscores.bed.gz,gnomAD_pLI,bed \\
--custom $REFERENCE_DIR/GWAScatalog/gwas_catalog_v1.0.2-associations_e98_r2020-02-08.vcf.gz,GWAScatalogue,vcf,overlap,0,PUBMEDID,LINK,DISEASE/TRAIT,STRONGEST_SNP-RISK_ALLELE \\
--custom $REFERENCE_DIR/OMIM/generated_on_2020-03-02/genemap2.vcf.gz,OMIM,vcf,overlap,0,Gene_Name,Gene_Symbols,MIM_Number,Comments,Phenotypes \\
--plugin CADD,$REFERENCE_DIR/cadd/hg38/whole_genome_SNVs.tsv.gz,/hpc/pmc_kuiper/References/cadd/hg38/InDels.tsv.gz \\
--plugin MaxEntScan,/<PATH>/maxEntScan/maxEntScan/,SWA,NCSS \\
--plugin SpliceAI,snv=/<PATH>/spliceAI/spliceai_scores.masked.snv.hg38.vcf.gz,indel=/hpc/pmc_kuiper/References/spliceAI/spliceai_scores.masked.indel.hg38.vcf.gz \\
--plugin GeneSplicer,/<PATH>/geneSplicer/geneSplicer/bin/linux/genesplicer,/hpc/pmc_kuiper/References/geneSplicer/geneSplicer/human/,context=200,cache_size=50 \\
--plugin LoF,loftee_path:$REFERENCE_DIR/vep92/plugin/,human_ancestor_fa:$REFERENCE_DIR/vep92/pluginReferences/loftee/hg38/human_ancestor.fa.gz,exonic_denovo_only:0 \\
--plugin ExACpLI,$REFERENCE_DIR/vep92/pluginReferences/ExACpLI_values.txt \\
--plugin dbscSNV,$REFERENCE_DIR/dbscSNV/hg38/dbscSNV1.1_GRCh38.txt.gz \\
--plugin dbNSFP,/<PATH>//dbNSFP/dbNSFP4.0a.gz,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,Interpro_domain

# convert output VCFs to tabular
perl /<PATH>/convertVEPannotatedVCFtoTable.pl \\
-vcf $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.vcf \\
-tab $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.txt \\
-type multisample \\
-popmax=gnomADg_AF_afr,gnomADg_AF_amr,gnomADg_AF_eas,gnomADg_AF_nfe,gnomADg_AF_sas,gnomADg_AF_oth

perl /<PATH>/convertVEPannotatedVCFtoTable.pl \\
-vcf $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.vcf \\
-tab $OUTPUTDIR/$SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.txt \\
-type multisample \\
-popmax=gnomADg_AF_afr,gnomADg_AF_amr,gnomADg_AF_eas,gnomADg_AF_nfe,gnomADg_AF_sas,gnomADg_AF_oth

cd $OUTPUTDIR

bgzip $SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.vcf
bgzip $SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.vcf

tabix $SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.vcf.gz
tabix $SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.vcf.gz

md5sum $SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.vcf.gz > $SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.vcf.gz.md5
md5sum $SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.vcf.gz > $SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.vcf.gz.md5

md5sum $SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.txt > $SAMPLE-merged-mutect2Calls.filterAnnotated.snp.vepAnnotated.txt.md5
md5sum $SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.txt > $SAMPLE-merged-mutect2Calls.filterAnnotated.indel.vepAnnotated.txt.md5

cd -


##--offline \\
##--cache \\
##--database \\
##--fork $NUM_FORKS \\
##--dir_cache $DIR_CACHE \\
##--assembly GRCh38 \\
##--input_file $OUTPUTDIR/$SAMPLE.filterAnnotated.PASSonly.chr$CHR.vcf \\
##--output_file $OUTPUTDIR/$SAMPLE.filterAnnotated.PASSonly.vepAnnotated.chr$CHR.vcf \\
##--force_overwrite \\
##--fasta $REF_GENOME \\
##--merged \\
##--vcf $flag_fields\\
##--pick \\
##--pick_order tsl,appris,rank \\
##--custom /<PATH>/gnomAD/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.gz,gnomADg,vcf,exact,0,AC,AF,AN,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH \\
##--custom /<PATH>/goNL/hg38/multisample.parents_only.info_only.vcf.gz,goNL,vcf,exact,0,AC,AF,AN \\
##--plugin CADD,/<PATH>/cadd/hg38/whole_genome_SNVs.tsv.gz \\
##--dir_plugins $PLUGIN_PATH


#Retrieve and check return code
returnCode=\$?
echo \"Return code \${returnCode}\"

if [ \"\${returnCode}\" -eq \"0\" ]
then

	echo -e \"Return code is zero, process was succesfull\n\n\"

else

	echo -e \"\nNon zero return code not making files final. Existing temp files are kept for debugging purposes\n\n\"
	#Return non zero return code
	exit 1

fi


#Write runtime of process to log file
endTime=\$(date +%s)
echo \"endTime: \$endTime\"


#Source: http://stackoverflow.com/questions/12199631/convert-seconds-to-hours-minutes-seconds-in-bash

num=\$endTime-\$startTime
min=0
hour=0
day=0
if((num>59));then
    ((sec=num%60))
    ((num=num/60))
    if((num>59));then
        ((min=num%60))
        ((num=num/60))
        if((num>23));then
            ((hour=num%24))
            ((day=num/24))
        else
            ((hour=num))
        fi
    else
        ((min=num))
    fi
else
    ((sec=num))
fi
echo \"Running time: \${day} days \${hour} hours \${min} mins \${sec} secs\"

""" >> $JOBFILE

done
