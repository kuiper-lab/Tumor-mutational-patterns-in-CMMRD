
SCHEDULER="SLURM"
BATCHDIR=""
PROJECTDIR=""
JOBDIR="$PROJECTDIR//$BATCHDIR/"
OUTPUTDIR="$PROJECTDIR//$BATCHDIR/"
SAMPLEFILE="$PROJECTDIR/samplelists/sampleListMutect2.$BATCHDIR.Converted.txt"
BEDDIR=""
CRAM_EXT=""

mkdir -p "$JOBDIR"
mkdir -p "$OUTPUTDIR"

# set global variables
REF_GENOME="REFDIR/Homo_sapiens_assembly38.fasta"
NUM_CORES=4
#NUM_THREADS=$(expr $NUM_CORES \* 2)
NUM_THREADS=4

#FUNCTION: Join array with local IFS, to avoid preserving non-default IFS
join_arr() {
    local IFS="$1"
    shift
    echo "$*"
}

while read line
do

    #Read line in array
    declare -a Array
    Array=($line)

    #Print all BAMs to be processed on a new line
    BAMstoProcess=$(for SAMPLE in ${Array[@]}; do printf "%s %s %s\n" "-I" $PROJECTDIR//${SAMPLE}${CRAM_EXT} "\\"; done)

    #Set samples for job name usage
    SAMPLES=$(join_arr "." "${Array[@]}")

    #According to input file format the first sampleID is the normal, the others tumors
    NORMAL=${Array[0]}

    echo "Generating jobs for samples: $SAMPLES"

#create job body
for BED in $( ls $BEDDIR/*.bed )
do

BEDFILE=$BED
BEDNAME=$( basename $BED .bed )
JOBNAME="mutect2.$SAMPLES.$BEDNAME.sh"
OUTPUTLOG="$JOBDIR/mutect2.$SAMPLES.$BEDNAME.sh.out"
ERRORLOG="$JOBDIR/mutect2.$SAMPLES.$BEDNAME.sh.err"
WALLTIME="23:59:00"
NUMTASKS=1
NUMCPUS=$NUM_CORES
MEM="16G"
TMPSPACE="300G"
JOBFILE="$JOBDIR/mutect2.$SAMPLES.$BEDNAME.sh"

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

#run gatk
java -Djava.io.tmpdir=\$TMPDIR -Xmx12g -jar \$GATK4 \\
Mutect2 \\
-R $REF_GENOME \\
--native-pair-hmm-threads $NUM_THREADS \\
$BAMstoProcess
-normal $NORMAL \\
-L $BEDFILE \\
--f1r2-tar-gz $OUTPUTDIR/$SAMPLES.$BEDNAME-f1r2.tar.gz \\
--dont-use-soft-clipped-bases true \\
--max-mnp-distance 0 \\
-O $OUTPUTDIR/$SAMPLES.$BEDNAME.vcf


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

done< $SAMPLEFILE
