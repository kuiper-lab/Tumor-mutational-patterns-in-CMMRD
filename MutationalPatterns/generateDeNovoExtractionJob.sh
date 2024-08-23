###
SCHEDULER="SLURM"

###PARAMETERS
SAMPLE="CMMRD_ControlCohort1_12sigs"

LIBRARY_LOCATION="<PATH_TO_R_LIBRARY"
MUT_MAT_LOCATION="<FILE_PATH>/<FILE>.RData"
PERFORM_ESTIMATE=FALSE
ESTIMATE_LOCATION="<FILE_PATH>/estimate_${SAMPLE}.pdf"
PERFORM_NMF=TRUE
NMF_LOCATION="<FILE_PATH>/nmf_res_sbs_mutation_matrix_${SAMPLE}.rdata"
NUMBER_CLUSTERS=12

R_CODE="<FILE_PATH>/denovoextraction.R"
JOB_DIR="<FILE_PATH>/jobs/denovosignatureExtraction/${SAMPLE}/"
mkdir -p ${JOB_DIR}

echo -e "Creating job for De Novo Extraction"


##########
####CREATE JOB HEADER

JOBNAME="$SAMPLE-denovoextraction"
OUTPUTLOG="$JOB_DIR/$SAMPLE-denovoextraction.sh.out"
ERRORLOG="$JOB_DIR/$SAMPLE-denovoextraction.sh.err"
WALLTIME="01:59:00"
NUMTASKS="1"
NUMCPUS=1
MEM="20G"
TMPSPACE="10G"
JOBFILE="$JOB_DIR/$SAMPLE-denovoextraction.sh"

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
#$ -N $JOBNAME
#$ -o $OUTPUTLOG
#$ -e $ERRORLOG
#$ -l h_rt=$WALLTIME
#$ -l h_vmem=$MEM
#$ -l tmpspace=$TMPSPACE
#$ -cwd


EOF

else
    echo "Type of scheduler not known: $SCHEDULER"
    exit
fi




echo "


set -e # exit if any subcommand or pipeline returns a non-zero status
set -u # exit if any uninitialised variable is used


startTime=\$(date +%s)
echo \"startTime: \$startTime\"

# load the required modules
module load R/4.1.2

echo \"Starting De Novo Extraction ..\"

# Execute R script
R --slave --no-save --no-restore --no-environ \\
--args $LIBRARY_LOCATION \\
$MUT_MAT_LOCATION \\
$PERFORM_ESTIMATE \\
$ESTIMATE_LOCATION \\
$PERFORM_NMF \\
$NMF_LOCATION \\
$NUMBER_CLUSTERS \\
< ${R_CODE}

echo \"Finished De Novo Extraction\"

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

" >> $JOBFILE
