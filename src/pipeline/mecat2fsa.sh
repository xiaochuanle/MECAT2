getOsMachineType() {
    OSTYPE=`uname`
    MACHINETYPE=`uname -m`

    if [ ${MACHINETYPE} == "x86_64" ]; then
        MACHINETYPE="amd64"
    fi

    if [ ${MACHINETYPE} == "Power Macintosh" ]; then	
        MACHINETYPE="ppc"
    fi

    if [ ${OSTYPE} == "SunOS" ]; then
        MACHINETYPE=`uname -p`
        if [ ${MACHINETYPE} == "sparc" ]; then
            if [ `/usr/bin/isainfo -b` == "64" ]; then
                MACHINETYPE=sparc64
            else
                MACHINETYPE=sparc32
            fi
        fi
    fi
    echo ${OSTYPE}-${MACHINETYPE}
}



CONFIG_FILE=$1
if [ "$2" =  "" ]; then
  steps=("0" "1" "2" "3")
else
  steps=(${2//,/ })
fi

basepath=$(cd `dirname $0`; pwd)
echo $basepath
PATH=$basepath/`getOsMachineType`/bin:$PATH

initialize() {
    mkdir -p $PROJECT


    dirCorrecting=$PROJECT/1-correct_raw_reads
    mkdir -p $dirCorrecting

    dirAligning=$PROJECT/2-align_corrected_reads
    mkdir -p $dirAligning


    dirAssemble=$PROJECT/3-assemble
    mkdir -p $dirAssemble

    while read line ; do
        #cat `echo ${line} | tr -d '\n'` >> $alignRRDir/rawreads.${line#*.}
        #cat  ${line}  >> $alignRRDir/rawreads.${line#*.}
        rawreads=$dirCorrecting/rawreads.${line#*.}
        break
    done < ${ONT_READ_LIST}

    corrected=$dirCorrecting/corrected.fasta
    candidates=$dirCorrecting/pm.can
    overlaps=$dirAligning/overlaps.m4
    contigs=$dirAssemble/p_ctg.fa
}


checkCode() {
   if [ $? != 0 ]; then
       echo $1
       exit 1
   fi
}

correctRawReads() {
    
    skip=1
    if [[ ! -f $rawreads ]]; then
        skip=0
    else
        while read line ; do
            if [[  ${line}  -nt $rawreads ]]; then
                skip=0
                break
            fi
        done < ${ONT_READ_LIST}
    fi

    if [[ $skip != 1 ]]; then
        rm -f $rawreads

        while read line ; do
            #cat `echo ${line} | tr -d '\n'` >> $alignRRDir/rawreads.${line#*.}
            cat  ${line}  >> $rawreads
        done < ${ONT_READ_LIST}
    fi

    if [[ ! -f $candidates || $rawreads -nt $candidates ]]; then
    
        mecat2pw -j 0 -d $rawreads  -o $candidates -w $dirCorrecting/wrk_dir -t ${THREADS}
        checkCode "Failed to align rawreads" 
    fi

    if [[ ! -f $corrected || $candidates -nt $corrected ]]; then
        mecat2cns -i 0 -t $THREADS $candidates $rawreads $corrected 
        checkCode "Failed to correct rawreads" 
    fi
}

alignCorrected() {
    if [[ ! -f $overlaps || $corrected -nt $overlaps ]]; then
        mecat2pw -j 1 -d $corrected -o $overlaps".tmp" -w $dirAligning/wrk_dir -t ${THREADS}
        checkCode "Falied to align corrected reads"
        cat $overlaps".tmp" | awk '{{ print $1+1, $2+1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12 }}' > $overlaps
        rm $overlaps".tmp"
    fi
}

assemble() {
    if [[ ! -f $contigs || $overlaps -nt $contigs ]]; then
        fsa_ol_filter $overlaps $dirAssemble/filter.m4 --thread_size=${THREADS} --output_directory=$dirAssemble --genome_size=$GENOME_SIZE $FSA_OL_FILTER_OPTIONS
        checkCode "Failed to filter overlaps"

        fsa_assemble $dirAssemble/filter.m4 --read_file=$corrected --thread_size=$THREADS --output_directory=$dirAssemble $FSA_ASSEMBLE_OPTIONS 
        checkCode "Failed to assembe reads"
    fi
}



### parse arguments
while read line ; do
    eval "${line}"
done < ${CONFIG_FILE}

initialize

for step in ${steps[@]}
do
    if [[ ${step} == "0" ]];then
        correctRawReads
    fi

    if [[ $step == "1" ]]; then
        alignCorrected
    fi

    if [[ $step == "2" ]]; then
        assemble
    fi
done


