#!/bin/sh

# script to submit 1 job for each sample to carry out the 
# qc, trimming, and mapping

OUTDIR=/rafalab/keegan/cfMeDIPseq/out
SLURMIO=/rafalab/keegan/cfMeDIPseq/_slurm
SCRIPTDIR=/rafalab/keegan/cfMeDIPseq/src/02_process
DATDIR=/rafalab/keegan/cfMeDIPseq/data
INDEX=/rafalab/keegan/ReferenceGenomes/human/bowtie2
export BOWTIE2_INDEXES=$INDEX
RCC=novogene-raw/C202SC18123014/raw_data
BLCA=bladder
CONTROL=healthy_control
URINE=urine/hwftp.novogene.com/data_release/C202SC19040590/raw_data
PATH=$PATH:$HOME/bin/FastQC:$HOME/bin/TrimGalore-0.6.0

export PATH
export SLURMIO
export SCRIPTDIR
export DATDIR
export OUTDIR
export INDEX

# Get index
cd $INDEX
if [ ! -f hg19.zip ]; then
  wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip
  unzip hg19.zip
fi

# loop through project directories
for GROUP in RCC BLCA CONTROL URINE; 
do 
    if [ $GROUP = "URINE" ]; then
      mkdir -p $OUTDIR/fastqc/$GROUP\_CONTROL
      mkdir -p $OUTDIR/trimgalore/$GROUP\_CONTROL
      mkdir -p $OUTDIR/bowtie2/$GROUP\_CONTROL
      mkdir -p $OUTDIR/sortedbam/$GROUP\_CONTROL
      mkdir -p $OUTDIR/sortedbam_dup/$GROUP\_CONTROL
      mkdir -p $OUTDIR/fastqc/$GROUP\_RCC
      mkdir -p $OUTDIR/trimgalore/$GROUP\_RCC
      mkdir -p $OUTDIR/bowtie2/$GROUP\_RCC
      mkdir -p $OUTDIR/sortedbam/$GROUP\_RCC
      mkdir -p $OUTDIR/sortedbam_dup/$GROUP\_RCC
    else
      mkdir -p $OUTDIR/fastqc/$GROUP
      mkdir -p $OUTDIR/trimgalore/$GROUP
      mkdir -p $OUTDIR/bowtie2/$GROUP
      mkdir -p $OUTDIR/sortedbam/$GROUP
      mkdir -p $OUTDIR/sortedbam_dup/$GROUP
    fi 
    export GROUP

    if [ $GROUP = "RCC" ]; then
      cd $DATDIR/$RCC
    elif [ $GROUP = "BLCA" ]; then
      cd $DATDIR/$BLCA
    elif [ $GROUP = "CONTROL" ]; then
      cd $DATDIR/$CONTROL
    elif [ $GROUP = "URINE" ]; then
      cd $DATDIR/$URINE
    else
      echo Group not found
      exit 1
    fi

    for fq in *_1.fq.gz;
    do
         id=$(sed 's/_1.fq.gz//g' <<< $fq)
  
         if [[ $GROUP =~ "URINE" ]]; then
            echo Check urine ids
            if [ $id = "S1" ] || [ $id = "S2" ] || [ $id = "S3" ] || [ $id = "S4" ] || [ $id = "S5" ] || [ $id = "S6" ] || [ $id = "S7" ] || [ $id = "S8" ] || [ $id = "S13" ] || [ $id = "S15" ] || [ $id = "S19" ] || [ $id = "S28" ] || [ $id = "S29" ] || [ $id = "S35" ] || [ $id = "S37" ] || [ $id = "S45" ] || [ $id = "S54" ] || [ $id = "S55" ] || [ $id = "S56" ]; then
              GROUP=URINE_CONTROL
            else
              GROUP=URINE_RCC
            fi
            export GROUP
         fi
         
         echo processing $GROUP sample $id
         RUN=TRUE
         # exclude certain ids (for urine samps)
         if [[ "$id" =~ ATKIDNRL|ATKIDTMR|ATPRONRL|ATPROTMR|Undetermined|SSIP01|SSIP04|SSIP07 ]]; then
           RUN=FALSE
         fi

         if [[ "$id" =~ S57|S58 ]]; then
           RUN=FALSE
           if [[ ! $GROUP =~ URINE ]]; then
             RUN=TRUE
           fi    
         fi
   
        if [[ $GROUP =~ URINE ]]; then
           if [[ "$id" =~ S65|S66|S67|S68  ]]; then
             RUN=FALSE
           fi    
         fi

        if [ -s $OUTDIR/sortedbam_dup/$GROUP/$id.sorted.bam ]; then
          RUN=FALSE
        fi
         
        if [[ $RUN = TRUE ]]; then
           export id
           sbatch -J $id\_$GROUP -n 2 -N 1 \
            -t 2-12:00 --mem 4G \
            -o $SLURMIO/map_%x_%j.out \
            -e $SLURMIO/map_%x_%j.err \
            --wrap="sh $SCRIPTDIR/20190402_mapone.sh" 
        fi    
     done  
done

GROUP=PRCA
cd $DATDIR/nepc

export GROUP
mkdir -p $OUTDIR/fastqc/$GROUP
mkdir -p $OUTDIR/sortedbam/$GROUP
mkdir -p $OUTDIR/sortedbam_dup/$GROUP
mkdir -p $OUTDIR/trimgalore/$GROUP
mkdir -p $OUTDIR/bowtie2/$GROUP

fqfiles=$(ls ./*/*_1.fq.gz)
for fq in $fqfiles;
do
	cd $DATDIR/nepc
	id=$(sed 's/_1.fq.gz//g' <<< $(basename $fq))
	id_short=$(sed 's/_HY2FCBBXX_L005//g' <<< $id)
	echo processing $GROUP sample $id_short

  cd $id_short 
  export id
  export id_short
  
  if [ ! -s $OUTDIR/sortedbam_dup/$GROUP/$id.sorted.bam ]; then
    sbatch -J $id_short\_$GROUP -n 2 -N 1 \
      -t 2-12:00 --mem 4G \
      -o $SLURMIO/map_%x_%j.out \
      -e $SLURMIO/map_%x_%j.err \
      --wrap="sh $SCRIPTDIR/20190402_mapone.sh"
  fi
done  


# run multiqc
module load python/3.5.6
mkdir -p $OUTDIR/multiqc


for GROUP in RCC BLCA CONTROL URINE_RCC URINE_CONTROL; 
do 
	multiqc $OUTDIR/fastqc/$GROUP/ --ignore *val* -f -o $OUTDIR/multiqc \
	    -n $GROUP\_fastqc_multiqc.html
	multiqc $OUTDIR/fastqc/$GROUP/*val* -f -o $OUTDIR/multiqc \
	    -n $GROUP\_fastqc_trimmed_multiqc.html
	multiqc $OUTDIR/trimgalore/$GROUP -f -o $OUTDIR/multiqc \
	    -n $GROUP\_trimgalore_multiqc.html
	multiqc $SLURMIO/ -m bowtie2 -f -o $OUTDIR/multiqc \
	    -n botwie_multiqc.html

  if [ $GROUP = "RCC" ]; then
    multiqc $SLURMIO/*RCC* -m bowtie2 -f -o $OUTDIR/multiqc \
        -n $GROUP\_botwie_multiqc.html
  elif [ $GROUP = "BLCA" ]; then
    multiqc $SLURMIO/*BLCA* -m bowtie2 -f -o $OUTDIR/multiqc \
        -n $GROUP\_botwie_multiqc.html
  elif [ $GROUP = "CONTROL" ]; then
    multiqc $SLURMIO/*CONTROL* -m bowtie2 -f -o $OUTDIR/multiqc \
        -n $GROUP\_botwie_multiqc.html
  elif [ $GROUP = "URINE_RCC" ]; then
    multiqc $SLURMIO/*URINE_RCC* -m bowtie2 -f -o $OUTDIR/multiqc \
        -n $GROUP\_botwie_multiqc.html
  elif [ $GROUP = "URINE_CONTROL" ]; then
    multiqc $SLURMIO/*URINE_CONTROL* -m bowtie2 -f -o $OUTDIR/multiqc \
        -n $GROUP\_botwie_multiqc.html
  else
    echo "Group not found"
      #exit 1
  fi
done



