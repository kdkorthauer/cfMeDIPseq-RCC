cd $cwd
pwd
echo $id
echo $fq
echo $OUTDIR
echo $GROUP
echo $PATH

module load miniconda3/4.6.14
eval "$(conda shell.bash hook)"
conda activate cfmedip

BOWTIE2_INDEXES=$INDEX

##################
# initial fastqc #
##################

if [ ! -s $OUTDIR/fastqc/$GROUP/$id\_1_fastqc.html ]; then
  fastqc --version
  fastqc -o $OUTDIR/fastqc/$GROUP $fq\_1.fq.gz

  if [[ $GROUP =~ "JAN2020" ]]; then
    mv $OUTDIR/fastqc/$GROUP/$(basename $fq)\_1_fastqc.html $OUTDIR/fastqc/$GROUP/$id\_1_fastqc.html
    mv $OUTDIR/fastqc/$GROUP/$(basename $fq)\_1_fastqc.zip $OUTDIR/fastqc/$GROUP/$id\_1_fastqc.zip
  fi
fi

if [ ! -s $OUTDIR/fastqc/$GROUP/$id\_2_fastqc.html ]; then
  fastqc -o $OUTDIR/fastqc/$GROUP $fq\_2.fq.gz

  if [[ $GROUP =~ "JAN2020" ]]; then
    mv $OUTDIR/fastqc/$GROUP/$(basename $fq)\_2_fastqc.html $OUTDIR/fastqc/$GROUP/$id\_2_fastqc.html
    mv $OUTDIR/fastqc/$GROUP/$(basename $fq)\_2_fastqc.zip $OUTDIR/fastqc/$GROUP/$id\_2_fastqc.zip
  fi
fi

##########################
# trim_galore and fastqc #
##########################

if [ ! -s $OUTDIR/trimgalore/$GROUP/$id\_2_val_2.fq.gz ]; then
  module load python/3.6.8
  trim_galore --version
  trim_galore --fastqc --fastqc_args "-o $OUTDIR/fastqc/$GROUP" \
      -o $OUTDIR/trimgalore/$GROUP \
      --paired $fq\_1.fq.gz $fq\_2.fq.gz

  if [[ $GROUP =~ "JAN2020" ]]; then
    mv $OUTDIR/trimgalore/$GROUP/$(basename $fq)\_1_val_1.fq.gz $OUTDIR/trimgalore/$GROUP/$id\_1_val_1.fq.gz
    mv $OUTDIR/trimgalore/$GROUP/$(basename $fq)\_2_val_2.fq.gz $OUTDIR/trimgalore/$GROUP/$id\_2_val_2.fq.gz
  fi 
fi

#######################
# mapping with bowtie #
#######################

if [ ! -s $OUTDIR/bowtie2/$GROUP/$id.bam ]; then
  bowtie2 --version
  bowtie2 -x hg19 -p 2 -S $OUTDIR/bowtie2/$GROUP/$id.sam \
    -1 $OUTDIR/trimgalore/$GROUP/$id\_1_val_1.fq.gz \
    -2 $OUTDIR/trimgalore/$GROUP/$id\_2_val_2.fq.gz

  samtools view -S -b $OUTDIR/bowtie2/$GROUP/$id.sam > \
    $OUTDIR/bowtie2/$GROUP/$id.bam

  rm $OUTDIR/bowtie2/$GROUP/$id.sam
fi


if [ -s $OUTDIR/bowtie2/$GROUP/$id.bam ]; then
  if [ -s $OUTDIR/bowtie2/$GROUP/$id.sam ]; then
    rm $OUTDIR/bowtie2/$GROUP/$id.sam
  fi
fi

###############################
# sam to bam, sort, and index #
###############################

if [ ! -s $OUTDIR/sortedbam/$GROUP/$id.sorted.bam ] || \
   [ -f  $OUTDIR/sortedbam/$GROUP/$id.sorted.bam.tmp.0000.bam ]; then
  rm -f $OUTDIR/sortedbam/$GROUP/$id.sorted.bam.*
  samtools sort $OUTDIR/bowtie2/$GROUP/$id.bam \
    -o $OUTDIR/sortedbam/$GROUP/$id.sorted.bam
fi


if [ ! -s $OUTDIR/sortedbam/$GROUP/$id.sorted.bam.bai ]; then
  samtools index $OUTDIR/sortedbam/$GROUP/$id.sorted.bam
fi

if [ ! -s $OUTDIR/sortedbam_dup/$GROUP/$id.sorted.bam ]; then
  # cleanup / remove extraneous files from previous runs
  rm -f $OUTDIR/sortedbam_dup/$GROUP/$id.positionsort.bam.*
  rm -f $OUTDIR/sortedbam_dup/$GROUP/$id.namesort.bam.*
  
  # The first sort can be omitted if the file is already name ordered
  samtools sort -n -o $OUTDIR/sortedbam_dup/$GROUP/$id.namesort.bam $OUTDIR/sortedbam/$GROUP/$id.sorted.bam

  # Add ms and MC tags for markdup to use later
  samtools fixmate -m $OUTDIR/sortedbam_dup/$GROUP/$id.namesort.bam $OUTDIR/sortedbam_dup/$GROUP/$id.fixmate.bam

  # Markdup needs position order
  samtools sort -o $OUTDIR/sortedbam_dup/$GROUP/$id.positionsort.bam $OUTDIR/sortedbam_dup/$GROUP/$id.fixmate.bam

  # Finally mark and remove duplicates
  samtools markdup -r -s $OUTDIR/sortedbam_dup/$GROUP/$id.positionsort.bam $OUTDIR/sortedbam_dup/$GROUP/$id.sorted.bam

  # and index
  samtools index $OUTDIR/sortedbam_dup/$GROUP/$id.sorted.bam

  # cleanup remove extraneous files
  rm $OUTDIR/sortedbam_dup/$GROUP/$id.namesort.bam
  rm $OUTDIR/sortedbam_dup/$GROUP/$id.fixmate.bam
  rm $OUTDIR/sortedbam_dup/$GROUP/$id.positionsort.bam
fi

