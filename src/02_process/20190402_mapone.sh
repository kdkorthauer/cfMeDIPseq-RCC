pwd
echo $id
echo $GROUP
echo $PATH

##################
# initial fastqc #
##################

if [ ! -s $OUTDIR/fastqc/$GROUP/$id\_1_fastqc.html ]; then
  fastqc --version
  fastqc -o $OUTDIR/fastqc/$GROUP $id\_1.fq.gz
fi

if [ ! -s $OUTDIR/fastqc/$GROUP/$id\_2_fastqc.html ]; then
  fastqc -o $OUTDIR/fastqc/$GROUP $id\_2.fq.gz
fi

##########################
# trim_galore and fastqc #
##########################

if [ ! -s $OUTDIR/trimgalore/$GROUP/$id\_2_val_2.fq.gz ]; then
  module load python/3.5.6
  trim_galore --version
  trim_galore --fastqc --fastqc_args "-o $OUTDIR/fastqc/$GROUP" \
      -o $OUTDIR/trimgalore/$GROUP \
      --paired $id\_1.fq.gz $id\_2.fq.gz 
fi

#######################
# mapping with bowtie #
#######################

if [ ! -s $OUTDIR/bowtie2/$GROUP/$id.sam ]; then
  module load bowtie
  bowtie2 --version
  bowtie2 -x hg19 -p 2 -S $OUTDIR/bowtie2/$GROUP/$id.sam \
    -1 $OUTDIR/trimgalore/$GROUP/$id\_1_val_1.fq.gz \
    -2 $OUTDIR/trimgalore/$GROUP/$id\_2_val_2.fq.gz
fi


###############################
# sam to bam, sort, and index #
###############################

if [ ! -s $OUTDIR/sortedbam/$GROUP/$id.sorted.bam ] || \
   [ -f  $OUTDIR/sortedbam/$GROUP/$id.sorted.bam.tmp.0000.bam ]; then
  rm -f $OUTDIR/sortedbam/$GROUP/$id.sorted.bam.*
  module load samtools
  samtools view -S -b $OUTDIR/bowtie2/$GROUP/$id.sam > \
    $OUTDIR/bowtie2/$GROUP/$id.bam
  samtools sort $OUTDIR/bowtie2/$GROUP/$id.bam \
    -o $OUTDIR/sortedbam/$GROUP/$id.sorted.bam
fi


if [ ! -s $OUTDIR/sortedbam/$GROUP/$id.sorted.bam.bai ]; then
  module load samtools
  samtools index $OUTDIR/sortedbam/$GROUP/$id.sorted.bam
fi

if [ ! -s $OUTDIR/sortedbam_dup/$GROUP/$id.sorted.bam ]; then
  module load samtools
  
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

