# check that fastq header contains last 5 barcode bases from spreadsheet

# control samples
spreadsheet=/n/irizarryfs01/kkorthauer/cfMeDIPseq/data/healthy_control/14\ healtly\ patients.csv 
datadir=/n/irizarryfs01/kkorthauer/cfMeDIPseq/data/healthy_control/

nrows=`cat "$(echo $spreadsheet)" | wc -l`
((nrows++))

for row in $( seq 2 $nrows )
do
  sample=`sed -n $row\p "$spreadsheet" | awk -F "\"*,\"*" '{print $3}'`
  barcode=`sed -n $row\p "$spreadsheet" | awk -F "\"*,\"*" '{print $8}'`
  barcode=${barcode:1}
  barcode=${barcode::-1}

  header=`zcat $datadir/$sample\_1.fq.gz | head -1`

  if [[ ! $header =~ $barcode ]];
  then
    echo "Mismatch for sample $sample (barcode $barcode)"
  else 
    echo $sample OK
  fi
done


# rcc samples
spreadsheet=/n/irizarryfs01/kkorthauer/cfMeDIPseq/data/RCC/42\ RCC\ fastq\ files.csv
datadir=/n/irizarryfs01/kkorthauer/cfMeDIPseq/data/novogene-raw/C202SC18123014/raw_data/

nrows=`cat "$(echo $spreadsheet)" | wc -l`
((nrows++))

for row in $( seq 2 $nrows )
do
  sample=`sed -n $row\p "$spreadsheet" | awk -F "\"*,\"*" '{print $3}'`
  barcode=`sed -n $row\p "$spreadsheet" | awk -F "\"*,\"*" '{print $7}'`
  barcode=${barcode:1}
  barcode=${barcode::-1}

  header=`zcat $datadir/$sample\_1.fq.gz | head -1`

  if [[ ! $header =~ $barcode ]];
  then
    echo "Mismatch for sample $sample (barcode $barcode)"
  else 
    echo $sample OK
  fi
done



# do the same for urine samples
spreadsheet=/n/irizarryfs01/kkorthauer/cfMeDIPseq/data/urine/Keegan\ -\ Urine\ samples.csv
datadir=/n/irizarryfs01/kkorthauer/cfMeDIPseq/data/urine/hwftp.novogene.com/data_release/C202SC19040590/raw_data/

nrows=`cat "$(echo $spreadsheet)" | wc -l`
((nrows++))

for row in $( seq 2 $nrows )
do
  sample=`sed -n $row\p "$spreadsheet" | awk -F "\"*,\"*" '{print $1}'`
  barcode=`sed -n $row\p "$spreadsheet" | awk -F "\"*,\"*" '{print $12}'`
  barcode=${barcode:1}
  
  header=`zcat $datadir/$sample\_1.fq.gz | head -1`

  if [[ ! $header =~ $barcode ]];
  then
    echo "Mismatch for sample $sample (barcode $barcode)"
  else 
    echo $sample OK
  fi
done


# do the same for urine CONTROL samples
spreadsheet=/n/irizarryfs01/kkorthauer/cfMeDIPseq/data/urine/Keegan\ -\ Urine\ samples\ -\ Normal.csv
datadir=/n/irizarryfs01/kkorthauer/cfMeDIPseq/data/urine/hwftp.novogene.com/data_release/C202SC19040590/raw_data/

nrows=`cat "$(echo $spreadsheet)" | wc -l`
((nrows++))

for row in $( seq 2 $nrows )
do
  sample=`sed -n $row\p "$spreadsheet" | awk -F "\"*,\"*" '{print $1}'`
  barcode=`sed -n $row\p "$spreadsheet" | awk -F "\"*,\"*" '{print $10}'`
  barcode=${barcode:1}
  
  header=`zcat $datadir/$sample\_1.fq.gz | head -1`

  if [[ ! $header =~ $barcode ]];
  then
    echo "Mismatch for sample $sample (barcode $barcode)"
  else 
    echo $sample OK
  fi
done