for fastq in data/*fastq.gz
do
  echo "bash process.sh $fastq" | qsub -V -l h_vmem=8g -cwd -N map.`basename $fastq`
done
