#!/usr/bin/env nextflow

params.Seqs='s3://pipe.scratch.3/resources/CHK*_R{1,2}_001.fastq.gz'
params.reffile='s3://pipe.scratch.3/resources/CHK22_ref1.fasta'



Channel
	.fromFilePairs(params.Seqs)
	.ifEmpty {error "Cannot find any reads matching: ${params.Seqs}"}
	.set { read_pairs_ch }
	
process map {


  memory '8G'

  input:
  tuple val(pair_id), path(Seqs) from read_pairs_ch
  path ref from params.reffile
  
  output:
  file "${pair_id}.sort.bam" into mapdir
  
  """
  bwa index $ref
  bwa aln $ref "${pair_id}_R1_001.fastq.gz" "${pair_id}_R2_001.fastq" > "${pair_id}.sai"
  bwa samse $ref "${pair_id}.sai" "${pair_id}_R1_001.fastq.gz" "${pair_id}_R2_001.fastq" > "${pair_id}.sam"
  samtools view -bS "${pair_id}.sam" > "${pair_id}.bam"
  samtools sort "${pair_id}.bam" -o "${pair_id}.sort.bam"
  """  

}




process pileup {

  memory '4G'
  
  input:
  path '*.bam' from mapdir.collect()
  
  output:
  file "p1_p2.mpileup" into pileupfile
  file "dir.txt" into bamcheck
  file "1.bam" into bam1
  file "2.bam" into bam2
  
  """
  ls *.bam > dir.txt
  samtools mpileup -B *.bam > p1_p2.mpileup
  """

}



process javasync {

  memory '8G'
  
  input:
  path pileup from pileupfile
  
  output:
  file "p1_p2.sync" into syncout
  
  """
  perl /popoolation2_1201/mpileup2sync.pl --input $pileup --output p1_p2.sync --fastq-type sanger 
  """

}

process fst {


  memory '8G'
  
  input:
  path syncfile from syncout
  
  output:
  file "p1_p2_w500.fst" into fstoutput
  
  
  """
  perl /popoolation2_1201/fst-sliding.pl --input $syncfile --output p1_p2_w500.fst --min-count 6 --min-coverage 6 --max-coverage 1000 --min-covered-fraction 0.75 --window-size 500 --step-size 500 --pool-size 500
  """

}
