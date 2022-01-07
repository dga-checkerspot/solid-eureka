#!/usr/bin/env nextflow

params.Seqs='s3://pipe.scratch.3/resources/CHK*_R{1,2}_001.fastq.gz'
params.reffile='s3://pipe.scratch.3/resources/CHK22_ref1.fasta'



Channel
	.fromFilePairs(params.Seqs)
	.ifEmpty {error "Cannot find any reads matching: ${params.Seqs}"}
	.set { read_pairs_ch }
	
//just testing with the R1 reads now
process map {


  memory '8G'

  input:
  tuple val(pair_id), path(Seqs) from read_pairs_ch
  path ref from params.reffile
  
  output:
  file "${pair_id}.bam" into mapdir
  
  """
  minimap2 -ax sr $ref "${pair_id}_R1_001.fastq.gz" | samtools sort -o "${pair_id}.bam"
  """  

}




process pileup {

  memory '4G'
  
  input:
  path '*.bam' from mapdir.collect()
  
  output:
  file "p1_p2.mpileup" into pileupfile
  
  """
  samtools mpileup -B *.bam > p1_p2.mpileup
  """

}



process javasync {

  memory '8G'
  
  input:
  path pileup from pileupfile
  
  output:
  file "p1_p2_java.sync" into syncout
  
  """
  java -jar /popoolation2_1201/mpileup2sync.jar --input $pileup --output p1_p2_java.sync
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
