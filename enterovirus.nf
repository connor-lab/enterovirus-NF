#!/usr/bin/env nextflow

params.ref = 'refs/D68_ref_genomes.fa'
params.bbnorm_target = '500'
params.output = 'enterovirus-NF.out'

Channel.fromFilePairs(params.inputfq, flat: true)
       .into{ inputReadsVapor; inputReadsClean; inputReadsMapping; inputReadsMappingRef }

Channel.fromPath(params.ref)
       .into{ refSeq ; refFile}


process VAPOR {

    tag { dataset_id }

    publishDir "${params.output}/vapor", pattern: "*.ref.{out,fa}", mode: 'copy'

    memory { 16.GB * task.attempt }

    errorStrategy 'retry'

    maxRetries 3

    input:
    file(ref) from refFile.collect()
    set dataset_id, file(fwd), file(rev) from inputReadsVapor

    output:
    set dataset_id, file("${dataset_id}.ref.fa") into vaporReference, vaporReferenceMapping
    file("*.ref.out") 

    script:
    """
    vapor.py -o ${dataset_id}.ref -fa $ref -fq $fwd
    """
}
   
    

process INDEX_REFERENCE {

    tag { dataset_id }
 
    conda 'bowtie2'
  
    input:
    set dataset_id, file(ref) from vaporReference

    output:
    set dataset_id, file("${dataset_id}*") into vaporReferenceIndex

    script:
    """
    bowtie2-build $ref ${dataset_id}.ref
    """
}


process CLEAN_INPUT_READS {
    tag { dataset_id }

    conda 'bowtie2 samtools picard'

    cpus 4

    input:
    set dataset_id, file(fwd), file(rev), file("*") from inputReadsClean.combine(vaporReferenceIndex, by: 0)

    output:
    set dataset_id, file('*.clean_1.fq.gz'), file('*.clean_2.fq.gz') into cleanReads

    script:
    """
    bowtie2 --very-sensitive-local --threads ${task.cpus} -x ${dataset_id}.ref -1 $fwd -2 $rev | samtools view -F4 -b - | samtools sort - | picard SamToFastq VALIDATION_STRINGENCY=LENIENT F=${dataset_id}.clean_1.fq.gz F2=${dataset_id}.clean_2.fq.gz I=/dev/stdin
    """
}

process SUBSAMPLE_READS {
    tag { dataset_id }

    conda 'bbmap'

    cpus 8

    input:
    set dataset_id, file(fwd), file(rev) from cleanReads

    output:
    set dataset_id, file('*.sampled_1.fq.gz'), file('*.sampled_2.fq.gz') into subsampledReadsAssembly, subsampledReadsMapping

    script:
    """
    #shuffle.sh -Xmx16g in=$fwd in2=$rev out=shuffled_1.fq.gz out2=shuffled_2.fq.gz
    bbnorm.sh t=${task.cpus} prefilter=t min=10 target=${params.bbnorm_target} in=${fwd} in2=${rev} out=${dataset_id}.sampled_1.fq.gz out2=${dataset_id}.sampled_2.fq.gz
    """
}


process ASSEMBLY {
    tag { dataset_id }

    publishDir "${params.output}/assembly", pattern: "*.fa", mode: 'copy'

    conda 'iva trimmomatic'

    cpus 4

    input:
    set dataset_id, file(fwd), file(rev) from subsampledReadsAssembly

    output:
    set dataset_id, file("${dataset_id}.iva.fa") optional true into assemblyPassCoverage, assemblyPassNormalization
    set dataset_id, file("${dataset_id}.assembly.failed") optional true into assemblyFail

    script:
    """
    if iva --threads ${task.cpus} -f $fwd -r $rev iva_assembly &> ${dataset_id}.assembly.failed ; then
      mv iva_assembly/contigs.fasta ${dataset_id}.iva.fa
    else
      sed -i "1s/^/${dataset_id}\t/" ${dataset_id}.assembly.failed
    fi
    """
}


process REMAP_READS_FOR_COVERAGE_ASSEMBLY {
    tag { dataset_id }

    publishDir "${params.output}/mapping", pattern: "${dataset_id}.assembly_mapping.sorted.bam", mode: 'copy'

    conda 'minimap2 samtools'

    cpus 4 

    input:
    set dataset_id, file(assembly), file(fwd), file(rev) from assemblyPassCoverage.join(inputReadsMapping, by: 0)

    output:
    file("${dataset_id}.assembly_mapping.sorted.bam")

    script:
    """
    minimap2 -ax sr $assembly $fwd $rev | samtools view -F4 -b - | samtools sort -o ${dataset_id}.assembly_mapping.sorted.bam -
    """
}

process REMAP_READS_FOR_COVERAGE_REF {
    tag { dataset_id }

    publishDir "${params.output}/mapping", pattern: "${dataset_id}.reference_mapping.sorted.bam", mode: 'copy'

    conda 'minimap2 samtools'

    cpus 4

    input:
    set dataset_id, file(assembly), file(fwd), file(rev) from vaporReferenceMapping.join(inputReadsMappingRef, by: 0)

    output:
    file("${dataset_id}.reference_mapping.sorted.bam")

    script:
    """
    minimap2 -ax sr $assembly $fwd $rev | samtools view -F4 -b - | samtools sort -o ${dataset_id}.reference_mapping.sorted.bam -
    """
}

process REMAP_NORMALIZED_READS {
    tag { dataset_id }

    publishDir "${params.output}/mapping", pattern: "${dataset_id}.assembly_mapping.normalized.sorted.bam", mode: 'copy'

    conda 'minimap2 samtools'

    cpus 4

    input:
    set dataset_id, file(assembly), file(fwd), file(rev) from assemblyPassNormalization.join(subsampledReadsMapping, by: 0)

    output:
    file("${dataset_id}.assembly_mapping.normalized.sorted.bam")

    script:
    """
    minimap2 -ax sr $assembly $fwd $rev | samtools view -F4 -b - | samtools sort -o ${dataset_id}.assembly_mapping.normalized.sorted.bam -
    """
}

