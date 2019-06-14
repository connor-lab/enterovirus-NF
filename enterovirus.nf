#!/usr/bin/env nextflow

params.ref = 'refs/D68_ref_genomes.fa'
params.shiverinitdir = 'refs/D68_shiver_init'
params.bbnorm_target = '5000'
params.output = 'enterovirus-NF.out'

Channel.fromFilePairs(params.inputfq, flat: true)
       .set{ inputReadsTrim }

Channel.fromPath(params.ref)
       .into{ refSeq ; installVapor}

Channel.fromPath(params.shiverinitdir, type: 'dir' )
       .set{ shiverInitDir }


process TRIM_READS {

    tag { dataset_id }

    publishDir "${params.output}/trimming", pattern: "${dataset_id}.trimmed_{1,2}.fq.gz", mode: 'copy'

    conda 'trim-galore'

    input:
    set dataset_id, file(fwd), file(rev) from inputReadsTrim

    output:
    set dataset_id, file("${dataset_id}.trimmed_1.fq.gz"), file("${dataset_id}.trimmed_2.fq.gz") into trimmedReadsVapor, trimmedReadsClean, trimmedReadsMapping, trimmedReadsMappingRef, trimmedReadsVariantCalling, trimmedReadsShiver

    script:
    """
    trim_galore --paired $fwd $rev
    mv ${dataset_id}_R1_001_val_1.fq.gz ${dataset_id}.trimmed_1.fq.gz
    mv ${dataset_id}_R2_001_val_2.fq.gz ${dataset_id}.trimmed_2.fq.gz
    """
}

process INSTALL_VAPOR {
    conda 'numpy'

    input:
    file(ref) from installVapor

    output:
    file("$ref") into refFile

    script:
    """
    pip install git+https://github.com/connor-lab/vapor
    """

}

process VAPOR {

    tag { dataset_id }

    publishDir "${params.output}/vapor", pattern: "*.ref.{out,fa}", mode: 'copy'

    conda 'numpy'
 
    memory { 4.GB * task.attempt }

    errorStrategy { task.attempt > 3 ? 'ignore' : 'retry' }

    maxRetries 3

    input:
    file(ref) from refFile.collect()
    set dataset_id, file(fwd), file(rev) from trimmedReadsVapor

    output:
    set dataset_id, file("${dataset_id}.ref.fa") into vaporReference, vaporReferenceMapping, vaporReferenceVariantCalling
    file("*.ref.out") 

    script:
    """
    vapor.py -o ${dataset_id}.ref -m 0.01 -c 2 -f 0.4 -fa $ref -fq $fwd
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
    set dataset_id, file(fwd), file(rev), file("*") from trimmedReadsClean.combine(vaporReferenceIndex, by: 0)

    output:
    set dataset_id, file('*.clean_1.fq.gz'), file('*.clean_2.fq.gz') into cleanReads, cleanReadsMapping

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
    bbnorm.sh t=${task.cpus} prefilter=t min=20 target=${params.bbnorm_target} in=${fwd} in2=${rev} out=${dataset_id}.sampled_1.fq.gz out2=${dataset_id}.sampled_2.fq.gz
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
    set dataset_id, file("${dataset_id}.iva.fa") optional true into assemblyShiver
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

process SHIVER {
    tag { dataset_id }

    publishDir "${params.output}/assembly", pattern: "${dataset_id}.shiver.fa", mode: 'copy'

    conda 'shiver'

    input:
    set dataset_id, file(assembly), file(forward), file(reverse), file(shiverinit) from assemblyShiver.join(trimmedReadsShiver).combine(shiverInitDir)

    output:
    set dataset_id, file("${dataset_id}.shiver.fa") into assemblyPassCoverage, assemblyPassNormalization 

    script:
    """
    zcat $forward | awk '{if (NR%4 == 1) {print \$1 "/1"} else print}' | gzip > ${dataset_id}_shiverRenamed_1.fq.gz
    zcat $reverse | awk '{if (NR%4 == 1) {print \$1 "/2"} else print}' | gzip > ${dataset_id}_shiverRenamed_2.fq.gz
    shiver_align_contigs.sh ${shiverinit} ${shiverinit}/config.sh ${assembly} ${dataset_id}
    if [ -f ${dataset_id}_cut_wRefs.fasta ]; then
      shiver_map_reads.sh ${shiverinit} ${shiverinit}/config.sh ${assembly} ${dataset_id} ${dataset_id}.blast ${dataset_id}_cut_wRefs.fasta ${dataset_id}_shiverRenamed_1.fq.gz ${dataset_id}_shiverRenamed_2.fq.gz
    else
      shiver_map_reads.sh ${shiverinit} ${shiverinit}/config.sh ${assembly} ${dataset_id} ${dataset_id}.blast ${dataset_id}_raw_wRefs.fasta ${dataset_id}_shiverRenamed_1.fq.gz ${dataset_id}_shiverRenamed_2.fq.gz
    fi
    seqtk seq -l0 ${dataset_id}_remap_consensus_MinCov_15_30.fasta | head -n2 | sed '/>/!s/-//g' | sed 's/\\?/N/g' | sed 's/_remap_consensus//g' | seqtk seq -l80 > ${dataset_id}.shiver.fa
    """
}

    
process REMAP_READS_FOR_COVERAGE_ASSEMBLY {
    tag { dataset_id }

    publishDir "${params.output}/mapping", pattern: "${dataset_id}.shiver.assembly_mapping.sorted.bam", mode: 'copy'

    conda 'minimap2 samtools'

    cpus 4 

    input:
    set dataset_id, file(assembly), file(fwd), file(rev) from assemblyPassCoverage.join(cleanReadsMapping, by: 0)

    output:
    file("${dataset_id}.shiver.assembly_mapping.sorted.bam")

    script:
    """
    minimap2 -ax sr $assembly $fwd $rev | samtools view -F4 -b - | samtools sort -o ${dataset_id}.shiver.assembly_mapping.sorted.bam -
    """
}

/*
process REMAP_READS_FOR_COVERAGE_REF {
    tag { dataset_id }

    publishDir "${params.output}/mapping", pattern: "${dataset_id}.reference_mapping.sorted.bam", mode: 'copy'

    conda 'minimap2 samtools'

    cpus 4

    input:
    set dataset_id, file(assembly), file(fwd), file(rev) from vaporReferenceMapping.join(trimmedReadsMappingRef, by: 0)

    output:
    file("${dataset_id}.reference_mapping.sorted.bam")

    script:
    """
    minimap2 -ax sr $assembly $fwd $rev | samtools view -F4 -b - | samtools sort -o ${dataset_id}.reference_mapping.sorted.bam -
    """
}
*/

process REMAP_NORMALIZED_READS {
    tag { dataset_id }

    publishDir "${params.output}/mapping", pattern: "${dataset_id}.shiver.assembly_mapping.normalized.sorted.bam", mode: 'copy'

    conda 'minimap2 samtools'

    cpus 4

    input:
    set dataset_id, file(assembly), file(fwd), file(rev) from assemblyPassNormalization.join(subsampledReadsMapping, by: 0)

    output:
    file("${dataset_id}.shiver.assembly_mapping.normalized.sorted.bam")

    script:
    """
    minimap2 -ax sr $assembly $fwd $rev | samtools view -F4 -b - | samtools sort -o ${dataset_id}.shiver.assembly_mapping.normalized.sorted.bam -
    """
}


process MAP_READS_VARIANT_CALLING {

    tag { dataset_id }

    conda 'minimap2 samtools'
    
    cpus 4

    input:
    set dataset_id, file(ref), file(fwd), file(rev) from vaporReferenceVariantCalling.join(trimmedReadsVariantCalling, by: 0)

    output:
    set dataset_id, file("${dataset_id}.bam"), file("${ref}") into variantCallingBAM

    script:
    """
    minimap2 -t ${task.cpus} -x sr -a $ref $fwd $rev | samtools view -@ 2 -b | samtools sort -@ 2 -o ${dataset_id}.bam
    """
}


process VARIANT_CALLING {
    tag { dataset_id }

    publishDir "${params.output}/variant-calling", pattern: "${dataset_id}.consensus.fa", mode: 'copy'
    publishDir "${params.output}/variant-calling", pattern: "${dataset_id}.vcf", mode: 'copy'

    conda 'samtools bcftools tabix lofreq'

    input:
    set dataset_id, file(bam), file(ref) from variantCallingBAM

    output:
    file("*.consensus.fa") 
    file("${dataset_id}.vcf")

    script:
    """
    samtools faidx $ref
    bcftools mpileup --redo-BAQ --min-MQ 20 -Ou -f ${ref} ${bam} | bcftools call --ploidy 1 -mv -Ov | bcftools view -q 0.1:nref -Oz -o ${dataset_id}.variants.vcf.gz
    zcat ${dataset_id}.variants.vcf.gz > ${dataset_id}.vcf
    tabix ${dataset_id}.variants.vcf.gz
    bcftools consensus -I -f ${ref} ${dataset_id}.variants.vcf.gz --output ${dataset_id}.consensus.fa
    sed -i "s/>/>${dataset_id}.consensus/" ${dataset_id}.consensus.fa
    """
}
