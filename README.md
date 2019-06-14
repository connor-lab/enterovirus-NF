# enterovirus-NF
Nextflow pipeline for _de novo_ assembly of enterovirus genomes 

### Usage:

```
nextflow run enterovirus.nf --inputfq 'reads/*_R{1,2}_001.fastq.gz' -profile ( slurm | local )

    REQUIRED

	--inputfq [Fileglob that matches your input reads, must be surrounded with 
                    'single-quotes' prevent shell expansion]

	-profile [If you are running this on a SLURM cluster, set to 'slurm' to submit jobs to 
	          the scheduler. If not, set 'local']

    OPTIONAL

	--output [Path to output directory, 'enterovirus-NF.out' by default]

	--bbnorm_target [Kmer normalization target for bbnorm, 5000 by default]
```

### Outputs:

| Filename | File content |
| :--- | :--- |
|vapor/{SAMPLENAME}.ref.out | [VAPOR](https://github.com/connor-lab/vapor) statistics |
|vapor/{SAMPLENAME}.ref.fa | [VAPOR](https://github.com/connor-lab/vapor) selected reference sequence |
|mapping/{SAMPLENAME}.assembly_mapping.sorted.bam | Reads mapped to sample assembly |
|mapping/{SAMPLENAME}.reference_mapping.sorted.bam | Reads mapped to reference sequence |
|mapping/{SAMPLENAME}.assembly_mapping.normalized.sorted.bam | Normalized reads mapped to assembly |
|assembly/{SAMPLENAME}.iva.fa | Sample assembly |
|assembly/{SAMPLENAME}.shiver.fa | Sample assembly post-processed with [Shiver](https://github.com/ChrisHIV/shiver) |



### Depends:
* Nextflow (https://www.nextflow.io/)
* Conda (Install miniconda: https://docs.conda.io/en/latest/miniconda.html)
* Conda channels bioconda and conda-forge (https://bioconda.github.io/#set-up-channels)
