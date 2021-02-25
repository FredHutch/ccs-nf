#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Containers
container__pbccs = "quay.io/biocontainers/pbccs:5.0.0--0"
container__samtools = "quay.io/biocontainers/samtools:0.1.18--hfb9b9cc_10"
container__bam2fastx = "quay.io/biocontainers/bam2fastx:1.3.1--he1c1bb9_0"
container__sequeltools = "quay.io/fhcrc-microbiome/sequeltools:latest"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/ccs-nf <ARGUMENTS>
    
    Required Arguments:
      --input_folder        Folder containing PacBio output data
      --output_folder       Folder to place analysis outputs

    Optional Arguents:
      --n_shards            Number of shards to use for extracting the CCS from CLR
                            Default: 10
      --min_passes          Minimum number of passes required for each CCS
                            Default: 3

    For more details on SequelTools, see https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03751-8.
    
    For more details on CCS, see https://ccs.how/ and https://github.com/PacificBiosciences/ccs

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
params.help = false
params.input_folder = null
params.prefix = null
params.output_folder = null
params.n_shards = 10
params.min_passes = 3
if (params.help || params.input_folder == null || params.output_folder == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 1
}

/////////////////////
// DEFINE FUNCTIONS /
/////////////////////

// Run QC with SequelTools
process qc_metrics {

  // Docker container to use
  container "${container__sequeltools}"
  label "mem_medium"
  errorStrategy 'finish'

  publishDir "${params.output_folder}/${prefix}/qc/" 
  
  input:
    tuple val(prefix), file(subreads_bam), file(subreads_pbi), file(scraps_bam), file(scraps_pbi)

  output:
    file "SequelToolsResults/*"

"""
#!/bin/bash

set -Eeuo pipefail

echo "${subreads_bam.name}" > subFiles.txt
echo "${scraps_bam.name}" > scrFiles.txt

# Run QC
SequelTools.sh \
    -v \
    -n ${task.cpus} \
    -t Q \
    -u subFiles.txt \
    -c scrFiles.txt

"""

}

// Extract CCS from CLR
process extract_ccs {

  // Docker container to use
  container "${container__pbccs}"
  label "mem_medium"
  errorStrategy 'finish'

  input:
    tuple val(prefix), file(subreads_bam), file(subreads_pbi), file(scraps_bam), file(scraps_pbi)
    each shard_ix

  output:
    tuple val(prefix), file("${prefix}.ccs.bam")

"""
#!/bin/bash

set -Eeuo pipefail

ccs ${subreads_bam} ${prefix}.ccs.bam -j ${task.cpus} --chunk ${shard_ix}/${params.n_shards} --minPasses ${params.min_passes}

"""

}

// Merge CCS BAM files across shards
process merge_ccs {

  // Docker container to use
  container "${container__samtools}"
  label "mem_medium"
  errorStrategy 'finish'

  publishDir "${params.output_folder}/${prefix}/ccs/" 
  
  input:
    tuple val(prefix), file("input.*.bam")

  output:
    tuple val(prefix), file("${prefix}.ccs.bam")

"""
#!/bin/bash

set -Eeuo pipefail

samtools merge -@${task.cpus} ${prefix}.ccs.bam input.*.bam

"""

}

// Extract CCS FASTQ from BAM
process extract_fastq {

  // Docker container to use
  container "${container__bam2fastx}"
  label "mem_medium"
  errorStrategy 'finish'

  publishDir "${params.output_folder}/${prefix}/ccs/" 
  
  input:
    tuple val(prefix), file(ccs_bam)

  output:
    tuple val(prefix), file("${prefix}*")

"""
#!/bin/bash

set -Eeuo pipefail

pbindex ${ccs_bam}
bam2fastq -o ${prefix} ${ccs_bam}

"""

}

// Start the workflow
workflow {

    // Get the input files ending with {subreads | scraps}.bam(.pbi)
    subreads_bam_ch = Channel.fromPath(
        "${params.input_folder}**.subreads.bam"
    ).map { it -> [it.name.replaceAll(/.subreads.bam/, ''), it]}
    
    subreads_pbi_ch = Channel.fromPath(
        "${params.input_folder}**.subreads.bam.pbi"
    ).map { it -> [it.name.replaceAll(/.subreads.bam.pbi/, ''), it]}

    scraps_bam_ch = Channel.fromPath(
        "${params.input_folder}**.scraps.bam"
    ).map { it -> [it.name.replaceAll(/.scraps.bam/, ''), it]}
    
    scraps_pbi_ch = Channel.fromPath(
        "${params.input_folder}**.scraps.bam.pbi"
    ).map { it -> [it.name.replaceAll(/.scraps.bam.pbi/, ''), it]}

    combined_ch = subreads_bam_ch.join(
        subreads_pbi_ch
    ).join(
        scraps_bam_ch
    ).join(
        scraps_pbi_ch
    )

    // Run QC with SequelTools
    qc_metrics(
        combined_ch
    )

    // Make a channel for the shard index
    shard_ix_ch = Channel.of(1..params.n_shards)

    // Extract the CCS into BAM format
    extract_ccs(
        combined_ch,
        shard_ix_ch
    )

    // Combine CCS BAM files across shards
    merge_ccs(
        extract_ccs.out.groupTuple()
    )

    // Extract the CCS into FASTQ format
    extract_fastq(
        merge_ccs.out
    )

}
