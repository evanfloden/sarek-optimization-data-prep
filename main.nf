#!/usr/bin/env nextflow
/*
 * Data Preparation Workflow for Sarek Optimization
 * Extracts chr21 subset from GIAB HG002 data for benchmarking
 *
 * Reviewed by Seqera AI - 2026-01-22
 */

nextflow.enable.dsl = 2

params.outdir = 's3://scidev-playground-eu-west-2/sarek-optimization/data'
params.chromosome = 'chr21'

// GIAB HG002 data URLs (GRCh38)
params.source_bam = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam'
params.source_bai = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam.bai'

// GIAB HG002 truth VCF (NISTv4.2.1 GRCh38)
params.truth_vcf = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz'
params.truth_vcf_tbi = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi'
params.truth_bed = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed'

// Validate chromosome parameter (using regex for reliability)
if (!(params.chromosome ==~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/)) {
    error "Invalid chromosome: ${params.chromosome}. Must be chr1-22, chrX, chrY, or chrM"
}

process EXTRACT_CHR21_BAM {
    tag "HG002_${params.chromosome}"
    publishDir "${params.outdir}/input", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    cpus 4
    memory '8 GB'
    time '4h'

    input:
    val chromosome

    output:
    path "HG002_${chromosome}.bam", emit: bam
    path "HG002_${chromosome}.bam.bai", emit: bai

    script:
    """
    # Extract chr21 reads from remote BAM (streaming with HTTP range requests)
    samtools view -@ ${task.cpus} -b -h \\
        "${params.source_bam}" \\
        ${chromosome} \\
        -o HG002_${chromosome}.bam

    # Validate BAM was created successfully
    if [ ! -s HG002_${chromosome}.bam ]; then
        echo "Error: Failed to extract BAM data for ${chromosome}" >&2
        exit 1
    fi

    # Index the subset BAM
    samtools index -@ ${task.cpus} HG002_${chromosome}.bam

    # Quick validation
    echo "Extracted reads: \$(samtools view -c HG002_${chromosome}.bam)"
    """
}

process EXTRACT_CHR21_VCF {
    tag "truth_${params.chromosome}"
    publishDir "${params.outdir}/truth", mode: 'copy'

    container 'quay.io/biocontainers/bcftools:1.21--h3a4d415_1'

    cpus 2
    memory '4 GB'
    time '1h'

    input:
    val chromosome

    output:
    path "HG002_${chromosome}_truth.vcf.gz", emit: vcf
    path "HG002_${chromosome}_truth.vcf.gz.tbi", emit: tbi

    script:
    """
    # Extract chromosome from truth VCF (bcftools handles remote URLs and index automatically)
    bcftools view \\
        -r ${chromosome} \\
        "${params.truth_vcf}" \\
        -Oz -o HG002_${chromosome}_truth.vcf.gz

    if [ ! -s HG002_${chromosome}_truth.vcf.gz ]; then
        echo "Error: Failed to extract VCF data for ${chromosome}" >&2
        exit 1
    fi

    # Index the subset VCF
    bcftools index -t HG002_${chromosome}_truth.vcf.gz

    # Quick validation
    echo "Extracted variants: \$(bcftools view -H HG002_${chromosome}_truth.vcf.gz | wc -l)"
    """
}

process EXTRACT_CHR21_BED {
    tag "bed_${params.chromosome}"
    publishDir "${params.outdir}/truth", mode: 'copy'

    container 'community.wave.seqera.io/library/bedtools:2.31.1--fb43a408b35255d6'

    cpus 1
    memory '2 GB'
    time '30m'

    input:
    val chromosome

    output:
    path "HG002_${chromosome}_highconf.bed", emit: bed

    script:
    """
    # Download full BED and filter to chromosome
    curl -sL "${params.truth_bed}" | \\
        grep "^${chromosome}\\s" > HG002_${chromosome}_highconf.bed || true

    # Check if we got any regions
    if [ ! -s HG002_${chromosome}_highconf.bed ]; then
        echo "Warning: No high-confidence regions found for ${chromosome}" >&2
        echo "Creating empty BED file"
        touch HG002_${chromosome}_highconf.bed
    fi

    # Report stats
    echo "High-confidence regions: \$(wc -l < HG002_${chromosome}_highconf.bed)"
    """
}

process BAM_TO_FASTQ {
    tag "fastq_${params.chromosome}"
    publishDir "${params.outdir}/input", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    cpus 4
    memory '8 GB'
    time '2h'

    input:
    path bam
    path bai

    output:
    path "HG002_${params.chromosome}_R1.fastq.gz", emit: r1
    path "HG002_${params.chromosome}_R2.fastq.gz", emit: r2

    script:
    """
    # Convert BAM to paired FASTQ
    samtools collate -@ ${task.cpus} -O -u ${bam} | \\
        samtools fastq -@ ${task.cpus} \\
        -1 HG002_${params.chromosome}_R1.fastq.gz \\
        -2 HG002_${params.chromosome}_R2.fastq.gz \\
        -0 /dev/null -s /dev/null -n

    # Validate outputs
    if [ ! -s HG002_${params.chromosome}_R1.fastq.gz ] || [ ! -s HG002_${params.chromosome}_R2.fastq.gz ]; then
        echo "Error: Failed to convert BAM to FASTQ" >&2
        exit 1
    fi

    # Report stats
    echo "R1 reads: \$(zcat HG002_${params.chromosome}_R1.fastq.gz | wc -l | awk '{print \$1/4}')"
    echo "R2 reads: \$(zcat HG002_${params.chromosome}_R2.fastq.gz | wc -l | awk '{print \$1/4}')"
    """
}

process COLLECT_STATS {
    tag "stats"
    publishDir "${params.outdir}", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    cpus 1
    memory '2 GB'
    time '30m'

    input:
    path bam
    path bai
    path r1
    path r2
    path vcf
    path bed

    output:
    path "data_summary.txt", emit: summary

    script:
    """
    echo "=== GIAB HG002 Chr21 Data Summary ===" > data_summary.txt
    echo "Generated: \$(date)" >> data_summary.txt
    echo "" >> data_summary.txt

    echo "=== BAM Stats ===" >> data_summary.txt
    samtools flagstat ${bam} >> data_summary.txt
    echo "" >> data_summary.txt

    echo "=== FASTQ Stats ===" >> data_summary.txt
    echo "R1 file: ${r1}" >> data_summary.txt
    echo "R1 size: \$(ls -lh ${r1} | awk '{print \$5}')" >> data_summary.txt
    echo "R2 file: ${r2}" >> data_summary.txt
    echo "R2 size: \$(ls -lh ${r2} | awk '{print \$5}')" >> data_summary.txt
    echo "" >> data_summary.txt

    echo "=== Truth VCF Stats ===" >> data_summary.txt
    echo "VCF file: ${vcf}" >> data_summary.txt
    echo "VCF size: \$(ls -lh ${vcf} | awk '{print \$5}')" >> data_summary.txt
    echo "" >> data_summary.txt

    echo "=== High Confidence BED Stats ===" >> data_summary.txt
    echo "BED file: ${bed}" >> data_summary.txt
    echo "Regions: \$(wc -l < ${bed})" >> data_summary.txt
    echo "" >> data_summary.txt

    echo "=== Output Locations ===" >> data_summary.txt
    echo "Input FASTQs: ${params.outdir}/input/" >> data_summary.txt
    echo "Truth data: ${params.outdir}/truth/" >> data_summary.txt

    cat data_summary.txt
    """
}

process CREATE_SAMPLESHEET {
    publishDir "${params.outdir}", mode: 'copy'

    container 'ubuntu:22.04'

    input:
    path r1
    path r2

    output:
    path "samplesheet.csv", emit: csv

    script:
    // Use the published output paths for the samplesheet
    def r1_path = "${params.outdir}/input/${r1.name}"
    def r2_path = "${params.outdir}/input/${r2.name}"
    """
    echo "patient,sample,lane,fastq_1,fastq_2" > samplesheet.csv
    echo "HG002,HG002,1,${r1_path},${r2_path}" >> samplesheet.csv

    echo "Created samplesheet:"
    cat samplesheet.csv
    """
}

workflow {
    // Extract chr21 data in parallel
    EXTRACT_CHR21_BAM(params.chromosome)
    EXTRACT_CHR21_VCF(params.chromosome)
    EXTRACT_CHR21_BED(params.chromosome)

    // Convert BAM to FASTQ for sarek input
    BAM_TO_FASTQ(EXTRACT_CHR21_BAM.out.bam, EXTRACT_CHR21_BAM.out.bai)

    // Create samplesheet for sarek
    CREATE_SAMPLESHEET(BAM_TO_FASTQ.out.r1, BAM_TO_FASTQ.out.r2)

    // Collect statistics
    COLLECT_STATS(
        EXTRACT_CHR21_BAM.out.bam,
        EXTRACT_CHR21_BAM.out.bai,
        BAM_TO_FASTQ.out.r1,
        BAM_TO_FASTQ.out.r2,
        EXTRACT_CHR21_VCF.out.vcf,
        EXTRACT_CHR21_BED.out.bed
    )
}