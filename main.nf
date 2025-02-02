// main.nf: RNA-seq Pipeline using Nextflow

// Define parameters
params.reads = "data/raw/*_{R1,R2}.fastq.gz"
params.genome_index = "path/to/genome/index"  // Update with your genome index path
params.gtf = "path/to/annotation.gtf"        // Update with your GTF file path
params.outdir = "results"

// Print pipeline parameters
log.info """
    RNA-seq Pipeline
    ================
    Reads: ${params.reads}
    Genome Index: ${params.genome_index}
    GTF: ${params.gtf}
    Output Directory: ${params.outdir}
"""

// Define processes

// Process: Run FastQC for quality control
process fastqc {
    input:
    tuple val(sample), path(reads)

    output:
    path("${sample}_fastqc.html"), emit: fastqc_html

    script:
    """
    fastqc $reads -o .
    """
}

// Process: Trim reads using Trimmomatic
process trim_reads {
    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_R1_paired.fastq.gz"), path("${sample}_R2_paired.fastq.gz"), emit: trimmed_reads

    script:
    """
    trimmomatic PE -phred33 $reads \
        ${sample}_R1_paired.fastq.gz ${sample}_R1_unpaired.fastq.gz \
        ${sample}_R2_paired.fastq.gz ${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process: Align reads using STAR
process star_align {
    input:
    tuple val(sample), path(reads)

    output:
    path("${sample}_Aligned.sortedByCoord.out.bam"), emit: bam

    script:
    """
    STAR --genomeDir ${params.genome_index} \
        --readFilesIn $reads \
        --outFileNamePrefix ${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --runThreadN ${task.cpus}
    """
}

// Process: Count reads using featureCounts
process feature_counts {
    input:
    path bam
    path gtf

    output:
    path("counts.txt"), emit: counts

    script:
    """
    featureCounts -a $gtf -o counts.txt $bam
    """
}

// Process: Differential Expression Analysis using DESeq2 (R script)
process deseq2 {
    input:
    path counts
    path sample_info

    output:
    path("differential_expression_results.csv"), emit: de_results

    script:
    """
    Rscript scripts/deseq2.R $counts $sample_info differential_expression_results.csv
    """
}

// Workflow definition
workflow {
    // Read input files
    Channel.fromFilePairs(params.reads, flat: true)
        .set { read_pairs }

    // Run FastQC
    fastqc(read_pairs)
        .view()

    // Trim reads
    trim_reads(read_pairs)
        .view()

    // Align reads
    star_align(trim_reads.out.trimmed_reads)
        .view()

    // Count reads
    feature_counts(star_align.out.bam, params.gtf)
        .view()

    // Differential Expression Analysis
    deseq2(feature_counts.out.counts, file("data/sample_info.csv"))
        .view()
}

// Output results to the specified directory
workflow.onComplete {
    println "Pipeline completed successfully! Results are in ${params.outdir}"
}
