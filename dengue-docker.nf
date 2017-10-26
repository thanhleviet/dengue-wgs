#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017 Thanh Le Viet

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// General configuration variables
params.help = false
params.pwd = "$PWD"
params.input = "/Users/thanhlv/Downloads/offline/WGS_Dengue/HCM/raw/tmp"
params.output = "den-hcm"
params.reads = "*_L001_R{1,2}_001.fastq.gz"
params.read_pairs = params.input + "/" + params.reads
params.out_dir = "$baseDir/" + params.output
params.refdb = "/Users/thanhlv/Downloads/offline/WGS_Dengue/Ref/"
params.primerdb = "$baseDir/" + "primer"
params.threads = 8
params.trim = true
params.normalize = false
params.min_dp = 10
params.min_baseq = 13

threads = params.threads
mf = 1
// Trimmomatic configuration variables
params.leading = 20
params.trailing = 20
params.slidingwindow = "3:30"
params.minlen = 30
params.adapters = "NexteraPE-PE.fa"

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen
adapters = params.adapters

// SCRIPTS = "$baseDir/scripts"
LB = "Nextera"

// Reference
ref_db = file(params.refdb)
ref_denv = file(ref_db + "/denv.fasta")
ref_denv1 = "DENV1_FJ687432.1.fasta"
ref_denv2 = "DENV2_FJ639718.1.fasta"
ref_denv3 = "DENV3_AB189127.1.fasta"
ref_denv4 = "DENV4_AY618993.1.fasta"

// Primer
primerdb = file(params.primerdb)
primer_denv1 = "primer_d1.fasta"
primer_denv2 = "primer_d2.fasta"
primer_denv3 = "primer_d3.fasta"
primer_denv4 = "primer_d4.fasta"


// Display help menu
if(params.help) {
    log.info ''
    log.info 'OUCRU - Dengue NGS Pipelines'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow denue-docker.nf -profile docker [options]'
    log.info ''
    log.info 'General Options: '
    log.info '    --input           DIR     Directory of paired FASTQ files'
    log.info '    --reads           STRING["*_L001_R{1,2}_001.fastq.gz"]  Pattern for recognizing paired reads'
    log.info '    --threads         INT     Number of threads to use for each process'
    log.info '    --output          DIR     Directory to write output files to'
    log.info '    --refdb           DIR     Directory of 4 reference fasta files (DENV1, DENV2, DENV3, DENV4)'
    log.info '    --primerdb        DIR     Directory of 4 primer fasta files (DENV1, DENV2, DENV3, DENV4)'
    log.info ''
    log.info 'Trimmomatic Options: '
    log.info '    --leading         INT[20]     Remove leading low quality or N bases'
    log.info '    --trailing        INT[20]     Remove trailing low quality or N bases'
    log.info '    --slidingwindow   INT[3:30]   Scan read with a sliding window'
    log.info '    --minlen          INT[30]     Drop reads below INT bases long'
    log.info ''
    log.info 'Notes: [defaul value]'
    log.info 'thanhlv@oucru.org'
    log.info ''
    return
}


Channel
    .fromFilePairs(params.read_pairs, flat: true)
    .ifEmpty { exit 1, "Read pairs could not be found: ${params.read_pairs}" }
    .set { trimmomatic_read_pairs }



/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process RunQC {
    publishDir "${params.out_dir}/${dataset_id}", mode: "copy"

    tag { dataset_id }

    maxForks mf

    input:
    set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

    output:
    set dataset_id, file("${dataset_id}_R1.fastq.gz"), file("${dataset_id}_R2.fastq.gz"), file("stats") into QC
    file("trimmomatic.config") into ch_trimmomatic_log

    """
    mkdir stats

    trimmomatic PE -threads ${threads} \
    $forward $reverse -baseout ${dataset_id} \
    ILLUMINACLIP:${TRIMMOMATIC}/${adapters}:2:30:10 \
    LEADING:${leading} \
    TRAILING:${trailing} \
    SLIDINGWINDOW:${slidingwindow} \
    MINLEN:${minlen} 2>&1 | tee ${dataset_id}_trimmomatic.log

    cat ${dataset_id}_1P | gzip -1 > ${dataset_id}_R1.fastq.gz
    cat ${dataset_id}_2P | gzip -1 > ${dataset_id}_R2.fastq.gz
    rm ${dataset_id}_1P
    rm ${dataset_id}_2P

    echo \"${adapters}:2:30:10
    LEADING:${leading}
    TRAILING:${trailing}
    SLIDINGWINDOW:${slidingwindow}
    MINLEN:${minlen}\" > trimmomatic.config

    fastqc $forward $reverse *_R*.fastq.gz -o stats

    mv ${dataset_id}_trimmomatic.log stats/${dataset_id}_trimmomatic.log
    """
    }

process BWA_mapping {
    publishDir "${params.out_dir}/${dataset_id}", mode: "copy", pattern: "*.{bam,bai,fasta,fai,dict,tsv,bed,pos}"

    tag {dataset_id}

    input:
    set dataset_id, file(forward), file(reverse), file("stats") from QC
    file(primerdb)
    file(ref_db)
    output:
    set dataset_id, file("${dataset_id}.bam"), file("${dataset_id}.bam.bai"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") into ch_bwa_mark_duplicate

    shell:
    '''
    bwa mem -t !{threads} !{ref_db}/denv.fasta \
     !{forward} \
     !{reverse} | samtools view -@ !{threads} -Sb - | samtools sort  -@ !{threads} -o tmp.bam -

    samtools index -@ !{threads} tmp.bam

    st=`samtools idxstats tmp.bam | sort -k3 -nr | awk 'NR==1 {print $1}' | tr -d '\n'`

    # Query the correct reference genome and primer
    !{SCRIPTS}/get_ref.py $st \
    !{ref_db}/!{ref_denv1} \
    !{ref_db}/!{ref_denv2} \
    !{ref_db}/!{ref_denv3} \
    !{ref_db}/!{ref_denv4} \
    !{primerdb}/!{primer_denv1} \
    !{primerdb}/!{primer_denv2} \
    !{primerdb}/!{primer_denv3} \
    !{primerdb}/!{primer_denv4}

    # Map to correct reference
    bwa mem -t !{threads} -M -R "@RG\\tID:!{dataset_id}\\tSM:!{dataset_id}\\tLB:!{LB}\\tPL:ILLUMINA" \
    new_ref.fasta \
    !{forward} \
    !{reverse} | samtools view -@ !{threads} -Sb - | samtools sort  -@ !{threads} -o !{dataset_id}.bam -

    samtools index -@ !{threads} !{dataset_id}.bam

    !{SCRIPTS}/primer_bed.py -p primer.fasta -r new_ref.fasta -o primer

    rm tmp.bam
    '''

}


process Mark_duplicate {

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}.bam"), file("${dataset_id}.bam.bai"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") from ch_bwa_mark_duplicate

    output:
    set dataset_id, file("${dataset_id}_dedup.bam"), file("${dataset_id}_dedup.bam.bai"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") into (ch_dedup_align_metrics, ch_dedup_PreLofreq)

    """
    ${SCRIPTS}/mark_primer.py -i ${dataset_id}.bam -o ${dataset_id}_dedup.bam -p primer.pos
    samtools index -@ ${threads} ${dataset_id}_dedup.bam
    """
}


process AligmentMetrics {
    publishDir "${params.out_dir}/${dataset_id}/stats", mode: "copy" //, pattern: "*.{txt,pdf}"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}_dedup.bam"), file("${dataset_id}_dedup.bam.bai"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict") from ch_dedup_align_metrics.take( 6 )

    output:
    file("${dataset_id}_alignment_metrics.txt")
    file("${dataset_id}_insert_metrics.txt")
    file("${dataset_id}_insert_size_histogram.pdf")

    """
    picard CollectAlignmentSummaryMetrics \
    R=new_ref.fasta \
    I=${dataset_id}_dedup.bam \
    O=${dataset_id}_alignment_metrics.txt

    picard CollectInsertSizeMetrics \
    INPUT=${dataset_id}_dedup.bam \
    OUTPUT=${dataset_id}_insert_metrics.txt \
    HISTOGRAM_FILE=${dataset_id}_insert_size_histogram.pdf
    """
}

process Realign_PreLofreq {
    publishDir "${params.out_dir}/${dataset_id}", mode: "copy", pattern: "*.{bam,bai,vcf}"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}_dedup.bam"), file("${dataset_id}_dedup.bam.bai"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") from ch_dedup_PreLofreq

    output:
    set dataset_id, file("${dataset_id}_realigned_reads.bam"), file("${dataset_id}_realigned_reads.bai"), file("${dataset_id}_pre_lofreq.vcf"), file("${dataset_id}_pre_lofreq_raw_snps.vcf"), file("${dataset_id}_pre_lofreq_raw_indels.vcf"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") into ch_PreLofreq_Select_SNPs_INDEL

    """
    # Indel Realignment
    gatk -nt $threads -T RealignerTargetCreator -R new_ref.fasta \
    -I ${dataset_id}_dedup.bam \
    -o realignment_targets.list

    gatk -T IndelRealigner -R new_ref.fasta \
    -I ${dataset_id}_dedup.bam \
    -targetIntervals realignment_targets.list \
    -o ${dataset_id}_realigned_reads.bam

    samtools index -@ ${threads} ${dataset_id}_realigned_reads.bam

    # Pre process BAM file
    lofreq viterbi -f new_ref.fasta ${dataset_id}_realigned_reads.bam | \
    lofreq alnqual -u - new_ref.fasta | \
    lofreq indelqual --dindel -f new_ref.fasta - | \
    samtools sort -@ $threads -o viterbi.bam -
    samtools index viterbi.bam

    # Variants calling
    lofreq call-parallel --call-indels --pp-threads $threads -l primer.bed -f new_ref.fasta -o ${dataset_id}_pre_lofreq.vcf viterbi.bam

    # Extract snps
    gatk -T SelectVariants \
    -R new_ref.fasta \
    -V ${dataset_id}_pre_lofreq.vcf \
    -selectType SNP \
    -o ${dataset_id}_pre_lofreq_raw_snps.vcf

    # Extract INDELS
    gatk -T SelectVariants \
    -R new_ref.fasta \
    -V ${dataset_id}_pre_lofreq.vcf \
    -selectType INDEL \
    -o ${dataset_id}_pre_lofreq_raw_indels.vcf
    """
}

process BQSR {
    publishDir "${params.out_dir}/${dataset_id}", mode: "copy", pattern: "*.{bam,bai,pdf}"
    publishDir "${params.out_dir}/${dataset_id}/stats", mode: "copy", pattern: "recal_data.table"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}_realigned_reads.bam"), file("${dataset_id}_realigned_reads.bai"), file("${dataset_id}_pre_lofreq.vcf"), file("${dataset_id}_pre_lofreq_raw_snps.vcf"), file("${dataset_id}_pre_lofreq_raw_indels.vcf"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") from ch_PreLofreq_Select_SNPs_INDEL

    output:
    set dataset_id, file("${dataset_id}_realigned_reads.bam"), file("${dataset_id}_realigned_reads.bai"), file("recal_data.table"), file("post_recal_data.table"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") into ch_bqsr_result
    file("${dataset_id}_recalibration_plots.pdf")
    """
    # Round 1
    gatk -T BaseRecalibrator \
    -R new_ref.fasta \
    -I ${dataset_id}_realigned_reads.bam \
    -knownSites ${dataset_id}_pre_lofreq_raw_snps.vcf \
    -knownSites ${dataset_id}_pre_lofreq_raw_indels.vcf \
    -o recal_data.table

    # Round 2
    gatk -T BaseRecalibrator \
    -R new_ref.fasta \
    -I ${dataset_id}_realigned_reads.bam \
    -knownSites ${dataset_id}_pre_lofreq_raw_snps.vcf \
    -knownSites ${dataset_id}_pre_lofreq_raw_indels.vcf \
    -BQSR recal_data.table \
    -o post_recal_data.table

    # AnalyseVariant
    gatk -T AnalyzeCovariates \
    -R new_ref.fasta \
    -before recal_data.table \
    -after post_recal_data.table \
    -plots ${dataset_id}_recalibration_plots.pdf
    """
}

process ApplyBQSR {
    publishDir "${params.out_dir}/${dataset_id}", mode: "copy", pattern: "*.{bam,bai}"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}_realigned_reads.bam"), file("${dataset_id}_realigned_reads.bai"), file("recal_data.table"), file("post_recal_data.table"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") from ch_bqsr_result

    output:
    set dataset_id, file("${dataset_id}_recal_reads.bam"), file("${dataset_id}_recal_reads.bam.bai"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") into ch_apply_bqrs_result

    """
    gatk -T PrintReads \
    -R new_ref.fasta \
    -I ${dataset_id}_realigned_reads.bam \
    -BQSR recal_data.table \
    -o ${dataset_id}_recal_reads.bam

    samtools index -@ $threads ${dataset_id}_recal_reads.bam
    """
}

process PostLofreq {
    publishDir "${params.out_dir}/${dataset_id}", mode: "copy", pattern: "*.{vcf,bedgraph,bed,pos,tsv,html}"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}_recal_reads.bam"), file("${dataset_id}_recal_reads.bam.bai"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("primer.pos"), file("primer.bed"), file("primer.tsv"), file("stats") from ch_apply_bqrs_result

    output:
    set dataset_id, file("${dataset_id}_recal_reads.bam"), file("${dataset_id}_post_lofreq.vcf"), file("${dataset_id}_post_lofreq_final.vcf"), file("exc_dup_${dataset_id}_post_lofreq_final.vcf"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("${dataset_id}_mask.vcf"), file("stats") into ch_post_lofreq
    file("${dataset_id}_recal_genomecov.bedgraph")
    file("${dataset_id}_recal_zero_genomecov.bed")
    """
    lofreq call-parallel \
    --call-indels \
    --pp-threads $threads \
    -f new_ref.fasta \
    -q ${params.min_baseq} \
    -Q ${params.min_baseq} \
    -s \
    -l primer.bed \
    -o ${dataset_id}_post_lofreq.vcf ${dataset_id}_recal_reads.bam

    bedtools genomecov -bga -ibam ${dataset_id}_recal_reads.bam > ${dataset_id}_recal_genomecov.bedgraph

    bedtools genomecov -d -ibam ${dataset_id}_recal_reads.bam > ${dataset_id}_recal_zero_genomecov.bed

    ${SCRIPTS}/bed2vcf.py ${dataset_id}_recal_zero_genomecov.bed new_ref.fasta ${dataset_id} ${params.min_dp}

    ${SCRIPTS}/vcf_filter.py -i ${dataset_id}_post_lofreq.vcf -o ${dataset_id}_post_lofreq_final.vcf
    """
}


process Consensus {
    publishDir "${params.out_dir}/${dataset_id}/consensus", mode: "copy", pattern: "*.{fasta}"

    tag {dataset_id}

    input:
    set dataset_id, file("${dataset_id}_recal_reads.bam"), file("${dataset_id}_post_lofreq.vcf"), file("${dataset_id}_post_lofreq_final.vcf"), file("exc_dup_${dataset_id}_post_lofreq_final.vcf"), file("new_ref.fasta"), file("new_ref.fasta.fai"), file("new_ref.dict"), file("${dataset_id}_mask.vcf"), file("stats") from ch_post_lofreq

    output:
    file("*.fasta")
    set val(dataset_id), file("stats"), file("${dataset_id}_recal_reads.bam"), file("${dataset_id}_post_lofreq_final.vcf") into ch_stats_report

    """
    gatk -T FastaAlternateReferenceMaker \
    -R new_ref.fasta \
    -o ${dataset_id}_consensus\
    -V ${dataset_id}_post_lofreq_final.vcf\
    -snpmask ${dataset_id}_mask.vcf \
    --snpmaskPriority

    ${SCRIPTS}/fix_consensus.py ${dataset_id}_consensus new_ref.fasta ${dataset_id}
    """
}

process Report {
    publishDir "${params.out_dir}/${dataset_id}", mode: "copy"

    tag {dataset_id}

    input:
    set val(dataset_id), file("stats"), file("${dataset_id}_recal_reads.bam"), file("${dataset_id}_post_lofreq_final.vcf") from ch_stats_report

    output:
    file("stats")
    file("report")

    """
    bcftools stats ${dataset_id}_post_lofreq_final.vcf > stats/bcftools.stats
    samtools stats ${dataset_id}_recal_reads.bam > stats/samtools.stats

    echo \"id: 'file-list'
    section_name: 'Output files'
    plot_type: 'html'
    data: |
        <table>
        <thead class=\"rowheader header\">
        <th>File</th>
        <th>Description</th>
        </thead>
        <tbody>
        <tr><th class=\"rowheader\">${dataset_id}.bam</th><td>Raw BAM file generated by BWA MEM</td></tr>
        <tr><th class=\"rowheader\">${dataset_id}_recal_reads.bam</th><td>BAM file after PCR mark duplicate</td></tr>
        <tr><th class=\"rowheader\">${dataset_id}_pre_lofreq.vcf</th><td>Raw VCF file prepared for Base Quality Recalibration</td></tr>
        <tr><th class=\"rowheader\">${dataset_id}_pre_lofreq_final.vcf</th><td>Final VCF file used for calling consensus fasta</td></tr>
        <tr><th class=\"rowheader\">${dataset_id}_mask.vcf</th><td>VCF file with locus Depth of Coverage < <b>${params.min_dp}</b>, used for replacing with N character when calling consensus fasta</td></tr>
        <tr><th class=\"rowheader\">consensus</th><td>Folder contains consensus fasta file</td></tr>
        </tbody>
        </table>
    \" > stats/file_list_mqc.yaml

    cp $SCRIPTS/workflow_mqc.yaml stats/
    multiqc -f -c $SCRIPTS/multiqc_config.yaml -o report -n report_${dataset_id} -d stats -dd -1
    """
}

workflow.onComplete {
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command Line:     $workflow.commandLine"
    log.info "Container:        $workflow.container"
    log.info "Duration:     $workflow.duration"
    log.info "Output Directory: $params.out_dir"
}
