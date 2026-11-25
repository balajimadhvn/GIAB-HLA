#!/usr/bin/env nextflow
/*
 * GIAB-HLA: HLA class I and II typing from whole-genome short reads, benchmarked against published types.
 *
 *   reads (MHC + HLA contigs) ─┬─ T1K      (HLA-A, -B, -C, -DQB1, -DRB1 ...)
 *                              └─ OptiType (HLA-A, -B, -C)
 *                                     │
 *              2-field calls, concordance vs truth, tool agreement, trio consistency, R report
 */
nextflow.enable.dsl = 2

process T1K {
    tag "${meta.id}"
    label 't1k'
    publishDir "${params.outdir}/t1k", mode: 'copy', pattern: '*_genotype.tsv'

    input:
    tuple val(meta), path(r1), path(r2)
    path t1k_index

    output:
    tuple val(meta), path("${meta.id}_genotype.tsv"), emit: calls

    script:
    """
    idx=${t1k_index}/t1k_hla_dna_seq.fa
    # the bundled test index is gzipped to keep the repo small
    if [ ! -s \$idx ]; then gunzip -c \$idx.gz > t1k_hla_dna_seq.fa; idx=t1k_hla_dna_seq.fa; fi
    run-t1k -1 ${r1} -2 ${r2} --preset hla -f \$idx \\
        -t ${task.cpus} -o ${meta.id} > t1k.log 2>&1
    """
}

process OPTITYPE {
    tag "${meta.id}"
    label 'optitype'
    publishDir "${params.outdir}/optitype", mode: 'copy', pattern: '*_result.tsv'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_result.tsv"), emit: calls

    script:
    """
    # OptiType reads its settings from config.ini next to the script; write a local one
    bindir=\$(dirname \$(command -v OptiTypePipeline.py))
    cat > config.ini <<CFG
[mapping]
razers3=\$(command -v razers3)
threads=${task.cpus}
[ilp]
solver=glpk
threads=1
[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
CFG
    # OptiType's shebang is "env python"; call its own environment's interpreter explicitly
    \$bindir/python \$bindir/OptiTypePipeline.py -i ${r1} ${r2} --dna -c config.ini -o out -p ${meta.id} > optitype.log 2>&1
    cp out/${meta.id}_result.tsv .
    """
}

process SUMMARISE {
    label 'base'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path t1k_calls
    path optitype_calls
    path samples
    path truth

    output:
    path "tables", emit: tables

    script:
    """
    hla_concordance.py --t1k ${t1k_calls} --optitype ${optitype_calls} \\
        --samples ${samples} --truth ${truth} --outdir tables
    """
}

process REPORT {
    label 'base'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tables
    path rmd

    output:
    path "hla_report.html"
    path "figures"

    script:
    """
    cp -L ${rmd} report.Rmd
    Rscript -e "rmarkdown::render('report.Rmd', output_file = 'hla_report.html', params = list(tables = '\$PWD/${tables}'), quiet = TRUE)"
    """
}

workflow {
    // relative FASTQ paths are read relative to the samplesheet
    def sheet = file(params.input, checkIfExists: true)
    def fq = { String p -> p.startsWith('/') ? file(p, checkIfExists: true) : file("${sheet.parent}/${p}", checkIfExists: true) }
    reads = Channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { r -> [ [id: r.sample], fq(r.fastq_1), fq(r.fastq_2) ] }

    T1K(reads, file(params.t1k_index, checkIfExists: true))
    OPTITYPE(reads)

    SUMMARISE(
        T1K.out.calls.map { it[1] }.collect(),
        OPTITYPE.out.calls.map { it[1] }.collect(),
        file(params.samples, checkIfExists: true),
        file(params.truth, checkIfExists: true)
    )
    REPORT(SUMMARISE.out.tables, file("${projectDir}/reports/hla_report.Rmd"))
}
