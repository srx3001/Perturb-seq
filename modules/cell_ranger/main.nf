#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process run_cellranger_multi {
	publishDir 'results', mode: 'copy'
    input:
        val analysis_id
		val csv_file_name
        path rundir
		path reference

    output:
        path "${analysis_id}/outs/per_sample_outs/${analysis_id}/count/sample_filtered_feature_bc_matrix", emit: matrix
        path "${analysis_id}/outs/per_sample_outs/${analysis_id}/metrics.csv"
        path "${analysis_id}/outs/per_sample_outs/${analysis_id}/web_summary.html"

    script:
    """
    /cellranger-8.0.1/cellranger multi --id=${analysis_id} --csv=data/${csv_file_name}
    """
}

workflow {
    // Define the input directory and results output
    ch_id = Channel.value(params.analysis_id)
    ch_csv = Channel.value(params.csv_file_name)
	ch_rundir =  Channel.fromPath(params.rundir)
	ch_ref =  Channel.fromPath(params.reference)
 

    // Call the run_seurat process
    run_cellranger_multi(ch_id, ch_csv, ch_rundir, ch_ref)
}
