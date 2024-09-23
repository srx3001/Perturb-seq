#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process run_seurat {
	publishDir 'results', mode: 'copy'
    input:
        path data_dir
		val seurat_results

    output:
        path "${seurat_results}"

    script:
    """
    Seurat.R --data_dir ${data_dir} --output_dir ${seurat_results}
    """
}

workflow {
    // Define the input directory and results output
    ch_data_dir =  Channel.fromPath(params.data_dir)
    ch_out_dir = Channel.value(params.seurat_results)

    // Call the run_seurat process
    run_seurat(ch_data_dir, ch_out_dir)
}
