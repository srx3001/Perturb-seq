nextflow.enable.dsl=2

include { printHeader; helpMessage } from './help' params ( params )
include { run_cellranger_multi } from './modules/cell_ranger/main.nf' params ( params )
include { run_seurat } from './modules/seurat/main.nf' params ( params )



if ( params.help ) {
    helpMessage()
    exit 0
}


workflow {
       
    ch_id = Channel.value(params.analysis_id)
    ch_csv = Channel.value(params.csv_file_name)
	ch_rundir =  Channel.fromPath(params.rundir)
	ch_ref =  Channel.fromPath(params.reference)
    ch_seurat_out_dir = Channel.value(params.seurat_results)
	
    cell_ranger_outs = run_cellranger_multi(ch_id, ch_csv, ch_rundir, ch_ref)
	matrix_output = cell_ranger_outs.matrix
	run_seurat(matrix_output, ch_seurat_out_dir)
               
}

// OnComplete
workflow.onComplete{
    println( "\nPipeline completed successfully.\n\n" )
}

// OnError
workflow.onError{
    println( "\nPipeline failed.\n\n" )
}