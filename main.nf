nextflow.enable.dsl=2

include {QC} from './modules/qc.nf'
include {NORMALIZE} from './modules/normalize.nf'
include {CLUSTER} from './modules/cluster.nf'

workflow{
    raw_data_ch = channel.fromPath("${params.data_path}", type: 'dir')
    QC(raw_data_ch)
    NORMALIZE(QC.out[0])
    CLUSTER(NORMALIZE.out[0])
}