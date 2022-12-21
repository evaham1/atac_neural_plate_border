// //ATAC channel = ch_atac
// ch_atac = [
//     [[meta:'HH5'], 'rds_atac'],
//     [[meta:'HH6'], 'rds_atac']
// ]

// //RNA channel = ch_rna
// ch_rna = [
//     [[meta:'HH5'], 'rds_rna'],
//     [[meta:'HH6'], 'rds_rna']

// ]

// Channel.from(ch_atac).set{ch_atac}
// Channel.from(ch_rna).set{ch_rna}

// // ch_atac.view()
// // ch_rna.view()

// // println 'channels now set'

// ch_atac
//     .concat(ch_rna)
//     .groupTuple(by:0)
//     .view()

// // // add full data to each channel
// // ch_atac
// //     .concat(ch_rna)
// //     .groupTuple(by:0)