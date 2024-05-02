//LOOP_CALL.out after:
//.map { row -> [ row[0], [ row[1].findAll { it =~ ".*rds_files" }[0] ] ] }
// loop_call_out = [
//     [[sample_id:'WE_HiChip_r1'], 'work/39/39b12/rds_files'],
//     [[sample_id:'WE_HiChip_r2'], 'work/48/93b63/rds_files']
// ]
loop_call_out = [
    [[sample_id:'WE_HiChip_r2'], 'work/48/93b63/rds_files']
]

//INTERSECT_BINS_PEAKS.out:
intersect_peaks_out = ['work/33/12gh56/FullData_PeakSet_bins_intersected.bed']

//INTERSECT_BINS_GENES.out:
intersect_genes_out = ['work/5a/9a9ed/tag_chroms_bins_intersected.bed']


// Make channels and check them
Channel.from(loop_call_out).set{loop_call_out}
Channel.from(intersect_peaks_out).set{intersect_peaks_out}
Channel.from(intersect_genes_out).set{intersect_genes_out}

// loop_call_out.view()
// intersect_peaks_out.view()
// intersect_genes_out.view()

// Try combining channels
//Channel.combine(loop_call_out, intersect_peaks_out, intersect_genes_out).view()

loop_call_out
    .combine(intersect_peaks_out)
    .combine(intersect_genes_out)
    .map{ [ it[0], [it[1], it[2], it[3]] ] }
    .view()
//[[sample_id:WE_HiChip_r2], [work/48/93b63/rds_files, work/33/12gh56/FullData_PeakSet_bins_intersected.bed, work/5a/9a9ed/tag_chroms_bins_intersected.bed]]

// println 'channels now set'

// ch_atac
//     .concat(ch_rna)
//     .groupTuple(by:0)
//     .view()

// // add full data to each channel
// ch_atac
//     .concat(ch_rna)
//     .groupTuple(by:0)