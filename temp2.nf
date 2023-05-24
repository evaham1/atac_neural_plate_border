//LOOP_CALL.out after:
//.map{it[1].findAll{it =~ /rds_files/}[0].listFiles()}
// loop_call_out = [
//     ['WE_HiChip_r1_HiCDC_output_filtered.txt', 'WE_HiChip_r1_HiCDC_output.txt.gz'],
//     ['NF_HiChip_r1_HiCDC_output.txt.gz', 'NF_HiChip_r1_HiCDC_output_filtered.txt'],
//     ['WE_HiChip_r2_HiCDC_output_filtered.txt', 'WE_HiChip_r2_HiCDC_output.txt.gz'],
//     ['WE_HiChip_r3_HiCDC_output.txt.gz', 'WE_HiChip_r3_HiCDC_output_filtered.txt'],
//     ['NF_HiChip_r2_HiCDC_output_filtered.txt', 'NF_HiChip_r2_HiCDC_output.txt.gz'],
//     ['NF_HiChip_r3_HiCDC_output.txt.gz', 'NF_HiChip_r3_HiCDC_output_filtered.txt']
// ]

loop_call_out = [
    ['[WE_HiChip_r1_HiCDC_output_filtered.txt', 'WE_HiChip_r1_HiCDC_output.txt.gz'], ['WE_HiChip_r2_HiCDC_output_filtered.txt', 'WE_HiChip_r2_HiCDC_output.txt.gz'], ['NF_HiChip_r1_HiCDC_output.txt.gz', 'NF_HiChip_r1_HiCDC_output_filtered.txt'], ['WE_HiChip_r3_HiCDC_output.txt.gz', 'WE_HiChip_r3_HiCDC_output_filtered.txt'], ['NF_HiChip_r2_HiCDC_output_filtered.txt', 'NF_HiChip_r2_HiCDC_output.txt.gz'], ['NF_HiChip_r3_HiCDC_output.txt.gz', 'NF_HiChip_r3_HiCDC_output_filtered.txt']]

// Make channels and check them
Channel.from(loop_call_out).set{loop_call_out}
//loop_call_out.view()


loop_call_out
    .flatMap().collect()
    .view()
    .map { [[sample_id:'AllSamples'], it] } //
    .view()





