//stages channel = ch_stages
ch_stages = [
    [[meta:'HH5'], 'HH5_RDS'],
    [[meta:'HH6'], 'HH6_RDS']
]

//RNA channel = ch_rna
ch_full = [
    [[meta:'FullData'], 'Full_RDS'],
]

// Make channels and check them
Channel.from(ch_stages).set{ch_stages}
Channel.from(ch_full).set{ch_full}

// ch_stages.view()
// ch_full.view()

// Want to combine so each stage also has a full data
ch_stages
    .combine(ch_full)
    .view() //[[meta:HH5], HH5_RDS, [meta:FullData], Full_RDS]
    .map{ [ it[0], [ it[1], it[3] ] ] }
    .view() //[[meta:HH5], [HH5_RDS, Full_RDS]]







