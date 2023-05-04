process GTF_TO_BED {

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path gtf_file

    output:
    path "$output"       , emit: bed

    script:
    output = gtf_file.toString() - ".gtf" + ".bed"
    """
    convert_gtf_to_bed.sh $gtf_file $output
    """
}

    // sed '/^#/d' "./input/tag_chroms.gtf" \
    //     | awk -F "\t" '{
    //         if ($3 == "gene") {
    //             split($9, a, "\"");
    //             gene_id = a[2];
    //             chr = "chr" $1;
    //             gene_name = (a[10] == "") ? gene_id : a[6];
    //             print chr "\t" $4-1 "\t" $5 "\t" gene_id "\t" gene_name "\t" $7;
    //         }
    //     }' \
    //     > "gtf_file.bed"

    
    // sed '/^#/d' ./input/* \
    // | awk -F "\t" '{
    //     if ($3=="gene") {
    //         split($9, a, "\"");
    //         gene_id=a[2];
    //         chr="chr"$1
    //         if (a[10]=="") { 
    //             gene_name=gene_id;
    //         } else { 
    //             gene_name=a[6];
    //         } 
    //     print chr"\t"$4-1"\t"$5"\t"gene_id"\t"gene_name"\t"$7;
    // }
    // }' \
    // > "gtf_file.bed"