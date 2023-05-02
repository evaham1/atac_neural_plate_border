process EDIT_VALIDPAIRS {
    input:
    tuple val(meta), path('input/*')

    output:
    tuple val(meta), file('*')       , emit: outs

    script:
    def filename = sh(script: 'basename ./input/*', returnStdout: true).trim()

    def filename = sh(script: 'basename /path/to/directory/*', returnStdout: true).trim()


    """
    basename ./input/*
    cp ./input/WE_HiChip_r1.allValidPairs input.txt
    """
    // # Detect input file name
    // filename=$(basename ./input/*)

    // # Rename input file to .txt
    // cp $filename input.txt

    // # Split the input file into smaller chunks
    // split -l 1000000 input.txt input_part

    // # Edit the second and fifth column of each chunk to add 'chr' to the chromosome name
    // for file in input_part*; do
    //     awk 'BEGIN{FS=OFS="\t"}{$2="chr"$2;$5="chr"$5}1' "${file}" > "${file}.edited" &
    // done

    // # Wait for all editing jobs to finish
    // wait

    // # Concatenate the edited chunks into a single output file
    // cat input_part*.edited > "edited_ValidPairs.txt"

    // # Remove the intermediate files
    // rm input.txt input_part*

    // # Convert the output file to tab-delimited format
    // sed -i 's/ /\t/g' "edited_ValidPairs.ValidPairs"
    //"""
}
