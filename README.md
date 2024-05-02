Outline of 'Downstream_analysis' pipeline:

PREPROCESSING (input: cellranger output directory)
    - EDIT_GTF (adds 'chr' to gtf)
    - PREPROCESS (sets up ArchR object from cellranger output: "/outs/fragments.tsv.gz")

FILTERING
    - FILTER (filters on nFrags)
    - SPLIT_STAGES
    - FILTER_CLUSTER_LOOP (filter clusters)
    - CLUSTER_POSTFILTER
            - GENE_SCORES_POSTFILTER
            - HEATMAP_GEX
    - PEAK_CALL_POSTFILTER
            - HEATMAP_PEAKS

FULL_PROCESSING
    - FILTER_FULL (use cell ids from stages to filter whole data)
    - CLUSTER_POSTFILTER
            - GENE_SCORES_POSTFILTER
            - HEATMAP_GEX
    - PEAK_CALL_POSTFILTER
            - HEATMAP_PEAKS
        
INTEGRATING
    - UNCON_INTEGRATE
            - CLUSTER_IDENTIFY 

TBC
