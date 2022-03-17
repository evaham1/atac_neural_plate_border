#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import groovy libs
import groovy.transform.Synchronized

workflow METADATA {
    take: file_path
    main:
        Channel
            .fromPath( file_path )
            .splitCsv(header:true)
            .map { row -> ProcessRow(row) }
            .set { metadata }
    emit:
        metadata
}

def ProcessRow(LinkedHashMap row, boolean flattenData = false) {
    def meta = [:]
    meta.sample_id = row.sample_id

    for (Map.Entry<String, ArrayList<String>> entry : row.entrySet()) {
        String key = entry.getKey();
        String value = entry.getValue();
    
        if(key != "sample_id" && key != "data") {
            meta.put(key, value)
        }
    }

    def array = []

    data = file(row.data, checkIfExists: true)

    if (data instanceof List) {
        array = [ meta, data ] // read files from glob list
    } else if (data instanceof Path){
        array = [ meta, [ data ] ] //read path
    } else {
        throw new Exception("data class not recognised")
    }

    return array
}