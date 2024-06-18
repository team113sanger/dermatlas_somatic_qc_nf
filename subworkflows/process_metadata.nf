workflow DERMATLAS_METADATA {
    take: 
    pair_identities
    patient_metadata

    main:

    pair_identities 
    | splitCsv(sep:"\t", header:['tumor', 'normal']) 
    | map{ meta -> 
        [meta + [pair_id: meta.normal+ "_" + meta.tumor]]
        }
    | flatMap { meta -> 
    [
        [meta["normal"][0], meta],
        [meta["tumor"][0],  meta]
    ]}
    | set { pair_id_ch }
    
    patient_metadata
    | splitCsv(sep:"\t",header : true)
    | map {meta -> meta.subMap("Sex", "Sanger DNA ID", "OK_to_analyse_DNA?", "Phenotype")} 
    | set{ patient_metadata_ch }


    pair_id_ch.join(patient_metadata_ch).set{ combined_metadata }

    emit:
        combined_metadata

}