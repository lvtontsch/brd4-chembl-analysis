SELECT m.chembl_id AS "Molecule ChEMBL ID", cp.mw_freebase AS "Molecular Weight", cp.num_ro5_violations AS "#RO5 Violations", ac.pchembl_value AS "pChEMBL Value", ac.standard_value AS "Standard Value", a.description AS "Assay Description", bao.label AS "BAO Label", a.assay_cell_type AS "Assay Cell Type", a.variant_id AS "Assay Variant Mutation", a.doc_id AS "Document ID", s.src_description AS "Source Description", doc.journal AS "Document Journal", doc.year AS "Document Year"
FROM activities ac
JOIN assays a ON a.assay_id = ac.assay_id 
JOIN molecule_dictionary m ON m.molregno = ac.molregno
LEFT JOIN compound_properties cp ON cp.molregno = ac.molregno
LEFT JOIN bioassay_ontology bao ON bao.bao_id = a.bao_format
LEFT JOIN docs doc ON doc.doc_id = a.doc_id
JOIN source s ON s.src_id = a.src_id
WHERE a.tid = 103454 
AND ac.standard_type = 'IC50'
AND ac.standard_relation = '='
AND ac.standard_value IS NOT NULL
AND ac.standard_units = 'nM'
AND (ac.data_validity_comment IS NULL OR ac.data_validity_comment != 'Potential transcription error')
AND ac.potential_duplicate = 0
AND a.assay_type = 'B'
AND m.chembl_id != 'CHEMBL3220211'
AND s.src_description IN ('BindingDB Patent Bioactivity Data','Scientific Literature')
AND NOT (s.src_description = 'BindingDB Patent Bioactivity Data' AND ac.standard_value IN (50500,550,5500,2750));