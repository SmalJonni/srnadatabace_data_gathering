SELECT
  pc.id, tax.id, rna.seq_short, pc.rna_type
FROM rnc_rna_precomputed AS pc
LEFT JOIN rnc_taxonomy AS tax
ON
  tax.id = pc.taxid
LEFT JOIN rna AS rna 
ON rna.upi = pc.upi
WHERE
  tax.lineage LIKE 'cellular organisms; Bacteria; %' -- grab only Bacteria
  AND pc.is_active = true    -- exclude sequences without active cross-references
  AND pc.rna_type = 'sRNA';
