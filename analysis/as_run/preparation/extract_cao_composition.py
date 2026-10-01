from pathlib import Path
import csv,openpyxl
root=Path('<workspace>/reanalysis_2026-09-08')
book=openpyxl.load_workbook(root/'inputs/Cao2025_source_data.xlsx',read_only=True,data_only=True)
ws=book['Supplementary Fig. 2']
clinical=openpyxl.load_workbook(root/'inputs/Cao2025_41467_2024_55419_MOESM4_ESM.xlsx',read_only=True,data_only=True)
clinical_rows=list(clinical.active.values)
donors={r[1]:r[0] for r in clinical_rows[1:] if r[0]}
mapping=list(csv.DictReader(Path('<workspace>/sdrive_audit_20260908/metadata_audit/GSE250130_barcode_mapping.csv').open()))
print('Mapping columns',list(mapping[0]))
rows=[]
for col in range(3,31):
    sample=ws.cell(16,col).value
    row={'source_sample':sample,'donor':donors[sample],'condition':'RIF' if sample.startswith('RIF') else 'Control','group':sample.rsplit('_',1)[0]}
    for r in range(17,25): row[ws.cell(r,2).value]=ws.cell(r,col).value
    assert abs(sum(row[k] for k in ['Stromal','Unciliated','Ciliated','Endothelial','NK/T','Myeloid','B','Mast'])-1)<1e-8
    rows.append(row)
with (root/'inputs/Cao2025_donor_composition.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('Extracted',len(rows),'donor compositions; source workbook Supplementary Fig. 2, C16:AD24')
