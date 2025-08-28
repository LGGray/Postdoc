import pandas as pd

# Read coordinates
coords = pd.read_csv('non_overlapping_coords.txt', sep=',')

# Adjust flanks
coords['Gene start (bp)'] = coords.apply(
    lambda row: max(0, row['Gene start (bp)'] - 5000) if row['Strand'] == '1' else row['Gene start (bp)'],
    axis=1
)
coords['Gene end (bp)'] = coords.apply(
    lambda row: row['Gene end (bp)'] + 5000 if row['Strand'] == '-1' else row['Gene end (bp)'],
    axis=1
)

# Rename and prepare BED fields
coords.rename(columns={
    'Chromosome/scaffold name': 'chrom',
    'Gene start (bp)': 'start',
    'Gene end (bp)': 'end',
    'Strand': 'strand',
    'Gene name': 'name'
}, inplace=True)

# Convert Ensembl strand to BED strand
coords['strand'] = coords['strand'].replace({'1': '+', '-1': '-'})

# Add a dummy score column for BED format
coords['score'] = 4

# Redorder bed format
coords = coords[['chrom', 'start', 'end', 'name', 'strand']]

coords['start'] = coords['start'].astype(int)
coords['end'] = coords['end'].astype(int)

# Save as BED
coords.to_csv('lncRNA_5kb_flank.bed', sep='\t', header=False, index=False)