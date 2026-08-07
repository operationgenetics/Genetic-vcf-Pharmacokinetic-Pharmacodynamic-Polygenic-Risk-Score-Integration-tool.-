import sqlite3

conn = sqlite3.connect('genomic_knowledgebase.db')
c = conn.cursor()

c.execute('DROP TABLE IF EXISTS disease_targeted_therapies')
c.execute('''
    CREATE TABLE disease_targeted_therapies (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        condition_trait TEXT,
        primary_drug TEXT,
        gene_checked TEXT,
        pgx_status TEXT,
        alternative_drug TEXT,
        clinical_rationale TEXT
    )
''')

enriched_therapies = [
    (
        'Anxiety Disorder',
        'Escitalopram',
        'CYP2C19',
        'CONTRAINDICATED',
        'Venlafaxine / Duloxetine',
        'Indicated for high polygenic risk of Anxiety. CYP2C19 Poor Metabolizer status causes reduced drug clearance and elevated serum concentration. CPIC recommends avoiding Escitalopram or switching to SNRI alternatives (Venlafaxine/Duloxetine).'
    ),
    (
        'Generalized Anxiety Disorder',
        'Escitalopram',
        'CYP2C19',
        'CONTRAINDICATED',
        'Buspirone / Duloxetine',
        'Indicated for elevated GAD polygenic burden. Impaired CYP2C19 metabolism significantly impairs primary SSRI clearance; re-routed to non-CYP2C19 dependent anxiolytics.'
    ),
    (
        'Panic Disorder',
        'Sertraline',
        'CYP2C19',
        'SUITABLE',
        'Clonazepam / SSRI Alternative',
        'First-line SSRI indicated for elevated Panic Disorder polygenic score. Normal CYP2C19 metabolic capacity ensures expected plasma clearance and therapeutic efficacy.'
    ),
    (
        'Post-Traumatic Stress Disorder',
        'Sertraline',
        'CYP2C19',
        'SUITABLE',
        'Prazosin / Venlafaxine',
        'First-line pharmacotherapy for high PTSD polygenic risk. Patient profile indicates normal hepatic metabolism, approving standard CPIC dosing protocols.'
    ),
    (
        'Schizoaffective Disorder',
        'Aripiprazole / Risperidone',
        'CYP2D6',
        'HIGH_RISK',
        'Clozapine / Olanzapine',
        'Indicated for high Schizoaffective polygenic risk. Patient is a CYP2D6 Poor Metabolizer (*10/*10), slowing elimination and increasing risk of extrapyramidal side effects. Recommends a 50% dose reduction or reassignment to Clozapine/Olanzapine.'
    ),
    (
        'Bipolar Disorder',
        'Carbamazepine',
        'HLA-B*15:02',
        'CONTRAINDICATED',
        'Valproate / Lamotrigine',
        'First-line mood stabilizer for elevated Bipolar polygenic risk. HLA-B*15:02 positivity carries high risk of severe cutaneous adverse reactions (SJS/TEN). Strictly contraindicated; reassigned to Valproate or Lamotrigine.'
    ),
    (
        'Major Depressive Disorder',
        'Escitalopram',
        'CYP2C19',
        'CONTRAINDICATED',
        'Sertraline / Mirtazapine',
        'Indicated for elevated MDD polygenic burden. CYP2C19 Poor Metabolizer profile inhibits clearance, increasing toxicity risk. CPIC guidelines advise switching to Sertraline or Mirtazapine.'
    ),
    (
        'Coronary Artery Disease',
        'Simvastatin',
        'SLCO1B1',
        'SUITABLE',
        'Rosuvastatin / Pravastatin',
        'Primary lipid-lowering therapy for elevated CAD risk. SLCO1B1 hepatic influx transporter function is normal (*1/*1), permitting standard simvastatin dosing without heightened myopathy risk.'
    ),
    (
        'Type 2 Diabetes',
        'Metformin',
        'SLC22A1',
        'SUITABLE',
        'GLP-1 Receptor Agonist',
        'First-line biguanide therapy for elevated T2D polygenic risk. SLC22A1 hepatic uptake transporter status is normal, ensuring standard therapeutic glycemic response.'
    )
]

c.executemany('''
    INSERT INTO disease_targeted_therapies 
    (condition_trait, primary_drug, gene_checked, pgx_status, alternative_drug, clinical_rationale)
    VALUES (?, ?, ?, ?, ?, ?)
''', enriched_therapies)

conn.commit()
conn.close()
print('[✔] Database successfully updated with comprehensive clinical rationales.')
