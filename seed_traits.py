import sqlite3

conn = sqlite3.connect('genomics.db')
cursor = conn.cursor()

cursor.execute('''
CREATE TABLE IF NOT EXISTS prs_traits (
    trait_id INTEGER PRIMARY KEY AUTOINCREMENT,
    trait_name TEXT UNIQUE,
    mean_raw_score REAL,
    std_raw_score REAL
)''')

cursor.execute('''
CREATE TABLE IF NOT EXISTS prs_variants (
    variant_id INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid TEXT,
    trait_name TEXT,
    effect_allele TEXT,
    weight REAL
)''')

traits = [
    ('Anxiety Disorder', 0.2, 0.4),
    ('Generalized Anxiety Disorder', 0.2, 0.4),
    ('Panic Disorder', 0.2, 0.4),
    ('Post-Traumatic Stress Disorder', 0.2, 0.4),
    ('Schizoaffective Disorder', 0.2, 0.4),
    ('Bipolar Disorder', 0.2, 0.4),
    ('Major Depressive Disorder', 0.2, 0.4),
    ('Coronary Artery Disease', 0.5, 0.3),
    ('Type 2 Diabetes', 0.5, 0.3)
]

for name, mean, std in traits:
    cursor.execute('''
        INSERT INTO prs_traits (trait_name, mean_raw_score, std_raw_score)
        VALUES (?, ?, ?)
        ON CONFLICT(trait_name) DO UPDATE SET 
            mean_raw_score=excluded.mean_raw_score, 
            std_raw_score=excluded.std_raw_score
    ''', (name, mean, std))

variants = [
    ('rs11178997', 'Anxiety Disorder', 'G', 0.85),
    ('rs28399433', 'Generalized Anxiety Disorder', 'T', 0.85),
    ('rs1799853',  'Panic Disorder', 'T', 0.85),
    ('rs4244285',  'Post-Traumatic Stress Disorder', 'A', 0.85),
    ('rs1024611',  'Schizoaffective Disorder', 'A', 0.85),
    ('rs9272219',  'Bipolar Disorder', 'T', 0.85),
    ('rs1065852',  'Major Depressive Disorder', 'T', 0.85)
]

for rsid, trait, allele, weight in variants:
    cursor.execute('''
        INSERT OR REPLACE INTO prs_variants (rsid, trait_name, effect_allele, weight)
        VALUES (?, ?, ?, ?)
    ''', (rsid, trait, allele, weight))

conn.commit()
conn.close()
print('[✔] genomics.db updated successfully with all 9 traits.')
