import sqlite3

conn = sqlite3.connect('genomics.db')
c = conn.cursor()

c.execute('DROP TABLE IF EXISTS prs_traits')
c.execute('DROP TABLE IF EXISTS prs_weights')

c.execute('''
CREATE TABLE prs_traits (
    trait TEXT PRIMARY KEY,
    mean_score REAL,
    std_score REAL
)''')

c.execute('''
CREATE TABLE prs_weights (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    rsid TEXT,
    trait TEXT,
    risk_allele TEXT,
    effect_weight REAL
)''')

traits_data = [
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

c.executemany('INSERT INTO prs_traits VALUES (?, ?, ?)', traits_data)

weights_data = [
    ('rs11178997', 'Anxiety Disorder', 'G', 0.85),
    ('rs28399433', 'Generalized Anxiety Disorder', 'T', 0.85),
    ('rs1799853',  'Panic Disorder', 'T', 0.85),
    ('rs4244285',  'Post-Traumatic Stress Disorder', 'A', 0.85),
    ('rs1024611',  'Schizoaffective Disorder', 'A', 0.85),
    ('rs9272219',  'Bipolar Disorder', 'T', 0.85),
    ('rs1065852',  'Major Depressive Disorder', 'T', 0.85),
    ('rs10757278', 'Coronary Artery Disease', 'G', 0.50),
    ('rs7903146',  'Type 2 Diabetes', 'T', 0.50)
]

c.executemany('INSERT INTO prs_weights (rsid, trait, risk_allele, effect_weight) VALUES (?, ?, ?, ?)', weights_data)

conn.commit()
conn.close()
print('[✔] genomics.db seeded: 9 matched traits and weights.')
