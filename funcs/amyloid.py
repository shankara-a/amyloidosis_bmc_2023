# -- import packages: --------------------------------------------------------------------

# -------------------------
# Dataset Features
# -------------------------
# Note changed {Age at diagnosis (y): Age}
demographics = ['Age', 'Sex', 'Ethnicity', 'Race', 'Education'] 


dates = [
    'Date of admission to center', 'Date of birth', 'Date of diagnosis', 
    'Date of death', 'Date of last visit', 'Cardiopulmonary symptom onset', 
    'Renal signs/symptom onset', 'Neurologic symptom onset'
]

amyloid_vars = [
    'Amyloid type','Abdominal fat pad CR staining',
    'Bone marrow CR staining', 'Kappa or lambda PCD', 'Kappa-sFLC',
    'Lambda-sFLC', 'kappa:lambda ratio', 'dFLC',
    'Primary organ', 'Secondary organ',
    'Bone marrow plasma cells (%)',
    'SIFE M-component', 'UIFE M-component'
]

amyloid_ros = [
    'Heart', 'Kidney', 'Liver',
    'GI Tract', 'Lung', 'Neurologic', '(Autonomic)', '(Peripheral)',
    'Soft tissue/MSK', 'Dermatologic',
]

amyloid_symptoms = [
    'Weight loss', 'Fatigue', 'Exertional dyspnea', 'Dyspnea at rest', 'Orthopnea/PND',
    'Leg swelling (symptom)', 'Peripheral edema on exam', 'ascites symptom',
    'Ascites on exam', 'Orthostasis', 'Dizzy spells', 'Syncope',
    'Ecchymosis', 'Petechiae', 'Periorbital ecchymosis', 'Anorexia',
    'Abdominal pain', 'Large tongue', 'Carpal tunnel', 'Hypertension',
    'Arrhythmia '
]
    
labs = [
    'WBC', 'Hemoglobin', 'MCV', 'Platelet count', 'Sed rate',
    'CRP', 'D-dimer', 'Factor X', 'Creatinine', '24-hr UTP', 'Albumin',
    'Alk phos', 'BNP', 'NT-proBNP', 'Troponin', 'Calcium', 'Uric acid',
    'LDH', 'Iron', 'TIBC',
]

staging = [
    'Mayo  2004 cardiac biomarker staging',
    'BU (BNP-based) cardiac staging', 
]

other = ['IVSd', 'LVEF', 'Floors climbed', 'Distance walked', ]

vitals = ['Systolic BP sitting', 'Diastolic BP sitting','Pulse sitting']

# Quantitative Independent Variables
# Includes amyloid lab values, BM plasma cell, cardiac biomarkers, lab values, floors climbed, distance walked
qvars = [
    'Kappa-sFLC', 'Lambda-sFLC', 'kappa:lambda ratio', 'dFLC', 'Bone marrow plasma cells (%)',
    'IVSd', 'LVEF', 'Floors climbed','Distance walked',
    'WBC', 'Hemoglobin', 'MCV', 'Platelet count', 'Sed rate',
    'CRP', 'D-dimer', 'Factor X', 'Creatinine', '24-hr UTP', 'Albumin',
    'Alk phos', 'BNP', 'NT-proBNP', 'Troponin', 'Calcium', 'Uric acid',
    'LDH', 'Iron', 'TIBC'
] + vitals

# Categorical Independent Variables
# Includes demographics, amyloid testing, amyloid ROS, and amyloid related symptoms
catvars = demographics + list(set(amyloid_vars) - set(qvars)) + amyloid_ros + amyloid_symptoms

# -------------------------
# Dataset Encodings
# -------------------------
enc_amyloid_type = {
    1: "AL",
    2: "MM_associated_AL",
    14: "WM_associated_AL",
    15: "B-cell_lymphoproliferative_associated_AL"
}

enc_sex = {1: "male", 2:"female"}

enc_ethnicity = {1: "hispanic", 2:"non_hispanic"}

enc_race = {
    1: "White",
    2: "Black",
    4: "Asian",
    5: "American_Indian_Alaska_Native",
    6: "Native_Hawaiian_Pacific",
    7: "Multiracial",
    8: "Unknown/other"
}

enc_education = {
    1: "graduate",
    2: "college",
    3: "15years",
    4: "14years",
    5: "13years",
    6: "high_school_or_GED",
    7: "less_than_highschool",
    10: "unknown"
}

# For both abdominal fat pad & BM
enc_staining = {
    1: "neg",
    2: "pos",
    3: "pos_CR1+",
    4: "pos_CR2+",
    5: "pos_CR3+",
    6: "not_done"
}

enc_ros = {
    1: "not_involved",
    2: "involved",
    3: "equivocal"
}

# For both primary and secondary organ field
enc_organ = {
    1: "cardiac", 
    2: "renal", 
    3: "hepatic", 
    4: "gi", 
    5: "pulm",
    6: "ans", 
    7: "pns", 
    8: "soft_tissue", 
    9: "other"
}

# For all descriptive symptoms (ie amyloid_symptoms)
enc_descriptors = {1: "no", 2: "yes", 3: "uncertain"}

# For SIFE and UIFE
enc_m_components = {1: "no", 2: "yes", 3: "not_done"}

enc_mayo_2004 = {
    1: "stage I",
    2: "stage II",
    3: "stage III",
    4: "stage IIIb"
}

enc_bu_2019 = {
    1: "stage I",
    2: "stage II",
    3: "stage III",
    4: "stage IIIb"
}

# -------------------------
# Treatment Encodings
# -------------------------

enc_treaments = """
9	Daratumumab, pomalidomide, dex (DPd)
10	Venetoclax
11	Melphalan, IV/SCT
12	Melphalan, cyclic oral
13	Melphalan, daily oral
14	Melphalan, other
15	Decadron
16	Revlimid
17	Velcade
18	Thalidomide
19	Mini-Melphalan
23	Mel-Dex
25	Lenalidomide
26	Pomalidomide
27	Ixazomib
28	Daratumumab + other anti-CD38 therapy
29	Carfilzomib
30	Bendamustine
31	CyBorD
32	VMD
36  Bortezomib, SCT
37	Vel/Mel/Dex
39	CyBorD daratumumab
38	Mel/Rev/Dex
51	CAPD (peritoneal dialysis)
52	Haemodialysis
53	Heart Transplant
54	Kidney Transplant
55	Liver Transplant
56	Pacemaker
57	Defibrillator
99	Other
58	Bendamustine + rituximab (R-benda)
59	cycloposhamide 
60	Revlimid, velcade, dex
61	bortezomib, dex, rituximab
62	Other rituximab-based regimen
63	dara + IMiD
64	BTK inhibitor
"""

enc_treatment_groups = """
HDM/SCT	11,36
Melphalan-based regimen	12,13,14,19,23
Proteasome inhibitor-based	17,27,29,31,32,37
Dara-CyBorD	39
Daratumumab-based	9,28,63
IMiD-based	16,18,25,26,38
RVD	60
Venetoclax-based	10
Cyclophosphamide	59
Bendamustine	30
Glucocorticoid monotherapy	15
R-benda	58
BDR	61
Other rituximab-based regimen	62
BTK inhibitor	64
Supportive	51,52,53,54,55,56,57,58
"""