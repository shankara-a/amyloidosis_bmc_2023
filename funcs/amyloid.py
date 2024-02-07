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
# Rename Dataset Vars
# -------------------------
# def fix_string(s):
#     s = s.strip()
#     for ch in [" ","-",":","%","(",")","/"]:
#         s = s.replace(ch,"_")
#     return s

ddict_unclean = {'Date_of_admission_to_center': 'Date of admission to center',
 'Date_of_birth': 'Date of birth',
 'Date_of_diagnosis': 'Date of diagnosis',
 'Date_of_death': 'Date of death',
 'Date_of_last_visit': 'Date of last visit',
 'Cardiopulmonary_symptom_onset': 'Cardiopulmonary symptom onset',
 'Renal_signs_symptom_onset': 'Renal signs/symptom onset',
 'Neurologic_symptom_onset': 'Neurologic symptom onset',
 'Amyloid_type': 'Amyloid type',
 'Age': 'Age',
 'Time_to_death_from_diagnosis__m_': 'Time to death from diagnosis (m)',
 'Sex': 'Sex',
 'Ethnicity': 'Ethnicity',
 'Race': 'Race',
 'Education': 'Education',
 'Abdominal_fat_pad_CR_staining': 'Abdominal fat pad CR staining',
 'Bone_marrow_CR_staining': 'Bone marrow CR staining',
 'Kappa_or_lambda_PCD': 'Kappa or lambda PCD',
 'Kappa_sFLC': 'Kappa-sFLC',
 'Lambda_sFLC': 'Lambda-sFLC',
 'kappa_lambda_ratio': 'kappa:lambda ratio',
 'dFLC': 'dFLC',
 'Heart': 'Heart',
 'Kidney': 'Kidney',
 'Liver': 'Liver',
 'GI_Tract': 'GI Tract',
 'Lung': 'Lung',
 'Neurologic': 'Neurologic',
 '_Autonomic_': '(Autonomic)',
 '_Peripheral_': '(Peripheral)',
 'Soft_tissue_MSK': 'Soft tissue/MSK',
 'Dermatologic': 'Dermatologic',
 'Primary_organ': 'Primary organ',
 'Secondary_organ': 'Secondary organ',
 'Weight_loss': 'Weight loss',
 'Fatigue': 'Fatigue',
 'Exertional_dyspnea': 'Exertional dyspnea',
 'Dyspnea_at_rest': 'Dyspnea at rest',
 'Orthopnea_PND': 'Orthopnea/PND',
 'Leg_swelling__symptom_': 'Leg swelling (symptom)',
 'Peripheral_edema_on_exam': 'Peripheral edema on exam',
 'ascites_symptom': 'ascites symptom',
 'Ascites_on_exam': 'Ascites on exam',
 'Orthostasis': 'Orthostasis',
 'Dizzy_spells': 'Dizzy spells',
 'Syncope': 'Syncope',
 'Ecchymosis': 'Ecchymosis',
 'Petechiae': 'Petechiae',
 'Periorbital_ecchymosis': 'Periorbital ecchymosis',
 'Anorexia': 'Anorexia',
 'Abdominal_pain': 'Abdominal pain',
 'Large_tongue': 'Large tongue',
 'Carpal_tunnel': 'Carpal tunnel',
 'Hypertension': 'Hypertension',
 'Arrhythmia': 'Arrhythmia ',
 'WBC': 'WBC',
 'Hemoglobin': 'Hemoglobin',
 'MCV': 'MCV',
 'Platelet_count': 'Platelet count',
 'Sed_rate': 'Sed rate',
 'CRP': 'CRP',
 'D_dimer': 'D-dimer',
 'Factor_X': 'Factor X',
 'Creatinine': 'Creatinine',
 '24_hr_UTP': '24-hr UTP',
 'Albumin': 'Albumin',
 'Alk_phos': 'Alk phos',
 'BNP': 'BNP',
 'NT_proBNP': 'NT-proBNP',
 'Troponin': 'Troponin',
 'Calcium': 'Calcium',
 'Uric_acid': 'Uric acid',
 'LDH': 'LDH',
 'Iron': 'Iron',
 'TIBC': 'TIBC',
 'Bone_marrow_plasma_cells____': 'Bone marrow plasma cells (%)',
 'SIFE_M_component': 'SIFE M-component',
 'UIFE_M_component': 'UIFE M-component',
 'Mayo__2004_cardiac_biomarker_staging': 'Mayo  2004 cardiac biomarker staging',
 'BU__BNP_based__cardiac_staging': 'BU (BNP-based) cardiac staging',
 'IVSd': 'IVSd',
 'LVEF': 'LVEF',
 'Floors_climbed': 'Floors climbed',
 'Distance_walked': 'Distance walked',
 'Systolic_BP_sitting': 'Systolic BP sitting',
 'Diastolic_BP_sitting': 'Diastolic BP sitting',
 'Pulse_sitting': 'Pulse sitting',
 'status': 'status',
 'time': 'time'}

ddict_clean = {'Date of admission to center': 'Date_of_admission_to_center',
 'Date of birth': 'Date_of_birth',
 'Date of diagnosis': 'Date_of_diagnosis',
 'Date of death': 'Date_of_death',
 'Date of last visit': 'Date_of_last_visit',
 'Cardiopulmonary symptom onset': 'Cardiopulmonary_symptom_onset',
 'Renal signs/symptom onset': 'Renal_signs_symptom_onset',
 'Neurologic symptom onset': 'Neurologic_symptom_onset',
 'Amyloid type': 'Amyloid_type',
 'Age': 'Age',
 'Time to death from diagnosis (m)': 'Time_to_death_from_diagnosis__m_',
 'Sex': 'Sex',
 'Ethnicity': 'Ethnicity',
 'Race': 'Race',
 'Education': 'Education',
 'Abdominal fat pad CR staining': 'Abdominal_fat_pad_CR_staining',
 'Bone marrow CR staining': 'Bone_marrow_CR_staining',
 'Kappa or lambda PCD': 'Kappa_or_lambda_PCD',
 'Kappa-sFLC': 'Kappa_sFLC',
 'Lambda-sFLC': 'Lambda_sFLC',
 'kappa:lambda ratio': 'kappa_lambda_ratio',
 'dFLC': 'dFLC',
 'Heart': 'Heart',
 'Kidney': 'Kidney',
 'Liver': 'Liver',
 'GI Tract': 'GI_Tract',
 'Lung': 'Lung',
 'Neurologic': 'Neurologic',
 '(Autonomic)': '_Autonomic_',
 '(Peripheral)': '_Peripheral_',
 'Soft tissue/MSK': 'Soft_tissue_MSK',
 'Dermatologic': 'Dermatologic',
 'Primary organ': 'Primary_organ',
 'Secondary organ': 'Secondary_organ',
 'Weight loss': 'Weight_loss',
 'Fatigue': 'Fatigue',
 'Exertional dyspnea': 'Exertional_dyspnea',
 'Dyspnea at rest': 'Dyspnea_at_rest',
 'Orthopnea/PND': 'Orthopnea_PND',
 'Leg swelling (symptom)': 'Leg_swelling__symptom_',
 'Peripheral edema on exam': 'Peripheral_edema_on_exam',
 'ascites symptom': 'ascites_symptom',
 'Ascites on exam': 'Ascites_on_exam',
 'Orthostasis': 'Orthostasis',
 'Dizzy spells': 'Dizzy_spells',
 'Syncope': 'Syncope',
 'Ecchymosis': 'Ecchymosis',
 'Petechiae': 'Petechiae',
 'Periorbital ecchymosis': 'Periorbital_ecchymosis',
 'Anorexia': 'Anorexia',
 'Abdominal pain': 'Abdominal_pain',
 'Large tongue': 'Large_tongue',
 'Carpal tunnel': 'Carpal_tunnel',
 'Hypertension': 'Hypertension',
 'Arrhythmia ': 'Arrhythmia',
 'WBC': 'WBC',
 'Hemoglobin': 'Hemoglobin',
 'MCV': 'MCV',
 'Platelet count': 'Platelet_count',
 'Sed rate': 'Sed_rate',
 'CRP': 'CRP',
 'D-dimer': 'D_dimer',
 'Factor X': 'Factor_X',
 'Creatinine': 'Creatinine',
 '24-hr UTP': '24_hr_UTP',
 'Albumin': 'Albumin',
 'Alk phos': 'Alk_phos',
 'BNP': 'BNP',
 'NT-proBNP': 'NT_proBNP',
 'Troponin': 'Troponin',
 'Calcium': 'Calcium',
 'Uric acid': 'Uric_acid',
 'LDH': 'LDH',
 'Iron': 'Iron',
 'TIBC': 'TIBC',
 'Bone marrow plasma cells (%)': 'Bone_marrow_plasma_cells____',
 'SIFE M-component': 'SIFE_M_component',
 'UIFE M-component': 'UIFE_M_component',
 'Mayo  2004 cardiac biomarker staging': 'Mayo__2004_cardiac_biomarker_staging',
 'BU (BNP-based) cardiac staging': 'BU__BNP_based__cardiac_staging',
 'IVSd': 'IVSd',
 'LVEF': 'LVEF',
 'Floors climbed': 'Floors_climbed',
 'Distance walked': 'Distance_walked',
 'Systolic BP sitting': 'Systolic_BP_sitting',
 'Diastolic BP sitting': 'Diastolic_BP_sitting',
 'Pulse sitting': 'Pulse_sitting',
 'status': 'status',
 'time': 'time'}

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

# enc_treatment_metagroups = {
#     "HDM/SCT":"HDM/SCT",
#     "Melphalan-based regimen":"Melphalan-based regimen",
#     "Proteasome inhibitor-based":"Proteasome inhibitor-based",
#     "Dara-CyBorD":"Daratumumab-based",
#     "Daratumumab-based":"Daratumumab-based",
#     "IMiD-based":"IMiD-based",


# }

# -------------------------
# Final Tables
# -------------------------
tableone_names = {
    'Kappa-sFLC':'K (mg/L)',
    'Lambda-sFLC':'L (mg/L)', 
    'kappa:lambda ratio':'K:L',
    'dFLC':'dFLC mg/L',
    'Bone marrow plasma cells (%)':'BMPC (%)', 
    'IVSd':'IVSd (mm)', 
    'LVEF':'LVEF (%)', 
    'WBC':'WBC /mm^3', 
    'Hemoglobin':'Hgb (g/dL)',
    'MCV':'MCV (um^3)', 
    'Platelet count':'PLT /mm^2', 
    'Sed rate': 'ESR (mm/hr)', 
    'D-dimer':'D-dimer (mg/L)', 
    'Factor X':'Factor X (%)', 
    'eGFR':'eGFR (mL/min/1.72m^2)',
    '24-hr UTP':'Proteinuria (mg/24h)', 
    'Albumin':'Albumin (g/dL)', 
    'Alk phos':'ALP (U/L)', 
    'BNP':'BNP (pg/mL)', 
    'Troponin':'Troponin (ng/mL)', 
    'Calcium': 'Calcium (mg/dL)',
    'Uric acid': 'Uric Acid (mg/dL)', 
    'LDH':'LDH (U/L)', 
    'Iron':'Fe (mcg/dL)', 
    'TIBC':'TIBC (mcg/dL)',
    "Kappa or lambda PCD":"Amyloidogenic LC"
}