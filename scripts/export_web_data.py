"""
Export de-identified data for the GitHub Pages interactive web app.

Outputs (all under docs/):
  docs/data/pca_data.json       — per-patient PCA coords + lab overlays (no PHI)
  docs/data/survival.json       — Kaplan-Meier curves per subgroup
  docs/data/model_meta.json     — scaler params + feature info for the calculator
  docs/model/subgroup_model.onnx — trained 8-feature XGBoost → ONNX
"""

import os, sys, json
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from xgboost import XGBClassifier
from lifelines import KaplanMeierFitter

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(BASE, "docs")
DATA = os.path.join(BASE, "data")

os.makedirs(os.path.join(DOCS, "data"), exist_ok=True)
os.makedirs(os.path.join(DOCS, "model"), exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load and filter data
# ---------------------------------------------------------------------------
print("Loading data...")
df = pd.read_csv(os.path.join(DATA, "processed", "AL_with_ccp_03.tsv"),
                 sep="\t", low_memory=False)

# AL only, with cluster assignment — set Code ID as index for later joins
df = df[df["Amyloid_type"] == "AL"].copy()
df = df.dropna(subset=["fna3_cluster_n"])
df = df.set_index("Code ID")
print(f"  {len(df)} AL patients with cluster labels")

PHI_COLS = [
    "Date_of_admission_to_center", "Date_of_birth",
    "Date_of_diagnosis", "Date_of_death", "Date_of_last_visit",
    "Cardiopulmonary_symptom_onset", "Renal_signs_symptom_onset",
    "Neurologic_symptom_onset", "Time_to_death_from_diagnosis__m_",
    "Amyloid_type", "Ethnicity", "Education",
]
df = df.drop(columns=[c for c in PHI_COLS if c in df.columns])

# ---------------------------------------------------------------------------
# 2. PCA on quantitative clustering variables
# ---------------------------------------------------------------------------
QVARS = [
    "Kappa_sFLC", "Lambda_sFLC", "kappa_lambda_ratio", "dFLC",
    "Bone_marrow_plasma_cells____",
    "IVSd", "LVEF",
    "Floors_climbed", "Distance_walked",
    "WBC", "Hemoglobin", "MCV", "Platelet_count",
    "Sed_rate", "CRP", "D_dimer", "Factor_X",
    "Creatinine", "eGFR", "24_hr_UTP", "Albumin", "Alk_phos",
    "BNP", "NT_proBNP", "Troponin",
    "Calcium", "Uric_acid", "LDH", "Iron", "TIBC",
    "Systolic_BP_sitting", "Diastolic_BP_sitting", "Pulse_sitting",
]
QVARS = [c for c in QVARS if c in df.columns]

print("Computing PCA...")
X_pca = df[QVARS].copy()
# Median imputation then z-score (matches paper's approach for visualisation)
imputer = SimpleImputer(strategy="median")
scaler  = StandardScaler()
X_imp   = scaler.fit_transform(imputer.fit_transform(X_pca))
pca     = PCA(n_components=2)
coords  = pca.fit_transform(X_imp)
print(f"  Variance explained: PC1={pca.explained_variance_ratio_[0]:.1%}, "
      f"PC2={pca.explained_variance_ratio_[1]:.1%}")

# Display name map
LABEL_MAP = {
    "Kappa_sFLC": "Kappa-sFLC (mg/L)",
    "Lambda_sFLC": "Lambda-sFLC (mg/L)",
    "kappa_lambda_ratio": "Kappa:Lambda ratio",
    "dFLC": "dFLC (mg/L)",
    "Bone_marrow_plasma_cells____": "BM Plasma Cells (%)",
    "IVSd": "IVSd (mm)",
    "LVEF": "LVEF (%)",
    "Floors_climbed": "Floors Climbed",
    "Distance_walked": "Distance Walked (ft)",
    "WBC": "WBC (/mm³)",
    "Hemoglobin": "Hemoglobin (g/dL)",
    "MCV": "MCV (µm³)",
    "Platelet_count": "Platelet Count (/mm²)",
    "Sed_rate": "ESR (mm/hr)",
    "CRP": "CRP (mg/L)",
    "D_dimer": "D-dimer (mg/L)",
    "Factor_X": "Factor X (%)",
    "Creatinine": "Creatinine (mg/dL)",
    "eGFR": "eGFR (mL/min/1.73m²)",
    "24_hr_UTP": "Proteinuria (mg/24h)",
    "Albumin": "Albumin (g/dL)",
    "Alk_phos": "Alkaline Phosphatase (U/L)",
    "BNP": "BNP (pg/mL)",
    "NT_proBNP": "NT-proBNP (pg/mL)",
    "Troponin": "Troponin (ng/mL)",
    "Calcium": "Calcium (mg/dL)",
    "Uric_acid": "Uric Acid (mg/dL)",
    "LDH": "LDH (U/L)",
    "Iron": "Iron (mcg/dL)",
    "TIBC": "TIBC (mcg/dL)",
    "Systolic_BP_sitting": "Systolic BP (mmHg)",
    "Diastolic_BP_sitting": "Diastolic BP (mmHg)",
    "Pulse_sitting": "Pulse (bpm)",
    "Age": "Age at Diagnosis (y)",
    "Mayo__2004_cardiac_biomarker_staging": "Mayo 2004 Cardiac Stage",
    "BU__BNP_based__cardiac_staging": "BU Cardiac Stage",
    "Renal Stage (Palladini)": "Renal Stage (Palladini)",
    "Sex": "Sex",
    "Race": "Race",
    "Primary_organ": "Primary Organ",
}

# Staging label maps (encoded integers → strings)
MAYO_MAP = {1: "Stage I", 2: "Stage II", 3: "Stage III", 4: "Stage IIIb"}
BU_MAP   = {1: "Stage I", 2: "Stage II", 3: "Stage III", 4: "Stage IIIb"}
SEX_MAP  = {1: "Male", 2: "Female"}
RACE_MAP = {1: "White", 2: "Black", 4: "Asian", 5: "Am. Indian/AK Native",
            6: "Native Hawaiian/Pacific", 7: "Multiracial", 8: "Unknown/Other"}
ORGAN_MAP = {1: "Cardiac", 2: "Renal", 3: "Hepatic", 4: "GI",
             5: "Pulmonary", 6: "ANS", 7: "PNS", 8: "Soft tissue", 9: "Other"}

# Overlay columns to include (quantitative values + categorical labels)
OVERLAY_COLS = QVARS + [
    "Age", "Sex", "Race", "Primary_organ",
    "Mayo__2004_cardiac_biomarker_staging",
    "BU__BNP_based__cardiac_staging",
    "Renal Stage (Palladini)",
]

print("Building PCA export...")
records = []
for i, (idx, row) in enumerate(df.iterrows()):
    rec = {
        "pc1": round(float(coords[i, 0]), 4),
        "pc2": round(float(coords[i, 1]), 4),
        "cluster": row["fna3_cluster_n"],
    }
    for col in OVERLAY_COLS:
        if col not in df.columns:
            continue
        val = row[col]
        # Apply label maps for categorical columns
        if col in ["Sex", "Race", "Primary_organ",
                   "Mayo__2004_cardiac_biomarker_staging",
                   "BU__BNP_based__cardiac_staging",
                   "Renal Stage (Palladini)"]:
            val = str(val) if pd.notna(val) else None
        else:
            val = round(float(val), 4) if pd.notna(val) else None
        rec[col] = val
    records.append(rec)

pca_export = {
    "variance_explained": {
        "pc1": round(float(pca.explained_variance_ratio_[0]), 4),
        "pc2": round(float(pca.explained_variance_ratio_[1]), 4),
    },
    "overlay_columns": [
        {"key": c, "label": LABEL_MAP.get(c, c),
         "type": "categorical" if c in ["Sex","Race","Primary_organ",
             "Mayo__2004_cardiac_biomarker_staging",
             "BU__BNP_based__cardiac_staging",
             "Renal Stage (Palladini)"] else "continuous"}
        for c in OVERLAY_COLS if c in df.columns
    ],
    "patients": records,
}

out_pca = os.path.join(DOCS, "data", "pca_data.json")
with open(out_pca, "w") as f:
    json.dump(pca_export, f, separators=(",", ":"))
print(f"  Wrote {out_pca} ({os.path.getsize(out_pca)//1024} KB)")

# ---------------------------------------------------------------------------
# 3. Kaplan-Meier survival curves per subgroup
# ---------------------------------------------------------------------------
print("Computing KM survival curves...")

# Use time_from_admission (years) and status (death=1)
kmf = KaplanMeierFitter()
CLUSTER_COLORS = {"Low": "#2ca02c", "Intermediate": "#ff7f0e", "High": "#d62728"}
CLUSTER_ORDER  = ["Low", "Intermediate", "High"]

survival_curves = {}
for cluster in CLUSTER_ORDER:
    mask = df["fna3_cluster_n"] == cluster
    t    = df.loc[mask, "time_from_admission"].dropna()
    e    = df.loc[mask, "status"].loc[t.index]
    kmf.fit(t, e, label=cluster)
    # Sample at regular time points up to 15 years
    timeline = np.linspace(0, 15, 300)
    sf  = kmf.survival_function_at_times(timeline)
    ci_df = kmf.confidence_interval_survival_function_
    # Interpolate CI at timeline points
    from scipy.interpolate import interp1d
    ci_times = ci_df.index.values
    ci_lo_fn = interp1d(ci_times, ci_df.iloc[:, 0].values, bounds_error=False, fill_value=(1.0, ci_df.iloc[-1, 0]))
    ci_hi_fn = interp1d(ci_times, ci_df.iloc[:, 1].values, bounds_error=False, fill_value=(1.0, ci_df.iloc[-1, 1]))
    ci_lo = ci_lo_fn(timeline)
    ci_hi = ci_hi_fn(timeline)
    survival_curves[cluster] = {
        "n": int(mask.sum()),
        "color": CLUSTER_COLORS[cluster],
        "timeline": [round(x, 3) for x in timeline.tolist()],
        "survival": [round(x, 4) for x in sf.values.tolist()],
        "ci_lower": [round(float(x), 4) for x in ci_lo.tolist()],
        "ci_upper": [round(float(x), 4) for x in ci_hi.tolist()],
        "median_survival": round(float(kmf.median_survival_time_), 2)
            if not np.isinf(kmf.median_survival_time_) else None,
        "survival_10yr": round(float(kmf.survival_function_at_times([10]).values[0]), 3),
    }
    print(f"  {cluster}: n={mask.sum()}, 10-yr survival="
          f"{survival_curves[cluster]['survival_10yr']:.1%}")

out_surv = os.path.join(DOCS, "data", "survival.json")
with open(out_surv, "w") as f:
    json.dump(survival_curves, f, separators=(",", ":"))
print(f"  Wrote {out_surv}")

# ---------------------------------------------------------------------------
# 4. Train 8-feature abbreviated model and export to ONNX
# ---------------------------------------------------------------------------
print("Training abbreviated model (8 features)...")

ABBR_FEATURES = [
    "Albumin", "BNP", "Diastolic_BP_sitting",
    "Sed_rate", "Systolic_BP_sitting", "TIBC",
    "Troponin", "X24_hr_UTP",
]
ABBR_DISPLAY = [
    "Albumin (g/dL)", "BNP (pg/mL)", "Diastolic BP (mmHg)",
    "ESR (mm/hr)", "Systolic BP (mmHg)", "TIBC (mcg/dL)",
    "Troponin (ng/mL)", "Proteinuria (mg/24h)",
]
ABBR_KEYS = [
    "albumin", "bnp", "diastolic_bp",
    "esr", "systolic_bp", "tibc",
    "troponin", "proteinuria",
]

# Load abbreviated dataset (already prepared)
Xi_abbr = pd.read_csv(os.path.join(DATA, "processed", "Xi_abbr.tsv"), sep="\t", index_col=0)

# Align with cluster labels and drop missing
y_ser = df["fna3_cluster_n"].map({"Low": 0, "Intermediate": 1, "High": 2})
common_idx = Xi_abbr.index.intersection(y_ser.dropna().index)
X_ab = Xi_abbr.loc[common_idx]
y_ab = y_ser.loc[common_idx]
print(f"  Training set: {len(X_ab)} patients, {X_ab.shape[1]} features")

# Build pipeline: median impute → standard scale → XGBoost
pipe = Pipeline([
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler",  StandardScaler()),
    ("clf",     XGBClassifier(
        n_estimators=300,
        max_depth=4,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        use_label_encoder=False,
        eval_metric="mlogloss",
        random_state=42,
        n_jobs=-1,
    )),
])
pipe.fit(X_ab.values, y_ab.values)

# Quick accuracy check
from sklearn.model_selection import cross_val_score
cv_acc = cross_val_score(pipe, X_ab.values, y_ab.values, cv=5, scoring="accuracy")
print(f"  5-fold CV accuracy: {cv_acc.mean():.3f} ± {cv_acc.std():.3f}")

# Export scaler params and feature info as JSON (for JS-side normalization reference)
sc = pipe.named_steps["scaler"]
imp = pipe.named_steps["imputer"]
model_meta = {
    "features": [
        {
            "key": ABBR_KEYS[i],
            "col": ABBR_FEATURES[i],
            "label": ABBR_DISPLAY[i],
            "mean": round(float(sc.mean_[i]), 4),
            "std": round(float(sc.scale_[i]), 4),
            "median": round(float(imp.statistics_[i]), 4),
        }
        for i in range(len(ABBR_FEATURES))
    ],
    "classes": ["Low", "Intermediate", "High"],
    "cv_accuracy": round(float(cv_acc.mean()), 3),
    "cluster_colors": CLUSTER_COLORS,
}
out_meta = os.path.join(DOCS, "data", "model_meta.json")
with open(out_meta, "w") as f:
    json.dump(model_meta, f, indent=2)
print(f"  Wrote {out_meta}")

# Export XGBoost booster to ONNX
# Use onnxmltools to register the XGBoost converter, then convert the full pipeline
print("  Exporting to ONNX...")
from skl2onnx import convert_sklearn, update_registered_converter
from skl2onnx.common.data_types import FloatTensorType
from skl2onnx.common.shape_calculator import calculate_linear_classifier_output_shapes
from onnxmltools.convert.xgboost.operator_converters.XGBoost import convert_xgboost

update_registered_converter(
    XGBClassifier,
    "XGBoostXGBClassifier",
    calculate_linear_classifier_output_shapes,
    convert_xgboost,
    options={"nocl": [True, False], "zipmap": [True, False, "columns"]},
)

onnx_model = convert_sklearn(
    pipe,
    initial_types=[("float_input", FloatTensorType([None, len(ABBR_FEATURES)]))],
    target_opset={"": 12, "ai.onnx.ml": 3},
    options={XGBClassifier: {"nocl": True, "zipmap": False}},
)
out_onnx = os.path.join(DOCS, "model", "subgroup_model.onnx")
with open(out_onnx, "wb") as f:
    f.write(onnx_model.SerializeToString())
print(f"  Wrote {out_onnx} ({os.path.getsize(out_onnx)//1024} KB)")

print("\nDone! All files written to docs/")
