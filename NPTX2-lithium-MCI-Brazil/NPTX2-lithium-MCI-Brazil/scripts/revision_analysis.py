#!/usr/bin/env python3
"""
================================================================================
NPTX2 / lithium MCI trial  --  complete reanalysis for the 2026 revision
================================================================================

Xiao M-F*, Kim K*, et al.
"Lithium upregulates NPTX2 in association with cognitive maintenance."

This script reproduces, end to end, every analysis and figure prepared in
response to the reviewers of Alzheimer's Research & Therapy.

Run order is top to bottom.  Every table is written as .xlsx and every figure
as vector .pdf under ./outputs/.

--------------------------------------------------------------------------------
REQUIRED INPUT FILES  (place beside this script, or edit INPUT below)
--------------------------------------------------------------------------------
  nulisa.xlsx       Alamar NULISA export.  Sheets:
                      'NULISA data'            long format, one row per
                                               sample x target, with columns
                                               PlateID, SampleName, SampleType,
                                               Target, SampleQC, LOD, NPQ
                      'CSF sample spreadsheet' baseline / one-year ID pairing
                                               (header on row 4)
                      'list of analytes'       target annotation
  clinical.xlsx     Trial clinical database.  Sheets:
                      'Banco de dados'         three timepoint blocks
                      're-organized by MX'     authoritative treatment arm
  prm.xlsx          PRM-MS protein-level table (Sample Number, Group, targets)
  master_old.xlsx   Previous merged analysis table, used only to verify that
                    this pipeline reproduces the original values

--------------------------------------------------------------------------------
METHOD SUMMARY
--------------------------------------------------------------------------------
  Analyte inclusion   NULISA detectability >= 50% above the plate-specific LOD
                      (vendor recommendation).  Values below LOD but non-zero
                      are retained.  APOE4 is a genotype indicator, not a
                      continuous analyte, and is analyzed separately.
  Sample inclusion    >= 70% of targets above LOD (vendor threshold for CSF).
  Scale               NULISA NPQ (log2) and PRM (log2) as provided.  ELISA is
                      log2-transformed so all three platforms are comparable.
  Multiple testing    Benjamini-Hochberg within platform (primary); pooled
                      across all analytes (sensitivity).
  Between-arm test    ANCOVA on one-year change, baseline value as covariate.
  Group difference    Fisher r-to-z plus a formal analyte x arm interaction.
  Network tests       Label permutation and size-matched down-sampling, with a
                      threshold-independent mean |Fisher z| measure.
================================================================================
"""

from __future__ import annotations
import os, re, warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# ==============================================================================
# CONFIGURATION
# ==============================================================================
INPUT = {"nulisa": "nulisa.xlsx", "clinical": "clinical.xlsx",
         "prm": "prm.xlsx", "master_old": "master_old.xlsx"}

OUT_T, OUT_F = "outputs/tables", "outputs/figures"
os.makedirs(OUT_T, exist_ok=True)
os.makedirs(OUT_F, exist_ok=True)

DETECTABILITY_MIN = 0.50      # analyte inclusion, vendor recommendation
SAMPLE_QC_MIN     = 0.70      # sample inclusion for CSF, vendor threshold
N_BOOT, N_PERM, N_DOWNSAMPLE = 2000, 5000, 1000
SEED = 0

mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                     "font.family": "DejaVu Sans", "axes.linewidth": 0.7,
                     "xtick.major.width": 0.7, "ytick.major.width": 0.7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.labelsize": 8})

# Palette validated for colour-vision deficiency separation
ARM_STYLE = {"Lithium": {"c": "#2E75B6", "m": "o"},
             "Placebo": {"c": "#C55A11", "m": "^"}}
C_NULL, C_NOM, C_FDR = "#B8BEC6", "#2E75B6", "#C0392B"
PLATFORM_MARKER = {"NULISA": "o", "PRM": "s", "ELISA": "D"}


def save_table(obj, name):
    path = f"{OUT_T}/{name}.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as w:
        if isinstance(obj, dict):
            for sheet, df in obj.items():
                df.to_excel(w, sheet_name=str(sheet)[:31], index=False)
        else:
            obj.to_excel(w, index=False)
    print("  saved:", path)
    return path


def save_fig(fig, name):
    path = f"{OUT_F}/{name}.pdf"
    fig.savefig(path, format="pdf", bbox_inches="tight")
    print("  saved:", path)
    return path


def platform_of(col: str) -> str:
    if col.endswith("_NULISA"):
        return "NULISA"
    if col.endswith("_PRM"):
        return "PRM"
    return "ELISA"


def fisher_ci(r, n):
    if n < 5 or abs(r) >= 1:
        return (np.nan, np.nan)
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    return tuple(np.tanh([z - 1.96 * se, z + 1.96 * se]))


def bh(pvals):
    p = pd.Series(pvals, dtype=float)
    q = pd.Series(np.nan, index=p.index)
    ok = p.notna()
    if ok.sum():
        q.loc[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return q.values


# ==============================================================================
# STEP 1  --  Sample keys, treatment arm, timepoint assignment
# ==============================================================================
print("\nSTEP 1  sample keys and treatment arm")

nul = pd.read_excel(INPUT["nulisa"], "NULISA data", header=0)

pair = pd.read_excel(INPUT["nulisa"], "CSF sample spreadsheet", header=3)
pair = pair.dropna(how="all").dropna(axis=1, how="all")
pair.columns = ["base_id", "y1_id"] + list(pair.columns[2:])
pair = pair[["base_id", "y1_id"]].dropna().astype(int)

# Timepoint comes ONLY from the pairing table.  Sample-name suffixes do not
# encode timepoint and must not be used.
nul["csf_id"] = nul["SampleName"].str.extract(r"CSF_(\d+)").astype(int)
samples = nul[["SampleName", "csf_id", "PlateID"]].drop_duplicates()
base_set, y1_set = set(pair.base_id), set(pair.y1_id)
samples["timepoint"] = samples.csf_id.map(
    lambda i: "base" if i in base_set else ("1y" if i in y1_set else "UNMAPPED"))
assert (samples.timepoint == "UNMAPPED").sum() == 0, "unmapped NULISA sample"

cl = pd.read_excel(INPUT["clinical"], "Banco de dados", header=0)
cl_id0 = pd.to_numeric(cl.iloc[:, 0], errors="coerce")

# participant_id: the baseline CSF id.  Participants without a follow-up sample
# (and therefore absent from the pairing table) map to themselves.
id2pid = {**{b: b for b in pair.base_id},
          **{y: b for b, y in zip(pair.base_id, pair.y1_id)}}
for b in cl_id0.dropna().astype(int):
    id2pid.setdefault(b, b)
samples["participant_id"] = samples.csf_id.map(id2pid)

# Treatment arm comes from the re-organized sheet, which is complete.
r2 = pd.read_excel(INPUT["clinical"], "re-organized by MX", header=1)
r2.columns = [str(c).strip() for c in r2.columns]
arm_map = pd.DataFrame({"participant_id": pd.to_numeric(r2.iloc[:, 0], errors="coerce"),
                        "arm_code": r2["1-Li_0-placebo"]}).dropna(subset=["participant_id"])
arm_map["participant_id"] = arm_map.participant_id.astype(int)
arm_map["arm"] = arm_map.arm_code.map({1: "Lithium", 0: "Placebo"})
print("  arm:", arm_map.arm.value_counts().to_dict())


# ==============================================================================
# STEP 2  --  NULISA quality control
# ==============================================================================
print("\nSTEP 2  NULISA quality control")

# LOD is reported per target AND per plate, so detectability is evaluated
# row-wise.  A missing LOD marks a target for which the vendor does not report
# one; such targets are retained with an explicit note.
nul["above_lod"] = np.where(nul.LOD.isna(), np.nan, (nul.NPQ > nul.LOD).astype(float))

sample_det = nul.groupby("SampleName")["above_lod"].mean()
qc = nul.groupby("Target").agg(det_rate=("above_lod", "mean"),
                               n_obs=("NPQ", "size"),
                               lod_reported=("LOD", lambda s: s.notna().any())).reset_index()
qc["floor_tie"] = nul.groupby("Target")["NPQ"].apply(lambda s: (s == s.min()).mean()).values


def analyte_decision(row):
    if row.Target == "APOE4":
        return "CATEGORICAL (genotype indicator)"
    if not row.lod_reported:
        return "INCLUDE (LOD not reported by vendor)"
    return "INCLUDE" if row.det_rate >= DETECTABILITY_MIN else "EXCLUDE (detectability < 50%)"


qc["decision"] = qc.apply(analyte_decision, axis=1)
keep_nulisa = qc.loc[qc.decision.str.startswith("INCLUDE"), "Target"].tolist()
print("  analytes:", qc.decision.value_counts().to_dict())
print("  samples failing the %.0f%% detectability threshold: %s"
      % (SAMPLE_QC_MIN * 100, list(sample_det[sample_det < SAMPLE_QC_MIN].index)))

nul_wide = nul.pivot_table(index="SampleName", columns="Target", values="NPQ")
nul_wide = nul_wide.join(samples.set_index("SampleName")[
    ["csf_id", "participant_id", "timepoint", "PlateID"]])
nul_wide["nulisa_sample_qc_pass"] = nul_wide.index.map(
    lambda s: sample_det[s] >= SAMPLE_QC_MIN)


# ==============================================================================
# STEP 3  --  Clinical table, PRM, merged metadata
# ==============================================================================
print("\nSTEP 3  merged metadata")

RENAME = {"ID amostra": "csf_id", "Sexo": "sex", "Data de Nasc.": "dob", "Idade": "age",
          "Date": "visit_date", "CDR": "CDR", "CDR_SB": "CDR_SB", "ADAS": "ADAS",
          "SLN": "SLN", "TMTA": "TMTA", "TMTB": "TMTB",
          "Delayed_rec": "delayed_rec", "Figure_rec": "figure_rec",
          "APOE": "apoe_genotype", "1-Li_0-placebo": "arm_code", "Litemia": "litemia",
          "NPTXR": "NPTXR_ELISA", "Ab pg/mL": "Ab_ELISA", "Tau pg/mL": "Tau_ELISA",
          "pTau pg/mL": "pTau_ELISA", "Conversao": "conversion", "CIBIC": "CIBIC",
          "Data coleta": "collect_date"}


def clean_block(df, tp):
    out, names, seen = df.copy(), [], {}
    for c in out.columns:
        n = re.sub(r"\.\d+$", "", str(c)).strip()
        if n.startswith("CSF NPTX2"):
            n = "NPTX2_ELISA"
        elif n.startswith("Raz") and "NPTX2" in n:
            n = "NPTX2_NPTXR_ratio"
        elif n.startswith("Raz"):
            n = "Ab_pTau_ratio"
        elif n.startswith("Perfil"):
            n = "AD_profile"
        elif n.startswith("Convers"):
            n = "conversion"
        n = RENAME.get(n, n)
        seen[n] = seen.get(n, 0) + 1
        names.append(n if seen[n] == 1 else f"{n}__{seen[n]}")
    out.columns = names
    out["timepoint"] = tp
    return out


long = pd.concat([clean_block(cl.iloc[:, s:e], tp)
                  for s, e, tp in [(0, 26, "base"), (26, 48, "1y"), (48, 70, "3y")]],
                 ignore_index=True)
long["csf_id"] = pd.to_numeric(long["csf_id"], errors="coerce")
long = long.dropna(subset=["csf_id"])
long["csf_id"] = long.csf_id.astype(int)
long["age"] = pd.to_numeric(long.age.astype(str).str.replace("~", "").str.strip(),
                            errors="coerce")
long["litemia"] = pd.to_numeric(long.litemia, errors="coerce")
long["participant_id"] = long.csf_id.map(id2pid)

# Participant-level variables are taken from the baseline block and propagated.
static = (long[long.timepoint == "base"]
          .set_index("participant_id")[["sex", "dob", "apoe_genotype", "AD_profile", "age"]]
          .rename(columns={"age": "age_base"}))
long = long.drop(columns=[c for c in ["sex", "dob", "apoe_genotype", "AD_profile"]
                          if c in long], errors="ignore")
long = long.merge(static.reset_index(), on="participant_id", how="left")
long["apoe4_carrier"] = long.apoe_genotype.map(
    lambda g: np.nan if pd.isna(g) else int("4" in str(g).replace("Ɛ", "")))

prm = pd.read_excel(INPUT["prm"])
prm.columns = [str(c).strip() for c in prm.columns]
prm["csf_id"] = pd.to_numeric(prm["Sample Number"], errors="coerce").astype("Int64")
prm_targets = [c for c in prm.columns if c not in ("Sample Number", "Group", "csf_id")]
prm2 = prm.rename(columns={c: c + "_PRM" for c in prm_targets}).drop(
    columns=["Group", "Sample Number"])

nul_keep = nul_wide[keep_nulisa + ["csf_id", "participant_id", "timepoint",
                                   "PlateID", "nulisa_sample_qc_pass"]].copy()
nul_keep = nul_keep.rename(columns={c: c + "_NULISA" for c in keep_nulisa})

meta = (long.drop(columns=["arm_code"], errors="ignore")
        .merge(arm_map[["participant_id", "arm"]], on="participant_id", how="left")
        .merge(nul_keep.reset_index().drop(columns=["participant_id", "timepoint"]),
               on="csf_id", how="left")
        .merge(prm2, on="csf_id", how="left"))

NUMCOLS = ["age_base", "CDR", "CDR_SB", "ADAS", "SLN", "TMTA", "TMTB", "delayed_rec",
           "figure_rec", "litemia", "NPTX2_ELISA", "NPTXR_ELISA", "Ab_ELISA",
           "Tau_ELISA", "pTau_ELISA", "AD_profile", "apoe4_carrier", "CIBIC"]
for c in NUMCOLS:
    if c in meta.columns:
        meta[c] = pd.to_numeric(meta[c], errors="coerce")
print("  meta:", meta.shape)
print(meta.groupby(["arm", "timepoint"]).size().to_string())


# ==============================================================================
# STEP 4  --  Derived variables:  scale, direction, change scores
# ==============================================================================
print("\nSTEP 4  derived variables")

# Cognitive outcomes are harmonised so that a positive change means worsening.
DIRECTION = {"CDR": +1, "CDR_SB": +1, "ADAS": +1, "TMTA": +1, "TMTB": +1,
             "SLN": -1, "delayed_rec": -1, "figure_rec": -1}
for c, sgn in DIRECTION.items():
    meta[c + "_h"] = pd.to_numeric(meta[c], errors="coerce") * sgn
COG_H = [c + "_h" for c in DIRECTION]

ELISA = [c for c in ["NPTX2_ELISA", "NPTXR_ELISA", "Ab_ELISA", "Tau_ELISA", "pTau_ELISA"]
         if c in meta.columns]
ANALYTES = ([c + "_NULISA" for c in keep_nulisa] +
            [c + "_PRM" for c in prm_targets] + ELISA)
ANALYTES = [c for c in ANALYTES if c in meta.columns]
for c in ANALYTES:
    meta[c] = pd.to_numeric(meta[c], errors="coerce")
print("  analytes: %d  (NULISA %d / PRM %d / ELISA %d)"
      % (len(ANALYTES), len(keep_nulisa), len(prm_targets), len(ELISA)))

b1 = meta[meta.timepoint.isin(["base", "1y"])]
piv = b1.pivot_table(index="participant_id", columns="timepoint",
                     values=ANALYTES + COG_H)
ar = pd.DataFrame(index=piv.index)
for c in ANALYTES + COG_H:
    if (c, "base") in piv.columns and (c, "1y") in piv.columns:
        ar["base_" + c] = piv[(c, "base")]
        ar["d_" + c] = piv[(c, "1y")] - piv[(c, "base")]

# ELISA concentrations are right-skewed; log2 puts them on the same footing as
# the other two platforms and normalises the residuals.
pv = b1.pivot_table(index="participant_id", columns="timepoint", values=ELISA)
add = pd.DataFrame(index=pv.index)
for c in ELISA:
    bb, yy = pv[(c, "base")].where(lambda s: s > 0), pv[(c, "1y")].where(lambda s: s > 0)
    add["base_log_" + c] = np.log2(bb)
    add["d_log_" + c] = np.log2(yy) - np.log2(bb)
    add["pct_" + c] = (yy - bb) / bb * 100

cov = (meta[meta.timepoint == "base"]
       .set_index("participant_id")[["arm", "age_base", "sex", "apoe4_carrier", "AD_profile"]])
cov["female"] = cov.sex.astype(str).str.upper().eq("F").astype(float)
cov.loc[cov.sex.isna(), "female"] = np.nan
cov["AD_pos"] = cov.AD_profile.eq(1).astype(float)
cov.loc[cov.AD_profile.isna(), "AD_pos"] = np.nan
lit = meta[meta.timepoint == "1y"].set_index("participant_id")["litemia"].rename("litemia_1y")

ar = (cov[["arm", "age_base", "female", "apoe4_carrier", "AD_pos"]]
      .join(ar).join(add).join(lit).reset_index())


def dcol(a):
    """Change-score column for an analyte, preferring the log2 ELISA version."""
    return ("d_log_" + a) if ("d_log_" + a) in ar.columns else ("d_" + a)


# Verify that the pipeline reproduces the previously merged values.
try:
    old = pd.read_excel(INPUT["master_old"])
    old.columns = [str(c).strip() for c in old.columns]
    chk = (meta[["csf_id", "NPTX2_NULISA"]]
           .merge(old[["Sample Number", "NPTX2_NULISA"]]
                  .rename(columns={"Sample Number": "csf_id",
                                   "NPTX2_NULISA": "old"}), on="csf_id").dropna())
    print("  reproduction check: n=%d, max abs difference = %.10f"
          % (len(chk), (chk.NPTX2_NULISA - chk.old).abs().max()))
except Exception as exc:
    print("  reproduction check skipped:", type(exc).__name__)


# ==============================================================================
# STEP 5  --  Cohort description and analysis sets
# ==============================================================================
print("\nSTEP 5  cohort description")

base = ar[ar.arm.notna()].copy()


def describe(v, kind, label):
    d = base[["arm", v]].dropna()
    L = d.loc[d.arm == "Lithium", v].astype(float)
    P = d.loc[d.arm == "Placebo", v].astype(float)
    if len(L) < 2 or len(P) < 2:
        return {"variable": label, "Lithium": "-", "n_Li": len(L),
                "Placebo": "-", "n_Pl": len(P), "test": "-", "P": np.nan}
    if kind == "num":
        p, test = stats.ttest_ind(L, P, equal_var=False).pvalue, "Welch t"
        fL = f"{L.mean():.1f} +/- {L.std():.1f}"
        fP = f"{P.mean():.1f} +/- {P.std():.1f}"
    else:
        tab = np.array([[(L == 1).sum(), (L == 0).sum()],
                        [(P == 1).sum(), (P == 0).sum()]])
        p, test = stats.fisher_exact(tab).pvalue, "Fisher"
        fL = f"{int((L == 1).sum())} ({100 * (L == 1).mean():.0f}%)"
        fP = f"{int((P == 1).sum())} ({100 * (P == 1).mean():.0f}%)"
    return {"variable": label, "Lithium": fL, "n_Li": len(L),
            "Placebo": fP, "n_Pl": len(P), "test": test, "P": round(p, 3)}


TABLE1_SPEC = [("age_base", "num", "Age, years"), ("female", "cat", "Female"),
               ("apoe4_carrier", "cat", "APOE e4 carrier"),
               ("AD_pos", "cat", "AD biomarker profile positive"),
               ("base_CDR_SB_h", "num", "CDR-SB"), ("base_ADAS_h", "num", "ADAS-Cog"),
               ("base_Ab_ELISA", "num", "CSF Ab42, pg/mL"),
               ("base_Tau_ELISA", "num", "CSF t-tau, pg/mL"),
               ("base_pTau_ELISA", "num", "CSF p-tau181, pg/mL"),
               ("base_log_NPTX2_ELISA", "num", "CSF NPTX2 (ELISA, log2)"),
               ("base_NPTX2_NULISA", "num", "CSF NPTX2 (NULISA)"),
               ("base_NPTX2_PRM", "num", "CSF NPTX2 (PRM-MS)")]
table1 = pd.DataFrame([describe(*s) for s in TABLE1_SPEC if s[0] in base.columns])
print(table1.to_string(index=False))

# Baseline balance across every analyte, to show whether any imbalance is
# systematic or consistent with chance.
rows = []
for a in ANALYTES:
    b = "base_" + a
    if b not in ar.columns:
        continue
    L = ar.loc[ar.arm == "Lithium", b].dropna()
    P = ar.loc[ar.arm == "Placebo", b].dropna()
    if len(L) < 5 or len(P) < 5:
        continue
    rows.append({"analyte": a, "platform": platform_of(a), "n_Li": len(L), "n_Pl": len(P),
                 "mean_Li": L.mean(), "mean_Pl": P.mean(),
                 "cohen_d": (L.mean() - P.mean()) / np.sqrt((L.var() + P.var()) / 2),
                 "p_welch": stats.ttest_ind(L, P, equal_var=False).pvalue,
                 "p_mwu": stats.mannwhitneyu(L, P).pvalue})
balance = pd.DataFrame(rows)
balance["q_platform"] = np.nan
for pl, g in balance.groupby("platform"):
    balance.loc[g.index, "q_platform"] = bh(g.p_welch)

# Analysis set sizes, which become the figure-legend numbers.
recs = []
for pl, col in [("NULISA", "NPTX2_NULISA"), ("PRM", "NPTX2_PRM"), ("ELISA", "NPTX2_ELISA")]:
    c = dcol(col)
    for (a, t), g in b1.groupby(["arm", "timepoint"]):
        recs.append({"platform": pl, "arm": a, "set": t, "n": g[col].notna().sum()})
    paired = ar.dropna(subset=[c]).groupby("arm").size()
    for a, v in paired.items():
        recs.append({"platform": pl, "arm": a, "set": "paired", "n": int(v)})
    for extra, lab in [("d_CDR_SB_h", "paired + CDR-SB"), ("litemia_1y", "paired + plasma Li")]:
        n2 = ar.dropna(subset=[c, extra]).groupby("arm").size()
        for a, v in n2.items():
            recs.append({"platform": pl, "arm": a, "set": lab, "n": int(v)})
n_by_analysis = pd.DataFrame(recs)


# ==============================================================================
# STEP 6  --  Correlations, age-adjusted partial correlations, FDR
# ==============================================================================
print("\nSTEP 6  correlations and multiple-testing correction")

D_ANALYTES = [dcol(a) for a in ANALYTES if dcol(a) in ar.columns]
OUTCOMES = ["litemia_1y"] + [c for c in ["d_" + x for x in COG_H] if c in ar.columns]


def residual(y, X):
    X = sm.add_constant(np.asarray(X, float), has_constant="add")
    return sm.OLS(np.asarray(y, float), X).fit().resid


rows = []
for arm in ["Lithium", "Placebo"]:
    A = ar[ar.arm == arm]
    for out in OUTCOMES:
        for c in D_ANALYTES:
            d = A[[out, c, "age_base"]].dropna()
            n = len(d)
            rec = {"arm": arm, "outcome": out, "analyte": c.replace("d_log_", "").replace("d_", ""),
                   "platform": platform_of(c), "n": n}
            if n >= 10:
                rec["r"], rec["p"] = stats.pearsonr(d[out], d[c])
                rec["rho"], rec["p_spearman"] = stats.spearmanr(d[out], d[c])
                try:
                    pr, _ = stats.pearsonr(residual(d[c], d[["age_base"]]),
                                           residual(d[out], d[["age_base"]]))
                    t = pr * np.sqrt((n - 3) / max(1e-12, 1 - pr ** 2))
                    rec["r_adjAge"], rec["p_adjAge"] = pr, 2 * stats.t.sf(abs(t), n - 3)
                except Exception:
                    rec["r_adjAge"] = rec["p_adjAge"] = np.nan
            else:
                for k in ["r", "p", "rho", "p_spearman", "r_adjAge", "p_adjAge"]:
                    rec[k] = np.nan
            rows.append(rec)
corr = pd.DataFrame(rows)
corr["q_pooled"] = corr["q_platform"] = np.nan
for (a, o), g in corr.groupby(["arm", "outcome"]):
    corr.loc[g.index, "q_pooled"] = bh(g.p)
    for pl, gp in g.groupby("platform"):
        corr.loc[gp.index, "q_platform"] = bh(gp.p)

summary = (corr.dropna(subset=["p"]).groupby(["arm", "outcome", "platform"])
           .agg(m=("p", "size"), nominal=("p", lambda s: (s < 0.05).sum()),
                q_pooled_05=("q_pooled", lambda s: (s < 0.05).sum()),
                q_platform_05=("q_platform", lambda s: (s < 0.05).sum()),
                q_platform_10=("q_platform", lambda s: (s < 0.10).sum()),
                min_q_platform=("q_platform", "min")).reset_index())
print(summary.to_string(index=False))


# ==============================================================================
# STEP 7  --  Within-arm change, between-arm ANCOVA, interaction
# ==============================================================================
print("\nSTEP 7  between-arm comparison and interaction")

rows = []
for c in D_ANALYTES:
    for arm in ["Lithium", "Placebo"]:
        s = ar.loc[ar.arm == arm, c].dropna()
        if len(s) < 8:
            continue
        rows.append({"analyte": c.replace("d_log_", "").replace("d_", ""),
                     "platform": platform_of(c), "arm": arm, "n": len(s),
                     "mean_change": s.mean(), "sd": s.std(),
                     "p_paired_t": stats.ttest_1samp(s, 0).pvalue,
                     "p_wilcoxon": stats.wilcoxon(s).pvalue})
paired_within = pd.DataFrame(rows)
paired_within["q_platform"] = np.nan
for (a, pl), g in paired_within.groupby(["arm", "platform"]):
    paired_within.loc[g.index, "q_platform"] = bh(g.p_paired_t)

COVS = ["age_base", "female", "apoe4_carrier"]


def ancova(dc, bc, covs=None):
    cols = ["arm", dc, bc] + (covs or [])
    s = ar[cols].dropna()
    if s.arm.nunique() < 2:
        return {}
    X = pd.concat([(s.arm == "Lithium").astype(float).rename("lithium"),
                   s[[bc] + (covs or [])].astype(float)], axis=1)
    m = sm.OLS(s[dc].astype(float), sm.add_constant(X)).fit()
    lo, hi = m.conf_int().loc["lithium"]
    return {"n": len(s), "beta": m.params["lithium"], "ci_lo": lo, "ci_hi": hi,
            "p": m.pvalues["lithium"]}


rows = []
for c in D_ANALYTES:
    bc = c.replace("d_log_", "base_log_").replace("d_", "base_", 1) \
        if c.startswith("d_log_") else "base_" + c[2:]
    if bc not in ar.columns:
        continue
    rec = {"analyte": c.replace("d_log_", "").replace("d_", ""), "platform": platform_of(c)}
    crude, adj = ancova(c, bc), ancova(c, bc, COVS)
    rec.update({f"{k}_crude": v for k, v in crude.items()})
    rec.update({f"{k}_adj": v for k, v in adj.items()})
    rows.append(rec)
arm_compare = pd.DataFrame(rows)
for src, dst in [("p_crude", "q_crude"), ("p_adj", "q_adj")]:
    if src in arm_compare:
        arm_compare[dst] = np.nan
        for pl, g in arm_compare.groupby("platform"):
            arm_compare.loc[g.index, dst] = bh(g[src])

rows = []
for out in [o for o in OUTCOMES if o != "litemia_1y"]:
    for c in D_ANALYTES:
        d = ar[["arm", out, c]].dropna()
        if d.arm.nunique() < 2 or len(d) < 15:
            continue
        g = (d.arm == "Lithium").astype(float)
        X = pd.DataFrame({"d": d[c].astype(float), "li": g,
                          "d_x_li": d[c].astype(float) * g})
        name = c.replace("d_log_", "").replace("d_", "")
        rec = {"outcome": out, "analyte": name, "platform": platform_of(c)}
        try:
            m = sm.OLS(d[out].astype(float), sm.add_constant(X)).fit()
            rec["beta_int"], rec["p_int"] = m.params["d_x_li"], m.pvalues["d_x_li"]
        except Exception:
            rec["beta_int"] = rec["p_int"] = np.nan
        sub = corr[(corr.outcome == out) & (corr.analyte == name)]
        rL = sub[sub.arm == "Lithium"]
        rP = sub[sub.arm == "Placebo"]
        if len(rL) and len(rP) and rL.r.notna().all() and rP.r.notna().all():
            r1, n1, r2, n2 = rL.r.iloc[0], rL.n.iloc[0], rP.r.iloc[0], rP.n.iloc[0]
            if n1 > 3 and n2 > 3:
                z = (np.arctanh(r1) - np.arctanh(r2)) / np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
                rec.update({"r_Li": r1, "r_Pl": r2, "fisher_z": z,
                            "p_fisher": 2 * stats.norm.sf(abs(z))})
        rows.append(rec)
interaction = pd.DataFrame(rows)
interaction["q_int_platform"] = np.nan
for (o, pl), g in interaction.groupby(["outcome", "platform"]):
    interaction.loc[g.index, "q_int_platform"] = bh(g.p_int)


# ==============================================================================
# STEP 8  --  Influence diagnostics
# ==============================================================================
print("\nSTEP 8  leave-one-out influence")


def loo_range(arm, out, col):
    d = ar[ar.arm == arm][[out, col]].dropna()
    full = stats.pearsonr(d[out], d[col])[0]
    rs = [stats.pearsonr(d.drop(i)[out], d.drop(i)[col])[0] for i in d.index]
    return {"arm": arm, "outcome": out, "analyte": col, "n": len(d),
            "r": full, "loo_min": min(rs), "loo_max": max(rs)}


loo = pd.DataFrame([loo_range(a, o, c) for a, o, c in
                    [("Lithium", "d_CDR_SB_h", "d_NPTX2_PRM"),
                     ("Lithium", "d_CDR_SB_h", "d_GRIA4_PRM"),
                     ("Lithium", "d_CDR_SB_h", "d_NPTX1_PRM"),
                     ("Lithium", "d_CDR_SB_h", "d_NPTXR_PRM"),
                     ("Placebo", "d_CDR_SB_h", "d_NPTX2_PRM")]
                    if c in ar.columns])
print(loo.round(3).to_string(index=False))


# ==============================================================================
# STEP 9  --  NPTX2 network:  correlation change and between-arm comparison
# ==============================================================================
print("\nSTEP 9  network analyses")


def rvec(M):
    """Correlation of column 0 with every other column, complete cases."""
    Z = (M - M.mean(0)) / M.std(0, ddof=1)
    return (Z[:, :1] * Z[:, 1:]).sum(0) / (len(M) - 1)


def build_block(platform, ref):
    cols = [a for a in ANALYTES if platform_of(a) == platform and a != ref]
    cols = [c for c in cols if ("base_" + c) in ar.columns and ("d_" + c) in ar.columns
            and ar["base_" + c].std() > 0]
    out = {}
    for arm in ["Lithium", "Placebo"]:
        A = ar[ar.arm == arm]
        B = A[["base_" + ref] + ["base_" + c for c in cols]].values
        Y = B + A[["d_" + ref] + ["d_" + c for c in cols]].values
        ok = ~np.isnan(np.hstack([B, Y])).any(axis=1)
        Bo, Yo = B[ok], Y[ok]
        keepc = [i for i in range(Bo.shape[1])
                 if Bo[:, i].std() > 0 and Yo[:, i].std() > 0]
        out[arm] = {"base": Bo[:, keepc], "y1": Yo[:, keepc], "n": int(ok.sum()),
                    "names": [([ref] + cols)[i] for i in keepc][1:]}
    return out


def edge_count(M):
    n = len(M)
    r = rvec(M)
    t = np.abs(r) * np.sqrt((n - 2) / np.clip(1 - r ** 2, 1e-12, None))
    return int((t > stats.t.ppf(0.975, n - 2)).sum()), float(
        np.abs(np.arctanh(np.clip(r, -0.999, 0.999))).mean())


A4, A5 = [], []
rng_global = np.random.default_rng(SEED)
for platform, ref in [("NULISA", "NPTX2_NULISA"), ("PRM", "NPTX2_PRM")]:
    D = build_block(platform, ref)
    names = D["Lithium"]["names"]
    for arm in ["Lithium", "Placebo"]:
        B, Y, n = D[arm]["base"], D[arm]["y1"], D[arm]["n"]
        rb, ry = rvec(B), rvec(Y)
        rng = np.random.default_rng(SEED)
        boot = np.array([rvec(Y[i]) - rvec(B[i])
                         for i in (rng.integers(0, n, (N_BOOT, n)))])
        p_boot = 2 * np.minimum((boot <= 0).mean(0), (boot >= 0).mean(0))
        A4.append(pd.DataFrame({"platform": platform, "arm": arm,
                                "analyte": D[arm]["names"], "n": n,
                                "r_base": rb, "r_1y": ry, "delta_r": ry - rb,
                                "ci_lo": np.percentile(boot, 2.5, axis=0),
                                "ci_hi": np.percentile(boot, 97.5, axis=0),
                                "p_boot": p_boot,
                                "q_boot": bh(np.clip(p_boot, 1e-4, 1))}))

    obs = {}
    for arm in ["Lithium", "Placebo"]:
        eb, zb = edge_count(D[arm]["base"])
        ey, zy = edge_count(D[arm]["y1"])
        obs[arm] = {"edges_base": eb, "edges_1y": ey, "d_edges": ey - eb,
                    "meanz_base": zb, "meanz_1y": zy, "d_meanz": zy - zb,
                    "n": D[arm]["n"]}

    nmin = min(D["Lithium"]["n"], D["Placebo"]["n"])
    rng = np.random.default_rng(SEED + 1)
    ds = {arm: np.array([edge_count(D[arm]["y1"][rng.choice(D[arm]["n"], nmin, False)])[0]
                         - edge_count(D[arm]["base"][rng.choice(D[arm]["n"], nmin, False)])[0]
                         for _ in range(N_DOWNSAMPLE)]) for arm in ["Lithium", "Placebo"]}
    p_ds = 2 * min((ds["Lithium"] <= ds["Placebo"]).mean(),
                   (ds["Lithium"] >= ds["Placebo"]).mean())

    allB = np.vstack([D["Lithium"]["base"], D["Placebo"]["base"]])
    allY = np.vstack([D["Lithium"]["y1"], D["Placebo"]["y1"]])
    nL, N = D["Lithium"]["n"], len(allB)
    obs_e = obs["Lithium"]["d_edges"] - obs["Placebo"]["d_edges"]
    obs_z = obs["Lithium"]["d_meanz"] - obs["Placebo"]["d_meanz"]
    rng = np.random.default_rng(SEED + 2)
    ne, nz = [], []
    for _ in range(N_PERM):
        p_ = rng.permutation(N)
        a, b = p_[:nL], p_[nL:]
        ea, za = edge_count(allY[a])[0] - edge_count(allB[a])[0], \
                 edge_count(allY[a])[1] - edge_count(allB[a])[1]
        eb_, zb_ = edge_count(allY[b])[0] - edge_count(allB[b])[0], \
                   edge_count(allY[b])[1] - edge_count(allB[b])[1]
        ne.append(ea - eb_)
        nz.append(za - zb_)
    A5.append({"platform": platform,
               "obs_diff_edges": obs_e,
               "p_perm_edges": float((np.abs(ne) >= abs(obs_e)).mean()),
               "obs_diff_meanz": obs_z,
               "p_perm_meanz": float((np.abs(nz) >= abs(obs_z)).mean()),
               "p_downsample": p_ds,
               **{f"{k}_{a}": v for a in obs for k, v in obs[a].items()}})

A4 = pd.concat(A4, ignore_index=True)
A5 = pd.DataFrame(A5)
print(A5.round(4).to_string(index=False))


# ==============================================================================
# STEP 10  --  Reconciliation with the originally reported network analysis
# ==============================================================================
print("\nSTEP 10  reconciliation of the original network result")

try:
    orig_list = [re.sub(r"_NULISA$", "", c) for c in old.columns if c.endswith("_NULISA")]
    W = nul_wide.reset_index().merge(arm_map[["participant_id", "arm"]],
                                     on="participant_id", how="left").dropna(subset=["arm"])

    def edges_pairwise(arm, tp, analytes, paired_only, complete_only):
        d = W[(W.arm == arm) & (W.timepoint == tp)]
        if paired_only:
            pb = set(W[(W.arm == arm) & (W.timepoint == "base")].participant_id)
            py = set(W[(W.arm == arm) & (W.timepoint == "1y")].participant_id)
            d = d[d.participant_id.isin(pb & py)]
        M = d[["NPTX2"] + [a for a in analytes if a in W.columns and a != "NPTX2"]].astype(float)
        if complete_only:
            M = M.dropna()
        res = []
        for c in M.columns[1:]:
            s = M[[M.columns[0], c]].dropna()
            if len(s) < 4 or s.iloc[:, 0].std() == 0 or s.iloc[:, 1].std() == 0:
                res.append(np.nan)
                continue
            res.append(stats.pearsonr(s.iloc[:, 0], s.iloc[:, 1])[1])
        return pd.Series(res, index=M.columns[1:]), len(M)

    ladder = []
    for label, analytes, paired_only, complete_only in [
            ("0. original conditions", orig_list, False, False),
            ("1. + paired participants only", orig_list, True, False),
            ("2. + QC analyte set", keep_nulisa, True, False),
            ("3. + complete cases only", keep_nulisa, True, True)]:
        for arm in ["Lithium", "Placebo"]:
            pb, nb = edges_pairwise(arm, "base", analytes, paired_only, complete_only)
            py, ny = edges_pairwise(arm, "1y", analytes, paired_only, complete_only)
            sb, sy = (pb < 0.05).fillna(False), (py < 0.05).fillna(False)
            gained, lost = int((~sb & sy).sum()), int((sb & ~sy).sum())
            tot = gained + lost
            ladder.append({"condition": label, "arm": arm, "n_base": nb, "n_1y": ny,
                           "base_sig": int(sb.sum()), "y1_sig": int(sy.sum()),
                           "gained": gained, "lost": lost,
                           "net": int(sy.sum()) - int(sb.sum()),
                           "McNemar_p": 1.0 if tot == 0 else
                           stats.binomtest(gained, tot, 0.5).pvalue})
    reconciliation = pd.DataFrame(ladder)
    print(reconciliation.round(4).to_string(index=False))
except Exception as exc:
    reconciliation = pd.DataFrame()
    print("  reconciliation skipped:", type(exc).__name__, exc)


# ==============================================================================
# STEP 11  --  Figures
# ==============================================================================
print("\nSTEP 11  figures")


def scatter_panel(axx, xcol, ycol, arms, xlab, ylab, title=None, legend=False):
    lines = []
    for arm in arms:
        s = ar[ar.arm == arm][[xcol, ycol]].dropna()
        if len(s) < 4:
            continue
        st = ARM_STYLE[arm]
        axx.scatter(s[xcol], s[ycol], s=22, facecolor=st["c"], edgecolor="white",
                    linewidth=0.6, marker=st["m"], zorder=3,
                    label=f"{arm} (n={len(s)})")
        r, p = stats.pearsonr(s[xcol], s[ycol])
        b, a0 = np.polyfit(s[xcol], s[ycol], 1)
        xs = np.linspace(s[xcol].min(), s[xcol].max(), 50)
        axx.plot(xs, a0 + b * xs, color=st["c"], lw=1.6, zorder=2)
        lines.append((st["c"], f"r = {r:+.3f}, P = {p:.3f}"))
    for i, (c, t) in enumerate(lines):
        axx.text(0.03, 0.96 - i * 0.085, t, transform=axx.transAxes,
                 fontsize=6.6, color=c, va="top")
    axx.set_xlabel(xlab)
    axx.set_ylabel(ylab)
    if title:
        axx.set_title(title, fontsize=8.5)
    axx.spines[["top", "right"]].set_visible(False)
    axx.grid(alpha=0.18, lw=0.5)
    if legend:
        axx.legend(fontsize=6.5, frameon=False, loc="best")


# Figure 1B-I  --  plasma lithium versus analyte change, lithium arm
EXPOSURE = [("NPTX2_ELISA", "CSF NPTX2 (ELISA, log2)"), ("pTau_ELISA", "CSF p-tau181 (ELISA)"),
            ("NEGR1_PRM", "CSF NEGR1 (PRM-MS)"), ("CRH_NULISA", "CSF CRH (NULISA)"),
            ("S100B_NULISA", "CSF S100B (NULISA)"), ("MME_NULISA", "CSF MME (NULISA)"),
            ("POSTN_NULISA", "CSF POSTN (NULISA)")]
fig, axes = plt.subplots(2, 4, figsize=(11.2, 5.4))
for axx, (a, lab) in zip(axes.ravel(), EXPOSURE):
    c = dcol(a)
    if c not in ar.columns:
        axx.axis("off")
        continue
    scatter_panel(axx, "litemia_1y", c, ["Lithium"],
                  "Plasma lithium at 1 year (mEq/L)", f"Change in {lab}")
    q = corr[(corr.arm == "Lithium") & (corr.outcome == "litemia_1y") &
             (corr.analyte == a)]["q_platform"]
    if len(q) and pd.notna(q.iloc[0]):
        axx.text(0.03, 0.80, f"q = {q.iloc[0]:.3f}", transform=axx.transAxes,
                 fontsize=6.6, color="#555", va="top")
axes.ravel()[-1].axis("off")
plt.tight_layout()
save_fig(fig, "Figure1B-I_plasma_lithium_panels")
plt.close(fig)

# Figure 1K  --  NPTX2 change versus CDR-SB change, stratified by arm
fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.4))
for axx, (a, lab) in zip(axes, [("NPTX2_ELISA", "ELISA (log2)"),
                                ("NPTX2_NULISA", "NULISA (NPQ)"),
                                ("NPTX2_PRM", "PRM-MS (log2)")]):
    scatter_panel(axx, dcol(a), "d_CDR_SB_h", ["Lithium", "Placebo"],
                  f"Change in CSF NPTX2, {lab}", "Change in CDR-SB",
                  title=lab, legend=(a == "NPTX2_ELISA"))
plt.tight_layout()
save_fig(fig, "Figure1K_NPTX2_vs_CDRSB_by_arm")
plt.close(fig)

# Figure 2A  --  the pre-specified NPTX2 axis
AXIS = [("NPTX2_PRM", "NPTX2"), ("NPTX1_PRM", "NPTX1"),
        ("NPTXR_PRM", "NPTXR"), ("GRIA4_PRM", "GluA4 (GRIA4)")]
AXIS = [(a, l) for a, l in AXIS if dcol(a) in ar.columns]
fig, axes = plt.subplots(1, len(AXIS), figsize=(3.25 * len(AXIS), 3.3))
for axx, (a, lab) in zip(np.atleast_1d(axes), AXIS):
    scatter_panel(axx, dcol(a), "d_CDR_SB_h", ["Lithium", "Placebo"],
                  f"Change in CSF {lab} (PRM-MS)", "Change in CDR-SB",
                  title=lab, legend=(a == AXIS[0][0]))
plt.tight_layout()
save_fig(fig, "Figure2A_NPTX_axis_scatter")
plt.close(fig)

# Figure 2B  --  forest plot with the between-arm comparison
rows = []
for a, lab in AXIS:
    c = dcol(a)
    rec = {"label": lab, "analyte": a}
    for arm in ["Lithium", "Placebo"]:
        s = ar[ar.arm == arm][[c, "d_CDR_SB_h"]].dropna()
        r, p = stats.pearsonr(s[c], s["d_CDR_SB_h"])
        lo, hi = fisher_ci(r, len(s))
        rec.update({f"r_{arm}": r, f"lo_{arm}": lo, f"hi_{arm}": hi,
                    f"n_{arm}": len(s), f"p_{arm}": p})
    z = (np.arctanh(rec["r_Lithium"]) - np.arctanh(rec["r_Placebo"])) / \
        np.sqrt(1 / (rec["n_Lithium"] - 3) + 1 / (rec["n_Placebo"] - 3))
    rec["p_fisher"] = 2 * stats.norm.sf(abs(z))
    ii = interaction[(interaction.outcome == "d_CDR_SB_h") & (interaction.analyte == a)]
    rec["p_int"] = ii.p_int.iloc[0] if len(ii) else np.nan
    rec["q_int"] = ii.q_int_platform.iloc[0] if len(ii) else np.nan
    rows.append(rec)
F2 = pd.DataFrame(rows)

fig, axx = plt.subplots(figsize=(7.4, 3.2))
yb = np.arange(len(F2))[::-1]
for arm in ["Lithium", "Placebo"]:
    off = 0.17 if arm == "Lithium" else -0.17
    st = ARM_STYLE[arm]
    axx.errorbar(F2[f"r_{arm}"], yb + off,
                 xerr=[F2[f"r_{arm}"] - F2[f"lo_{arm}"], F2[f"hi_{arm}"] - F2[f"r_{arm}"]],
                 fmt=st["m"], ms=6, color=st["c"], ecolor=st["c"], elinewidth=1.5,
                 capsize=2.5, mec="white", mew=0.6, zorder=3,
                 label=f"{arm} (n={int(F2[f'n_{arm}'].iloc[0])})")
axx.axvline(0, color="#888", lw=0.8, ls="--")
axx.set_yticks(yb)
axx.set_yticklabels(F2.label, fontsize=8)
axx.set_xlabel("Pearson r (change in protein vs change in CDR-SB), 95% CI")
axx.set_xlim(-0.95, 0.95)
for i, row in F2.iterrows():
    star = "**" if row.p_fisher < 0.01 else ("*" if row.p_fisher < 0.05 else "ns")
    axx.text(0.90, yb[i], f"between-arm P = {row.p_fisher:.3f} {star}",
             fontsize=6.4, ha="right", va="center", color="#444")
axx.spines[["top", "right"]].set_visible(False)
axx.grid(axis="x", alpha=0.18, lw=0.5)
axx.legend(fontsize=7, frameon=False, loc="lower left")
plt.tight_layout()
save_fig(fig, "Figure2B_NPTX_axis_forest")
plt.close(fig)


def volcano(outcome, arm, title, fname, ntop=8):
    d = corr[(corr.arm == arm) & (corr.outcome == outcome)].dropna(subset=["p"]).copy()
    d["nlp"] = -np.log10(d.p)
    fig, a = plt.subplots(figsize=(5.4, 4.4))
    for pl, g in d.groupby("platform"):
        for sub, col in [(g[g.q_platform >= 0.10], C_NULL),
                         (g[(g.p < 0.05) & (g.q_platform >= 0.10)], C_NOM),
                         (g[g.q_platform < 0.10], C_FDR)]:
            if len(sub):
                a.scatter(sub.r, sub.nlp, s=26, marker=PLATFORM_MARKER[pl],
                          facecolor=col, edgecolor="white", linewidth=0.5, zorder=3)
    a.axhline(-np.log10(0.05), ls="--", lw=0.8, c="#999")
    a.axvline(0, lw=0.6, c="#ccc")
    for _, row in d.nsmallest(ntop, "p").iterrows():
        a.annotate(row.analyte.rsplit("_", 1)[0], (row.r, row.nlp), fontsize=6.2,
                   xytext=(4, 3), textcoords="offset points", color="#333")
    handles = [Line2D([], [], marker=PLATFORM_MARKER[k], ls="", mfc="#888",
                      mec="white", ms=6, label=k) for k in PLATFORM_MARKER] + \
              [Line2D([], [], marker="o", ls="", mfc=c, mec="white", ms=6, label=l)
               for c, l in [(C_NULL, "ns"), (C_NOM, "P < 0.05"), (C_FDR, "q < 0.10")]]
    a.legend(handles=handles, fontsize=6.2, frameon=False, ncol=2, loc="upper left")
    a.set_xlabel("Pearson r")
    a.set_ylabel(r"$-\log_{10}$ P")
    a.set_title(title, fontsize=8.5)
    a.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    save_fig(fig, fname)
    plt.close(fig)


volcano("litemia_1y", "Lithium",
        "Plasma lithium vs CSF analyte change (lithium arm)",
        "Figure1J_volcano_plasma_lithium")
volcano("d_CDR_SB_h", "Lithium",
        "CDR-SB change vs CSF analyte change (lithium arm)",
        "Figure1J_volcano_CDRSB_lithium")
volcano("d_CDR_SB_h", "Placebo",
        "CDR-SB change vs CSF analyte change (placebo arm)",
        "S_Figure_volcano_CDRSB_placebo")

# Participant flow, computed from the data rather than hard-coded
flow = [("Archived CSF from the randomized trial\n(low-dose lithium in amnestic MCI)",
         f"n = {(ar.arm == 'Lithium').sum()}", f"n = {(ar.arm == 'Placebo').sum()}"),
        ("Paired baseline and 1-year CSF\nELISA and PRM-MS",
         f"n = {ar[ar.arm == 'Lithium']['d_NPTX2_PRM'].notna().sum()}",
         f"n = {ar[ar.arm == 'Placebo']['d_NPTX2_PRM'].notna().sum()}"),
        ("Paired baseline and 1-year CSF\nNULISA (1 sample not returned by vendor)",
         f"n = {ar[ar.arm == 'Lithium']['d_NPTX2_NULISA'].notna().sum()}",
         f"n = {ar[ar.arm == 'Placebo']['d_NPTX2_NULISA'].notna().sum()}"),
        ("NPTX2 change and CDR-SB change\n(association with cognition)",
         f"n = {ar[(ar.arm == 'Lithium')][['d_NPTX2_PRM', 'd_CDR_SB_h']].notna().all(1).sum()}",
         f"n = {ar[(ar.arm == 'Placebo')][['d_NPTX2_PRM', 'd_CDR_SB_h']].notna().all(1).sum()}"),
        ("NPTX2 change and plasma lithium\n(exposure-response, lithium arm only)",
         f"n = {ar[(ar.arm == 'Lithium')][['d_NPTX2_PRM', 'litemia_1y']].notna().all(1).sum()}",
         "not applicable")]
fig, axx = plt.subplots(figsize=(7.5, 8.6))
axx.axis("off")
axx.set_xlim(0, 10)
axx.set_ylim(0, len(flow) * 2 + 0.6)
for i, (label, nl, np_) in enumerate(flow):
    y = (len(flow) - i) * 2 - 1.2
    axx.add_patch(FancyBboxPatch((0.25, y), 5.6, 1.15, boxstyle="round,pad=0.10",
                                 fc="#F2F5FB", ec="#2E75B6", lw=1.0))
    axx.text(3.05, y + 0.575, label, ha="center", va="center", fontsize=8.4)
    for j, (t, c) in enumerate([(nl, "#2E75B6"), (np_, "#C55A11")]):
        axx.add_patch(FancyBboxPatch((6.25 + j * 1.85, y + 0.22), 1.65, 0.72,
                                     boxstyle="round,pad=0.06", fc="white", ec=c, lw=1.0))
        axx.text(7.08 + j * 1.85, y + 0.58, t, ha="center", va="center",
                 fontsize=7.6, color=c)
    if i < len(flow) - 1:
        axx.annotate("", xy=(3.05, y - 0.22), xytext=(3.05, y),
                     arrowprops=dict(arrowstyle="-|>", color="#2E75B6", lw=1.0))
axx.text(7.08, len(flow) * 2 - 0.22, "Lithium", ha="center", fontsize=8.6,
         weight="bold", color="#2E75B6")
axx.text(8.93, len(flow) * 2 - 0.22, "Placebo", ha="center", fontsize=8.6,
         weight="bold", color="#C55A11")
plt.tight_layout()
save_fig(fig, "S_Figure_participant_flow")
plt.close(fig)


# ==============================================================================
# STEP 12  --  Write every table
# ==============================================================================
print("\nSTEP 12  writing tables")

save_table(meta, "meta_master")
save_table(ar, "analysis_ready")
save_table(table1, "Table1_baseline_characteristics")
save_table({"analyte_QC": qc,
            "sample_QC": sample_det.rename("detectability").reset_index(),
            "n_by_analysis": n_by_analysis,
            "baseline_balance": balance},
           "S_Table_QC_and_sample_disposition")
save_table({"summary": summary, "correlations": corr,
            "arm_comparison": arm_compare, "interaction": interaction,
            "within_arm_change": paired_within, "leave_one_out": loo},
           "S_Table_association_results")
save_table({"A4_correlation_change": A4, "A5_network_comparison": A5,
            "reconciliation_ladder": reconciliation},
           "S_Table_network_analysis")
save_table({"figure2_values": F2}, "S_Table_figure2_values")

# ==============================================================================
# STEP 13  --  Regression-to-the-mean diagnostics for the CSF NPTX2 decline
# ==============================================================================
# Baseline CSF NPTX2 is higher in the lithium arm, so a reviewer may ask whether
# the larger one-year decline in that arm is simply regression to the mean.
# Four independent diagnostics are reported.
print("\nSTEP 13  regression-to-the-mean diagnostics")

NPTX2_COLS = [("ELISA", "NPTX2_ELISA"), ("NULISA", "NPTX2_NULISA"), ("PRM", "NPTX2_PRM")]
NPTX2_COLS = [(p, c) for p, c in NPTX2_COLS if dcol(c) in ar.columns]

# (13a) Within each arm, is the baseline value correlated with the subsequent
#       change?  Under pure regression to the mean this correlation is strongly
#       negative.  A correlation near zero argues against it.
rtm_within = []
for pl, c in NPTX2_COLS:
    d_c = dcol(c)
    b_c = d_c.replace("d_log_", "base_log_") if d_c.startswith("d_log_") else "base_" + c
    if b_c not in ar.columns:
        continue
    for arm in ["Lithium", "Placebo"]:
        d = ar.loc[ar.arm == arm, [b_c, d_c]].dropna()
        r, p = (stats.pearsonr(d[b_c], d[d_c]) if len(d) >= 10 else (np.nan, np.nan))
        rtm_within.append({"platform": pl, "arm": arm, "n": len(d),
                           "baseline_col": b_c, "change_col": d_c,
                           "r_baseline_vs_change": r, "p": p})
rtm_within = pd.DataFrame(rtm_within)
print(rtm_within.round(4).to_string(index=False))

# (13b) Reliability of the baseline measurement itself.  Regression to the mean
#       is driven by measurement error, so a baseline that reproduces across
#       three independent platforms leaves little room for it.
rel = []
for i in range(len(NPTX2_COLS)):
    for j in range(i + 1, len(NPTX2_COLS)):
        (p1, c1), (p2, c2) = NPTX2_COLS[i], NPTX2_COLS[j]
        b1c = ("base_log_" + c1) if ("base_log_" + c1) in ar.columns else ("base_" + c1)
        b2c = ("base_log_" + c2) if ("base_log_" + c2) in ar.columns else ("base_" + c2)
        if b1c not in ar.columns or b2c not in ar.columns:
            continue
        d = ar[[b1c, b2c]].dropna()
        r, p = (stats.pearsonr(d[b1c], d[b2c]) if len(d) >= 10 else (np.nan, np.nan))
        rel.append({"pair": "%s vs %s" % (p1, p2), "n": len(d), "r": r, "p": p})
rel = pd.DataFrame(rel)
print(rel.round(4).to_string(index=False))

# (13c) ANCOVA on the change score adjusting for a composite baseline averaged
#       across platforms.  A composite carries less measurement error than any
#       single platform, so if the effect is real it should survive, and if it
#       is regression to the mean it should shrink.
zs = lambda s: (s - s.mean()) / s.std(ddof=1)
bcols = [(("base_log_" + c) if ("base_log_" + c) in ar.columns else ("base_" + c))
         for _, c in NPTX2_COLS]
bcols = [c for c in bcols if c in ar.columns]
ar["base_NPTX2_composite"] = pd.concat([zs(ar[c]) for c in bcols], axis=1).mean(axis=1)

rtm_ancova = []
for pl, c in NPTX2_COLS:
    d_c = dcol(c)
    b_c = d_c.replace("d_log_", "base_log_") if d_c.startswith("d_log_") else "base_" + c
    for adj_lab, adj in [("own baseline", b_c), ("composite baseline", "base_NPTX2_composite")]:
        if adj not in ar.columns:
            continue
        d = ar[[d_c, adj, "arm"]].dropna().copy()
        d["li"] = (d.arm == "Lithium").astype(int)
        if len(d) < 12:
            continue
        X = sm.add_constant(d[["li", adj]].astype(float).values, has_constant="add")
        m = sm.OLS(d[d_c].astype(float).values, X).fit()
        lo, hi = m.conf_int()[1]
        rtm_ancova.append({"platform": pl, "adjustment": adj_lab, "n": len(d),
                           "beta_lithium": m.params[1], "ci_low": lo,
                           "ci_high": hi, "p": m.pvalues[1]})
rtm_ancova = pd.DataFrame(rtm_ancova)
print(rtm_ancova.round(4).to_string(index=False))

# (13d) Batch check.  If treatment arms were unevenly distributed across NULISA
#       plates, a plate effect could masquerade as a treatment effect.
plate_tab = pd.crosstab(meta.PlateID, meta.arm)
try:
    chi2, p_plate, dof, _ = stats.chi2_contingency(plate_tab.values)
except Exception:
    chi2 = p_plate = dof = np.nan
print("  plate x arm chi-square: chi2 = %.3f, df = %s, P = %.3f" % (chi2, dof, p_plate))
plate_out = plate_tab.reset_index()
plate_out["chi2"] = chi2
plate_out["df"] = dof
plate_out["p"] = p_plate

save_table({"baseline_vs_change": rtm_within,
            "baseline_reliability": rel,
            "ancova_own_vs_composite": rtm_ancova,
            "plate_by_arm": plate_out},
           "S_Table_RTM_diagnostics")


# ==============================================================================
# STEP 14  --  Gene Ontology enrichment of the cognition-associated analytes
# ==============================================================================
# The measured panel is not a random sample of the genome: it is composed
# largely of synaptic and neuronal proteins.  Testing against a whole-genome
# background would therefore return synaptic terms whatever the result, so the
# primary analysis uses the measured panel as the background.  The whole-genome
# version is reported alongside it for comparison only.
print("\nSTEP 14  Gene Ontology enrichment")

# Analyte labels that are not official gene symbols.
SYMBOL_FIX = {"GluA4": "GRIA4", "GluA1": "GRIA1", "GluA2": "GRIA2", "GluA3": "GRIA3",
              "pTau": "MAPT", "Tau": "MAPT", "Ab": "APP", "aSyn": "SNCA",
              "NfL": "NEFL", "pTDP43": "TARDBP", "TDP43": "TARDBP",
              "GFAP": "GFAP", "Abeta40": "APP", "Abeta42": "APP",
              "IL1B": "IL1B", "CD200R1": "CD200R1", "14-3-3": "YWHAZ",
              "YWHAZ": "YWHAZ", "S100B": "S100B", "PDGFRB": "PDGFRB"}


def gene_symbol(analyte):
    """Strip the platform suffix from an analyte label and return a gene symbol."""
    base = re.sub(r"_(NULISA|PRM|ELISA)$", "", str(analyte))
    return SYMBOL_FIX.get(base, base).upper()


CDR_OUTCOME = "d_CDR_SB_h"
go_src = corr[(corr.arm == "Lithium") & (corr.outcome == CDR_OUTCOME)].dropna(subset=["p"])
gl_input = sorted({gene_symbol(a) for a in go_src.loc[go_src.p < 0.05, "analyte"]})
gl_bg = sorted({gene_symbol(a) for a in go_src["analyte"]})
print("  input genes: %d   panel background genes: %d" % (len(gl_input), len(gl_bg)))
print("  input:", ", ".join(gl_input))

GO_LIBS = ["GO_Biological_Process_2023", "GO_Cellular_Component_2023",
           "GO_Molecular_Function_2023"]


def fisher_enrich(genes, background, library_dict, min_overlap=2):
    """Fisher exact enrichment of `genes` within `background` for one library."""
    genes, background = set(genes) & set(background), set(background)
    N = len(background)
    out = []
    for term, members in library_dict.items():
        ann = set(m.upper() for m in members) & background
        k = len(ann & genes)
        if k < min_overlap or not ann:
            continue
        table = [[k, len(genes) - k], [len(ann) - k, N - len(genes) - len(ann) + k]]
        odds, p = stats.fisher_exact(table, alternative="greater")
        out.append({"term": term, "overlap": k, "n_input": len(genes),
                    "n_term_in_background": len(ann), "n_background": N,
                    "odds_ratio": odds, "p": p,
                    "overlap_genes": ";".join(sorted(ann & genes))})
    out = pd.DataFrame(out)
    if len(out):
        out["q_BH"] = bh(out.p)
        out = out.sort_values("p").reset_index(drop=True)
    return out


go_tables, go_ok = {}, False
try:
    import gseapy
    for lib in GO_LIBS:
        lib_dict = gseapy.get_library(name=lib, organism="Human")
        panel = fisher_enrich(gl_input, gl_bg, lib_dict)
        panel.insert(0, "background", "measured panel (%d genes)" % len(gl_bg))
        genome_bg = sorted({g.upper() for v in lib_dict.values() for g in v})
        genome = fisher_enrich(gl_input, genome_bg, lib_dict)
        genome.insert(0, "background", "whole genome (%d genes)" % len(genome_bg))
        go_tables[lib.replace("GO_", "").replace("_2023", "")[:31]] = \
            pd.concat([panel, genome], ignore_index=True)
        print("  %s: %d terms (panel background), top P = %s"
              % (lib, len(panel), "%.4g" % panel.p.iloc[0] if len(panel) else "na"))
    go_ok = True
except Exception as exc:
    print("  GO enrichment skipped:", type(exc).__name__, exc)
    print("  (gseapy needs internet access to Enrichr; run this step where that is allowed)")

if go_ok:
    go_tables["input_genes"] = pd.DataFrame({"gene": gl_input})
    go_tables["background_genes"] = pd.DataFrame({"gene": gl_bg})
    save_table(go_tables, "S_Table_GO_enrichment")

    # Figure S3:  top terms under the measured-panel background.
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.6))
    show = [k for k in ["Cellular_Component", "Biological_Process"] if k in go_tables]
    for axx, key in zip(axes, show):
        t = go_tables[key]
        t = t[t.background.str.startswith("measured")].nsmallest(8, "p").iloc[::-1]
        lab = [re.sub(r"\s*\(GO:\d+\)$", "", s) for s in t.term]
        lab = [s if len(s) <= 46 else s[:43] + "..." for s in lab]
        axx.barh(range(len(t)), -np.log10(t.p), color="#2E75B6", height=0.62)
        axx.set_yticks(range(len(t)))
        axx.set_yticklabels(lab, fontsize=7)
        for i, (k, m) in enumerate(zip(t.overlap, t.n_term_in_background)):
            axx.text(-np.log10(t.p.iloc[i]) + 0.05, i, "%d/%d" % (k, m),
                     va="center", fontsize=6.2, color="#444")
        axx.axvline(-np.log10(0.05), ls="--", lw=0.8, c="#999")
        axx.set_xlabel(r"$-\log_{10}$ P (uncorrected)")
        axx.set_title(key.replace("_", " "), fontsize=8.5)
        axx.spines[["top", "right"]].set_visible(False)
        axx.set_xlim(0, max(2.2, -np.log10(t.p.min()) * 1.35))
    plt.tight_layout()
    save_fig(fig, "FigureS3_GO_enrichment")
    plt.close(fig)

# ==============================================================================
# STEP 15  --  Low-QC sample sensitivity analysis (NULISA)
# ==============================================================================
# One NULISA sample fell just below the vendor's 70% sample-detectability
# threshold for CSF.  It was retained because excluding it would have broken a
# longitudinal pair.  This step repeats the NULISA NPTX2 analyses with that
# participant removed, so the effect of the decision can be inspected directly.
print("\nSTEP 15  low-QC sample sensitivity analysis")

low_samples = sample_det[sample_det < SAMPLE_QC_MIN]
if len(low_samples) == 0:
    print("  no sample falls below the %.0f%% threshold; nothing to test."
          % (SAMPLE_QC_MIN * 100))
else:
    low_name = low_samples.index[0]
    low_row = samples.loc[samples.SampleName == low_name].iloc[0]
    low_pid = low_row.participant_id
    low_arm = arm_map.set_index("participant_id").arm.get(low_pid, "unknown")
    print("  low-QC sample: %s  (participant %s, %s arm, %s, detectability %.1f%%)"
          % (low_name, low_pid, low_arm, low_row.timepoint, low_samples.iloc[0] * 100))

    D_NP = dcol("NPTX2_NULISA")
    B_NP = ("base_log_NPTX2_NULISA" if "base_log_NPTX2_NULISA" in ar.columns
            else "base_NPTX2_NULISA")

    def nulisa_panel(tag, A):
        rec = {"analysis_set": tag}
        for arm in ["Lithium", "Placebo"]:
            s = A.loc[A.arm == arm, D_NP].dropna()
            rec["n_" + arm] = len(s)
            rec["mean_change_" + arm] = s.mean() if len(s) else np.nan
            rec["P_paired_t_" + arm] = stats.ttest_1samp(s, 0).pvalue if len(s) > 2 else np.nan
            rec["P_wilcoxon_" + arm] = stats.wilcoxon(s).pvalue if len(s) > 2 else np.nan
        L = A[A.arm == "Lithium"]
        for lab, out in [("CDRSB", "d_CDR_SB_h"), ("plasma_Li", "litemia_1y")]:
            if out not in A.columns:
                continue
            d = L[[D_NP, out]].dropna()
            rec["n_" + lab] = len(d)
            if len(d) >= 10:
                rec["r_" + lab], rec["P_" + lab] = stats.pearsonr(d[D_NP], d[out])
            else:
                rec["r_" + lab] = rec["P_" + lab] = np.nan
        if B_NP in A.columns:
            d = A[[D_NP, B_NP, "arm"]].dropna()
            X = sm.add_constant(pd.concat(
                [(d.arm == "Lithium").astype(float).rename("li"),
                 d[B_NP].astype(float)], axis=1), has_constant="add")
            m = sm.OLS(d[D_NP].astype(float), X).fit()
            lo, hi = m.conf_int().loc["li"]
            rec.update({"n_ANCOVA": len(d), "beta_ANCOVA": m.params["li"],
                        "ci_lo_ANCOVA": lo, "ci_hi_ANCOVA": hi,
                        "P_ANCOVA": m.pvalues["li"]})
        return rec

    sens = pd.DataFrame([nulisa_panel("all samples (primary)", ar),
                         nulisa_panel("excluding participant %s" % low_pid,
                                      ar[ar.participant_id != low_pid])])
    cols = [c for c in sens.columns if c != "analysis_set"]
    print(sens.set_index("analysis_set")[cols].T.round(4).to_string())

    # Plate composition of the longitudinal pairs, reported alongside, because a
    # plate-by-arm test alone does not address plate effects on within-person
    # change.  Pairs split across plates carry plate variation as noise rather
    # than as a systematic between-arm difference.
    pl = (nul_wide.reset_index()[["SampleName", "participant_id", "PlateID"]]
          .dropna(subset=["participant_id"]).drop_duplicates())
    n_samp = pl.groupby("participant_id").SampleName.nunique()
    n_plate = pl.groupby("participant_id").PlateID.nunique()
    complete = n_samp[n_samp == 2].index
    split = int((n_plate[complete] > 1).sum())
    plate_pairs = pd.DataFrame([{
        "participants_with_NULISA_data": int(len(n_samp)),
        "complete_longitudinal_pairs": int(len(complete)),
        "pairs_split_across_plates": split,
        "pairs_on_one_plate": int(len(complete)) - split}])
    print("\n  complete NULISA pairs: %d, of which %d have their two timepoints "
          "on different plates" % (len(complete), split))

    save_table({"nulisa_lowQC_sensitivity": sens,
                "plate_composition_of_pairs": plate_pairs,
                "sample_detectability": sample_det.rename("detectability").reset_index()},
               "S_Table_lowQC_sensitivity")

print("\nDone.  Tables in %s, figures in %s" % (OUT_T, OUT_F))
