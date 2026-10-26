"""Extract Tables S2 and S3 of Kroenke et al. (2016) from the supplement PDF.

Input:  empirical/data/raw/kroenke2016/kroenke2016_supplement.pdf, the
        Supplemental Digital Content of Psychosomatic Medicine 78(6), 716-727
        (PMC4927366); not redistributed, see empirical/README.md.
Output: empirical/kroenke2016/table_s2_correlations.csv (the 16 x 16 item
        correlations, all trials pooled) and table_s3_loadings.csv
        (standardized loadings per trial: one-factor, two-factor, bifactor
        general factor).

Run from the repository root:  python3 -I empirical/kroenke2016/extract_tables.py
Requires pdftotext (poppler).
"""
import csv
import pathlib
import re
import subprocess

PDF = pathlib.Path("empirical/data/raw/kroenke2016/kroenke2016_supplement.pdf")
OUT = pathlib.Path("empirical/kroenke2016")
ITEMS = [f"phq{i}" for i in range(1, 10)] + [f"gad{i}" for i in range(1, 8)]
TRIALS = ["SCOPE", "ESCAPE", "INCPAD"]

text = subprocess.run(["pdftotext", "-layout", str(PDF), "-"], check=True,
                      capture_output=True, text=True).stdout
s2 = text[text.index("Table S2."):text.index("Table S3.")]
s3 = text[text.index("Table S3."):]
num = re.compile(r"-?\d*\.\d+")

# Table S2: one row per item, 15 correlations and "--" on the diagonal.
rows = {}
for line in s2.splitlines():
    parts = line.split()
    if parts and parts[0] in ITEMS and len(parts) == 17:
        rows[parts[0]] = [1.0 if v == "--" else float(v) for v in parts[1:]]
assert list(rows) == ITEMS, "Table S2: expected 16 item rows"
# The printed matrix is not exactly symmetric. Keep the lower triangle (row
# item later in the list) and record every disagreement.
asym = []
for i, a in enumerate(ITEMS):
    for j in range(i):
        b = ITEMS[j]
        if rows[a][j] != rows[b][i]:
            asym.append((b, a, rows[b][i], rows[a][j]))
        rows[b][i] = rows[a][j]
with open(OUT / "table_s2_correlations.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["item"] + ITEMS)
    for a in ITEMS:
        w.writerow([a] + [f"{v:.2f}" for v in rows[a]])
with open(OUT / "table_s2_asymmetries.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["item_1", "item_2", "upper_triangle", "lower_triangle_used"])
    w.writerows(asym)
print(f"Table S2: {len(asym)} asymmetric pair(s); lower triangle used")

# Table S3: per item, nine loadings in the order SCOPE (1F, 2F, bifactor g),
# ESCAPE (same), INCPAD (same). The source spells some labels "pqh-9".
loadings = []
for line in s3.splitlines():
    m = re.match(r"\s*(phq|pqh|gad)-(9|7)\s+(\d)\s+(.*)", line)
    if not m:
        continue
    vals = [float(v) for v in num.findall(m.group(4))]
    assert len(vals) == 9, f"Table S3: expected 9 loadings in: {line!r}"
    item = ("phq" if m.group(1) in ("phq", "pqh") else "gad") + m.group(3)
    for k, trial in enumerate(TRIALS):
        one, two, g = vals[3 * k:3 * k + 3]
        loadings.append([trial, item, "depression" if item.startswith("phq") else "anxiety",
                         one, two, g])
assert len(loadings) == 48, "Table S3: expected 16 items x 3 trials"
with open(OUT / "table_s3_loadings.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["trial", "item", "block", "one_factor", "two_factor", "bifactor_general"])
    w.writerows(loadings)
print("wrote", OUT / "table_s2_correlations.csv", "and", OUT / "table_s3_loadings.csv")
