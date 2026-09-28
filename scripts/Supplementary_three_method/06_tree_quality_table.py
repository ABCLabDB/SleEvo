#!/usr/bin/env python3
"""
Supplementary Table 13 — gene-tree quality and cross-criterion robustness.

Merges, per phenotype:
  * MP tree quality   <- v2/tables/SupplementaryTable_MP_tree_statistics_v2.tsv   (median over 143 genes)
  * ML tree quality   <- ML_based/tables/SupplementaryTable_ML_tree_statistics.tsv (median over 143 genes)
  * ML node support   <- ML_based/trees/<tag>/*_ML.nwk, parsing the "SH-aLRT/UFBoot" node labels
  * reproducibility   <- 3method_comparison/tables/SupplementaryTable_reproducibility_by_phenotype.tsv

Run 05_three_method_comparison.R first: it writes the reproducibility table this script reads.
"""
import csv, os, re, statistics as st

MP  = "/disk2/bijsy/Test/Sleep_association_Test/Result/Maximum_Parsimony_for_Cladistics"
V2  = os.path.join(MP, "v2", "tables")
MLT = os.path.join(MP, "ML_based", "tables")
CMP = os.path.join(MP, "3method_comparison", "tables")

def rd(p): return list(csv.DictReader(open(p, encoding="utf-8-sig"), delimiter="\t"))
def fl(v):
    try:
        x = float(v); return x if x == x else None
    except Exception:
        return None

def ml_support(tag):
    """Mean SH-aLRT / UFBoot over every internal node of every ML tree of this phenotype."""
    d = os.path.join(MP, "ML_based", "trees", tag)
    alrt, boot, per = [], [], []
    for f in sorted(os.listdir(d)):
        if not f.endswith("_ML.nwk"):
            continue
        t = open(os.path.join(d, f)).read()
        m = re.findall(r"\)([\d.]+)/([\d.]+):", t)      # node label = SH-aLRT/UFBoot
        if not m:
            continue
        a = [float(x) for x, _ in m]; b = [float(y) for _, y in m]
        alrt += a; boot += b
        per.append(sum(1 for x, y in zip(a, b) if x >= 80 and y >= 95) / len(a) * 100)
    return dict(
        SHaLRT_mean  = round(st.mean(alrt), 1),
        UFBoot_mean  = round(st.mean(boot), 1),
        pct_UFBoot95 = round(sum(1 for x in boot if x >= 95) / len(boot) * 100, 1),
        pct_both     = round(st.mean(per), 1))

PH = [("Sleep duration", "Sleep_duration"), ("NREM ratio", "NREM_ratio"),
      ("Sleep timing", "Sleep_timing"),     ("Sleep frequency", "Sleep_frequency")]

mpstat = rd(os.path.join(V2,  "SupplementaryTable_MP_tree_statistics_v2.tsv"))
mlstat = rd(os.path.join(MLT, "SupplementaryTable_ML_tree_statistics.tsv"))
rep    = {r["Phenotype"]: r for r in rd(os.path.join(CMP, "SupplementaryTable_reproducibility_by_phenotype.tsv"))}

rows = []
for lab, tag in PH:
    M = [r for r in mpstat if r["Phenotype"] == tag]
    L = [r for r in mlstat if r["Phenotype"] == tag]
    md  = lambda S, k: st.median([fl(r[k]) for r in S if fl(r[k]) is not None])
    sup = ml_support(tag)
    R   = rep[lab]
    rows.append(dict(
        Phenotype = lab, N_genes = 143, N_species = int(md(M, "N_tips")),
        # ---- MP tree quality (median across the 143 gene trees) ----
        MP_parsimony_score            = int(md(M, "MP_parsimony_score")),
        MP_CI                         = round(md(M, "CI"), 3),
        MP_RI                         = round(md(M, "RI"), 3),
        MP_RC                         = round(md(M, "RC"), 3),
        MP_bootstrap_mean             = round(md(M, "Bootstrap_mean"), 1),
        MP_pct_nodes_BS70             = round(md(M, "Pct_nodes_bs70"), 1),
        MP_equally_parsimonious_trees = int(md(M, "N_best_trees_recovered")),
        MP_strict_consensus_resolution= round(md(M, "Strict_consensus_resolution"), 1),
        # ---- ML tree quality ----
        ML_parsimony_informative_sites = int(md(L, "N_parsimony_informative")),
        ML_logL                        = round(md(L, "logL"), 0),
        ML_SHaLRT_mean                 = sup["SHaLRT_mean"],
        ML_UFBoot_mean                 = sup["UFBoot_mean"],
        ML_pct_nodes_UFBoot95          = sup["pct_UFBoot95"],
        ML_pct_nodes_SHaLRT80_UFBoot95 = sup["pct_both"],
        # ---- cross-criterion robustness ----
        Criterion              = R["Criterion"],
        NJ_significant         = int(R["NJ_significant"]),
        MP_significant         = int(R["MP_significant"]),
        ML_significant         = int(R["ML_significant"]),
        Retained_in_MP         = int(R["Replicated_in_MP"]),
        MP_retention_pct       = float(R["MP_retention_pct"]),
        Retained_in_ML         = int(R["Replicated_in_ML"]),
        ML_retention_pct       = float(R["ML_retention_pct"]),
        Supported_by_all_three = int(R["All_three"]),
        Min_possible_overlap   = int(R["Min_possible_overlap_NJ_MP"]),
        Spearman_NJ_MP = float(R["Spearman_NJ_MP"]),
        Spearman_NJ_ML = float(R["Spearman_NJ_ML"]),
        Spearman_MP_ML = float(R["Spearman_MP_ML"])))

out = os.path.join(CMP, "SupplementaryTable_tree_quality_and_reproducibility.tsv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
print("saved:", out)
for r in rows:
    print(f"  {r['Phenotype']:<16} n={r['N_species']:<3} PS={r['MP_parsimony_score']:<6} "
          f"CI={r['MP_CI']:<6} RI={r['MP_RI']:<6} UFBoot={r['ML_UFBoot_mean']}%")
