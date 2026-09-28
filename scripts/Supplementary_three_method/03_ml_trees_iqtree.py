#!/usr/bin/env python3
"""ML 계통수 (IQ-TREE) — 표현형별 143 유전자 x 4 = 572 트리

기존 두 기준과 동일 조건:
  · 같은 MUSCLE 정렬, 같은 표현형별 종 부분집합
  · 표현형 필터 -> 트리 생성 순서

IQ-TREE 설정
  -m MFP        ModelFinder Plus (유전자별 최적 치환모델 선택)
  -B 1000       ultrafast bootstrap
  --alrt 1000   SH-aLRT
  -T 1          유전자당 1스레드, 유전자 간 병렬 (572개이므로 이게 효율적)
"""
import os, csv, subprocess, shutil, tempfile, time, json
from concurrent.futures import ProcessPoolExecutor

ROOT = "/disk2/bijsy/Test/Sleep_association_Test"
ML   = os.path.join(ROOT, "Result/Maximum_Parsimony_for_Cladistics/ML_based")
IQ   = os.path.join(ML, "conda_iqtree/bin/iqtree3")
DATA = "/disk2/bijsy/Test/Muscle"
META = os.path.join(ROOT, "METADATA_fixedVersion.txt")
GENES= os.path.join(ROOT, "supplementaryS2.tsv")
TR   = os.path.join(ML, "trees"); RAW = os.path.join(ML, "iqtree_raw")
CK   = os.path.join(ML, "checkpoint")
for d in (TR, RAW, CK, os.path.join(ML,"tables")): os.makedirs(d, exist_ok=True)

NPROC = int(os.environ.get("N_CORES", "128"))
UFBOOT = int(os.environ.get("UFBOOT", "1000"))
ALRT   = int(os.environ.get("ALRT", "1000"))
SEED   = 319
ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

PH = [("Sleep_duration",  "Total_sleep_time_per_day",        "num"),
      ("NREM_ratio",      "Percentage_of_NREM_time_per_day", "num"),
      ("Sleep_timing",    "_Sleep_timing_per_day",           "chr"),
      ("Sleep_frequency", "Number_of_sleep_times_per_day",   "chr")]

def norm_sp(x): return x.replace(" ", "_").replace("'", "_")

meta = list(csv.DictReader(open(META, encoding="utf-8-sig"), delimiter="\t"))
for r in meta:
    r["sym"] = r["Species_symbol_name_ensembl"].strip().lower().replace(" ", "_")
    r["sp"]  = norm_sp(r["Species_name_ensembl"])

def species_of(col, typ):
    out = {}
    for r in meta:
        v = (r.get(col) or "").strip()
        if typ == "num":
            try: float(v)
            except: continue
        else:
            if v in ("", "NA"): continue
        out[r["sym"]] = r["sp"]
    return out

def read_fa(p):
    rec=[]; n=None; b=[]
    for l in open(p):
        l=l.rstrip("\n")
        if l.startswith(">"):
            if n: rec.append((n,"".join(b)))
            n=l[1:].strip(); b=[]
        elif n: b.append(l)
    if n: rec.append((n,"".join(b)))
    return rec

def prep(gene, keep):
    src = os.path.join(DATA, f"{gene}_muscle.fasta")
    if not os.path.exists(src): raise FileNotFoundError("FASTA 없음")
    out=[]; seen=set()
    for h,s in read_fa(src):
        parts = h.split(":")
        sym = (parts[1] if len(parts)>=2 else parts[0]).strip().lower().replace(" ","_")
        if sym in keep and sym not in seen:
            seen.add(sym); out.append((keep[sym], s))
    if len(out) < 5: raise ValueError(f"종 {len(out)} < 5")
    return out

def one(task):
    gene, tag, keep = task
    ck = os.path.join(CK, f"{tag}__{gene}.json")
    if os.path.exists(ck):
        try: return json.load(open(ck))
        except Exception: pass
    t0 = time.time()
    wd = tempfile.mkdtemp(prefix=f"iq_{tag}_{gene}_")
    try:
        recs = prep(gene, keep)
        fa = os.path.join(wd, "aln.fasta")
        with open(fa,"w") as fh:
            for s,q in recs: fh.write(f">{s}\n{q}\n")
        pre = os.path.join(wd, gene)
        cmd = [IQ, "-s", fa, "-m", "MFP", "-B", str(UFBOOT), "--alrt", str(ALRT),
               "-T", "1", "--prefix", pre, "--seed", str(SEED), "--quiet", "-redo"]
        r = subprocess.run(cmd, capture_output=True, env=ENV, timeout=36000)
        tf = pre + ".treefile"
        if not os.path.exists(tf):
            return dict(status="ERROR", Phenotype=tag, Gene=gene,
                        Error=(r.stderr.decode()[-200:] or f"rc={r.returncode}"))
        os.makedirs(os.path.join(TR, tag), exist_ok=True)
        shutil.copy(tf, os.path.join(TR, tag, f"{gene}_ML.nwk"))
        kd = os.path.join(RAW, tag); os.makedirs(kd, exist_ok=True)
        for ext in (".iqtree", ".log", ".contree"):
            if os.path.exists(pre+ext): shutil.copy(pre+ext, os.path.join(kd, gene+ext))
        # .iqtree 파싱
        info = dict(model=None, lnL=None, nsite=None, npars_inf=None,
                    nconst=None, ndistinct=None, AIC=None, BIC=None)
        txt = open(pre+".iqtree", errors="ignore").read()
        import re
        def grab(pat, cast=str):
            m = re.search(pat, txt)
            return cast(m.group(1)) if m else None
        info["model"]     = grab(r"Best-fit model according to BIC:\s*(\S+)")
        info["lnL"]       = grab(r"Log-likelihood of the tree:\s*(-?[\d.]+)", float)
        info["nsite"]     = grab(r"Input data:\s*\d+\s+sequences with\s+(\d+)\s+nucleotide", int)
        info["nconst"]    = grab(r"Number of constant sites:\s*(\d+)", int)
        info["npars_inf"] = grab(r"Number of parsimony informative sites:\s*(\d+)", int)
        info["ndistinct"] = grab(r"Number of distinct site patterns:\s*(\d+)", int)
        info["AIC"]       = grab(r"Akaike information criterion \(AIC\) score:\s*([\d.]+)", float)
        info["BIC"]       = grab(r"Bayesian information criterion \(BIC\) score:\s*([\d.]+)", float)
        res = dict(status="OK", Phenotype=tag, Gene=gene, N_tips=len(recs),
                   Alignment_length=len(recs[0][1]),
                   Best_model=info["model"], logL=info["lnL"],
                   N_sites=info["nsite"], N_constant_sites=info["nconst"],
                   N_parsimony_informative=info["npars_inf"],
                   N_distinct_patterns=info["ndistinct"],
                   AIC=info["AIC"], BIC=info["BIC"],
                   seconds=round(time.time()-t0,1), Error=None)
        json.dump(res, open(ck,"w"))
        print(f"[OK] {tag:<16} {gene:<12} tips={res['N_tips']:>2} "
              f"model={res['Best_model']} lnL={res['logL']} ({res['seconds']:.0f}s)", flush=True)
        return res
    except Exception as e:
        res = dict(status="ERROR", Phenotype=tag, Gene=gene, Error=f"{type(e).__name__}: {e}")
        json.dump(res, open(ck,"w"))
        print(f"[ERR] {tag:<16} {gene:<12} {res['Error'][:80]}", flush=True)
        return res
    finally:
        shutil.rmtree(wd, ignore_errors=True)

if __name__ == "__main__":
    gl = list(csv.DictReader(open(GENES, encoding="utf-8-sig"), delimiter="\t"))
    genes = sorted({r["Gene_symbol"] for r in gl
                    if os.path.exists(os.path.join(DATA, r["Gene_symbol"]+"_muscle.fasta"))})
    tasks = []
    for tag, col, typ in PH:
        keep = species_of(col, typ)
        print(f"  {tag:<16} {len(keep)}종")
        tasks += [(g, tag, keep) for g in genes]
    print(f"유전자 {len(genes)} | 트리 {len(tasks)} | 코어 {NPROC} | UFBoot {UFBOOT} | SH-aLRT {ALRT}\n", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=NPROC) as ex:
        res = list(ex.map(one, tasks))
    ok = [r for r in res if r.get("status")=="OK"]
    print(f"\n완료 {len(ok)} / {len(res)} | {(time.time()-t0)/60:.1f}분", flush=True)

    cols = ["Phenotype","Gene","N_tips","Alignment_length","N_sites","N_constant_sites",
            "N_parsimony_informative","N_distinct_patterns","Best_model","logL","AIC","BIC","seconds"]
    with open(os.path.join(ML,"tables","SupplementaryTable_ML_tree_statistics.tsv"),"w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t",extrasaction="ignore")
        w.writeheader(); w.writerows(ok)
    for tag,_,_ in PH:
        sub=[r for r in ok if r["Phenotype"]==tag]
        if sub:
            with open(os.path.join(ML,"tables",f"ML_tree_statistics_{tag}.tsv"),"w",newline="") as fh:
                w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t",extrasaction="ignore")
                w.writeheader(); w.writerows(sub)
    bad=[r for r in res if r.get("status")!="OK"]
    if bad:
        with open(os.path.join(ML,"tables","ML_tree_errors.tsv"),"w",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=["Phenotype","Gene","Error"],delimiter="\t",extrasaction="ignore")
            w.writeheader(); w.writerows(bad)
    from collections import Counter
    print("\n=== 표현형별 ===")
    for tag,_,_ in PH:
        sub=[r for r in ok if r["Phenotype"]==tag]
        if not sub: continue
        mods=Counter(r["Best_model"] for r in sub).most_common(3)
        print(f"  {tag:<16} {len(sub):>3}개 | 최빈 모델 {mods}")
    print("[DONE]")
