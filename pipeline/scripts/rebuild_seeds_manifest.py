#!/usr/bin/env python3
"""rebuild_seeds_manifest.py — 用 UniProt 实时查询重建种子 manifest（修复 reviewed 字段）

背景：旧 manifest 的 `reviewed` 字段因 `"reviewed" in entryType` 子串误判（"unreviewed"
也含 "reviewed"）而不可信。本脚本对现有 78 条 accession 逐条查询 UniProt：
  - entryType（"UniProtKB reviewed (Swiss-Prot)" vs "... unreviewed (TrEMBL)"）→ reviewed 真值
  - references.citationCrossReferences → 证据 PMID / DOI
  - retrieval_date（查询日期）+ split（train/validation，80/20 按 accession md5 取模）

输入: pipeline/seeds/seeds_manifest.tsv（旧版）
输出: pipeline/seeds/seeds_manifest.tsv（覆盖，新增 evidence/retrieval_date/split 列）
      pipeline/seeds/seeds_stats.json（更新）

依赖: 仅标准库（urllib）。本地或服务器均可运行。
用法: python pipeline/scripts/rebuild_seeds_manifest.py [--manifest pipeline/seeds/seeds_manifest.tsv]
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request

BASE = "https://rest.uniprot.org/uniprotkb"


def get_entry(acc: str, retries: int = 3):
    url = f"{BASE}/{acc}"
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None  # 已删除条目
            if attempt == retries - 1:
                print(f"  [WARN] {acc}: HTTP {e.code}", file=sys.stderr)
                return None
            time.sleep(1 + attempt)
        except Exception as e:
            if attempt == retries - 1:
                print(f"  [WARN] {acc}: {e}", file=sys.stderr)
                return None
            time.sleep(1 + attempt)
    return None


def is_reviewed(entry) -> bool:
    et = str(entry.get("entryType", "") or "")
    return et.startswith("UniProtKB reviewed")


def extract_evidence(entry) -> str:
    pmids, dois = set(), set()
    for ref in entry.get("references") or []:
        for xr in (ref.get("citation") or {}).get("citationCrossReferences") or []:
            db = (xr.get("database") or "").lower()
            xid = xr.get("id") or ""
            if db == "pubmed" and xid:
                pmids.add(xid)
            elif db == "doi" and xid:
                dois.add(xid)
    parts = []
    if pmids:
        parts.append("pmid:" + ";".join(sorted(pmids)))
    if dois:
        parts.append("doi:" + ";".join(sorted(dois)))
    return ";".join(parts)


def assign_split(accession: str) -> str:
    h = int(hashlib.md5(accession.encode("utf-8")).hexdigest(), 16)
    return "train" if h % 10 < 8 else "validation"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="pipeline/seeds/seeds_manifest.tsv")
    ap.add_argument("--outdir", default="pipeline/seeds")
    args = ap.parse_args()

    rows = []
    with open(args.manifest, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    print(f"现有 manifest: {len(rows)} 条")

    today = dt.date.today().isoformat()
    out_rows = []
    n_reviewed = 0
    for i, r in enumerate(rows):
        acc = r["accession"]
        entry = get_entry(acc)
        if entry is None:
            reviewed = r.get("reviewed", "false")
            evidence = ""
        else:
            reviewed = "true" if is_reviewed(entry) else "false"
            evidence = extract_evidence(entry)
        if reviewed == "true":
            n_reviewed += 1
        out_rows.append({
            "accession": acc,
            "query_group": r.get("query_group", ""),
            "reviewed": reviewed,
            "organism": r.get("organism", ""),
            "gene": r.get("gene", ""),
            "ec": r.get("ec", ""),
            "protein_name": r.get("protein_name", ""),
            "length": r.get("length", ""),
            "lineage": r.get("lineage", ""),
            "evidence": evidence,
            "retrieval_date": today,
            "split": assign_split(acc),
        })
        if (i + 1) % 10 == 0:
            print(f"  ... {i + 1}/{len(rows)}")
        time.sleep(0.3)  # 礼貌限速

    out_path = os.path.join(args.outdir, "seeds_manifest.tsv")
    fields = ["accession", "query_group", "reviewed", "organism", "gene", "ec",
              "protein_name", "length", "lineage", "evidence", "retrieval_date", "split"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)

    n_train = sum(1 for r in out_rows if r["split"] == "train")
    n_val = sum(1 for r in out_rows if r["split"] == "validation")
    stats = {
        "total": len(out_rows),
        "reviewed": n_reviewed,
        "unreviewed": len(out_rows) - n_reviewed,
        "train": n_train,
        "validation": n_val,
        "retrieval_date": today,
    }
    with open(os.path.join(args.outdir, "seeds_stats.json"), "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
        f.write("\n")

    print(f"\n[DONE] reviewed={n_reviewed}, unreviewed={len(out_rows) - n_reviewed}, "
          f"train/val={n_train}/{n_val}")
    print(f"  manifest -> {out_path}")
    print(f"  stats -> {os.path.join(args.outdir, 'seeds_stats.json')}")


if __name__ == "__main__":
    main()
