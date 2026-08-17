#!/usr/bin/env python3
"""
02b_collect_seeds_extra.py — 扩充 PHB 降解种子序列
来源：
  1) PAZy API（substrate=PHB / PHA，含 FASTA）
  2) ESTHER/PhaDED 家族页成员（Esterase_phb_PHAZ、PHAZ7_phb_depolymerase）
  3) UniProt unreviewed 关键物种酶（protein_name 检索，补齐文献酶）
输出（合并 02_collect_seeds.py 的 seeds_curated.faa）：
  data/seeds/seeds_all.faa         全量
  data/seeds/seeds_family.faa      按家族整理（header 含 family 标签）
  data/seeds/seeds_family_manifest.tsv
"""
import argparse
import csv
import os
import re
import sys
import time

import requests

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
PAZY_BASE = "https://api.pazy.eu/api"


def fetch_pazy(substrate: str) -> list[dict]:
    """PAZy 蛋白列表（含序列元数据）"""
    url = f"{PAZY_BASE}/proteins/"
    out = []
    params = {"substrate": substrate, "page_size": "100"}
    try:
        while url:
            r = requests.get(url, params=params, timeout=60)
            r.raise_for_status()
            data = r.json()
            if isinstance(data, list):
                out.extend(data)
                break
            out.extend(data.get("results", data.get("items", [])))
            url = data.get("next")
            params = None
        return out
    except Exception as e:
        print(f"  [PAZy {substrate}] ERR: {e}", file=sys.stderr)
        return []


def fetch_pazy_fasta(substrate: str) -> str:
    try:
        r = requests.get(f"{PAZY_BASE}/proteins/fasta/", params={"substrate": substrate}, timeout=120)
        r.raise_for_status()
        return r.text
    except Exception as e:
        print(f"  [PAZy fasta {substrate}] ERR: {e}", file=sys.stderr)
        return ""


def fetch_esther_members(family: str) -> list[str]:
    """从 ESTHER 家族页提取成员 accession"""
    try:
        r = requests.get(f"https://bioweb.supagro.inrae.fr/ESTHER/family/{family}", timeout=30)
        r.raise_for_status()
        # UniProt/ESTHER accession 模式（A2VBKQ 等）+ 页面可能含 SWISS-ACC 链接
        accs = set(re.findall(r"(?:uniprot\.org/(?:uniprot|entry)/|ACCESSION[^\n]{0,20}|>)([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][0-9]{2}){1,2}[0-9])", r.text))
        # 兜底：常见 6-10 字符 accession
        accs |= set(re.findall(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9])\b", r.text))
        accs |= set(re.findall(r"\b([A-NR-Z][0-9](?:[A-Z][0-9]{2}){1,2}[0-9])\b", r.text))
        return sorted(accs)
    except Exception as e:
        print(f"  [ESTHER {family}] ERR: {e}", file=sys.stderr)
        return []


def fetch_uniprot(query: str, max_results: int = 100) -> list[dict]:
    params = {"query": query, "format": "json", "size": str(min(max_results, 500))}
    try:
        r = requests.get(UNIPROT_STREAM, params=params, timeout=60)
        r.raise_for_status()
        return r.json().get("results", [])
    except Exception as e:
        print(f"  [UniProt {query[:40]}] ERR: {e}", file=sys.stderr)
        return []


def hit_to_fasta(h: dict, family: str) -> tuple[str, dict]:
    acc = h.get("primaryAccession", "")
    seq = (h.get("sequence") or {}).get("value", "")
    pdesc = h.get("proteinDescription") or {}
    ecs = [e.get("value", "") for e in pdesc.get("ecNumbers") or []]
    pname = (pdesc.get("recommendedName") or {}).get("fullName", {}).get("value", "")
    org = (h.get("organism") or {}).get("scientificName", "")
    genes = []
    for g in h.get("genes") or []:
        if g.get("geneName"):
            genes.append(g["geneName"].get("value", ""))
    gene_str = ";".join(dict.fromkeys(genes))
    reviewed = "true" if "reviewed" in str(h.get("entryType", "")) else "false"
    header = f">{acc}|{family}|{org}|{gene_str}|EC:{';'.join(ecs)}|rev:{reviewed}|{pname[:60]}"
    meta = {"accession": acc, "family": family, "organism": org, "gene": gene_str,
            "ec": ";".join(ecs), "reviewed": reviewed, "protein_name": pname}
    return header, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/seeds")
    ap.add_argument("--existing", default="data/seeds/seeds_curated.faa")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    fasta = []
    manifest = []
    seen = set()

    def add(header, meta, seq):
        acc = meta["accession"]
        if not acc or acc in seen or not seq:
            return
        seen.add(acc)
        fasta.append(header)
        fasta.append(seq)
        manifest.append(meta)

    # ---------- 1. PAZy ----------
    print("[*] PAZy PHB / PHA / PCL substrates")
    for sub in ["PHB", "PHA", "PCL", "3HB"]:
        entries = fetch_pazy(sub)
        print(f"    {sub}: {len(entries)} entries")
        for e in entries:
            acc = e.get("protein_id") or e.get("uniprot_id") or e.get("accession") or e.get("id", "")
            seq = e.get("sequence", "")
            name = e.get("protein_name") or e.get("name") or ""
            org = e.get("organism") or ""
            header = f">{acc}|PAZy_{sub}|{org}|||rev:false|{name}"
            meta = {"accession": str(acc), "family": f"PAZy_{sub}", "organism": org,
                    "gene": "", "ec": "", "reviewed": "false", "protein_name": str(name)}
            add(header, meta, seq)
        # FASTA 直接抓取（如果列表无序列）
        if not any(e.get("sequence") for e in entries):
            fa = fetch_pazy_fasta(sub)
            for h, s in re.findall(r"^>(.*?)\n([^>]+)", fa, re.M | re.S):
                acc = h.split("|")[0].split()[0]
                header = f">{acc}|PAZy_{sub}|{h[:80]}"
                meta = {"accession": acc, "family": f"PAZy_{sub}", "organism": h[:60],
                        "gene": "", "ec": "", "reviewed": "false", "protein_name": h[:60]}
                add(header, meta, s.strip())
        time.sleep(0.5)

    # ---------- 2. ESTHER 家族 ----------
    print("[*] ESTHER/PhaDED families")
    esther_families = {
        "Esterase_phb_PHAZ": "ePhaZ_ESTHER",
        "PHAZ7_phb_depolymerase": "iPhaZ_PHAZ7",
    }
    for fam, tag in esther_families.items():
        accs = fetch_esther_members(fam)
        print(f"    {fam}: {len(accs)} member accessions")
        for i in range(0, len(accs), 50):
            batch = accs[i:i + 50]
            q = " OR ".join(f"accession:{a}" for a in batch)
            hits = fetch_uniprot(q)
            for h in hits:
                header, meta = hit_to_fasta(h, tag)
                seq = (h.get("sequence") or {}).get("value", "")
                add(header, meta, seq)
        time.sleep(0.5)

    # ---------- 3. UniProt unreviewed 关键物种 ----------
    print("[*] UniProt key species (incl. unreviewed)")
    species_queries = [
        ("Rhodospirillum rubrum", "iPhaZ"),
        ("Paucimonas lemoignei", "ePhaZ"),
        ("Comamonas acidovorans", "ePhaZ"),
        ("Thermus thermophilus", "ePhaZ"),
        ("Bacillus thuringiensis", "iPhaZ"),
        ("Sinorhizobium meliloti", "iPhaZ"),
        ("Streptomyces ascomycinicus", "ePhaZ"),
        ("Ralstonia pickettii", "ePhaZ"),
        ("Pseudomonas lemoignei", "ePhaZ"),
        ("Cupriavidus necator", "iPhaZ"),
    ]
    for org, fam in species_queries:
        q = f'(protein_name:"depolymerase" OR protein_name:"PHA depolymerase" OR protein_name:"poly(3-hydroxybutyrate) depolymerase") AND organism_name:"{org}"'
        hits = fetch_uniprot(q)
        print(f"    {org}: {len(hits)} hits")
        for h in hits:
            header, meta = hit_to_fasta(h, fam)
            seq = (h.get("sequence") or {}).get("value", "")
            add(header, meta, seq)
        time.sleep(0.3)

    # ---------- 合并已有 curated ----------
    if os.path.exists(args.existing):
        cur = open(args.existing).read().strip()
        for block in cur.split("\n>"):
            if not block.strip():
                continue
            lines = block.split("\n")
            h = lines[0].lstrip(">")
            seq = "".join(lines[1:]).strip()
            acc = h.split("|")[0]
            if acc and acc not in seen and seq:
                seen.add(acc)
                fasta.append(">" + h)
                fasta.append(seq)
                manifest.append({"accession": acc, "family": h.split("|")[1] if len(h.split("|")) > 1 else "unclassified",
                                 "organism": "", "gene": "", "ec": "", "reviewed": "true", "protein_name": h})

    # ---------- 输出 ----------
    allfa = os.path.join(args.outdir, "seeds_all.faa")
    with open(allfa, "w", encoding="utf-8") as f:
        f.write("\n".join(fasta) + "\n")

    famfa = os.path.join(args.outdir, "seeds_family.faa")
    man = os.path.join(args.outdir, "seeds_family_manifest.tsv")
    with open(famfa, "w", encoding="utf-8") as f:
        f.write("\n".join(fasta) + "\n")
    with open(man, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "family", "organism", "gene", "ec", "reviewed", "protein_name"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(manifest)

    # 家族统计
    from collections import Counter
    cnt = Counter(m["family"] for m in manifest)
    print("\n[DONE] total seeds:", len(seen))
    for fam, n in cnt.most_common():
        print(f"  {fam}: {n}")
    print("  FASTA:", allfa)
    print("  manifest:", man)


if __name__ == "__main__":
    main()
