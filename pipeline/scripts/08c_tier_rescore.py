#!/usr/bin/env python3
"""08c_tier_rescore.py — 用 curated 金标准 HMM 对命中做三级重评分（tier1/tier2/tier3）

Python 版，替代 bash 版 08c_tier_rescore.sh：cwd 相对路径、无引号陷阱、fail-closed。
  tier1: curated HMM E<1e-20（严格）
  tier2: curated HMM E<1e-10（中等）
  tier3: 现有 validated（宽模型 + 通用验证）

流程（每个核心家族）：
  validated.faa → hmmsearch(tier2, E<1e-10) → tier2.ids + tier2.faa
  tier2.faa    → hmmsearch(tier1, E<1e-20) → tier1.ids + tier1.faa（从 validated.faa 提取）

用法（服务器 T141，cwd=工作区根）:
  ~/miniconda3/envs/phb_gtdb/bin/python scripts/08c_tier_rescore.py
"""
import os
import shutil
import subprocess
import sys
import argparse

HMMSEARCH = shutil.which("hmmsearch") or os.path.expanduser(
    "~/miniconda3/envs/phb_gtdb/bin/hmmsearch")

SEQDIR = "data/screen/family_seqs"
TIERDIR = "data/screen/tiers"
CURATED = {
    "ePhaZ": "data/hmms/ePhaZ.hmm",
    "iPhaZ": "data/hmms/iPhaZ.hmm",
    "OH": "data/hmms/OH.hmm",
    "ArchPhaZ_patatin": "data/hmms/v2/ArchPhaZ_patatin.hmm",
    "ArchPhaZ_hydrolase": "data/hmms/v2/ArchPhaZ_hydrolase.hmm",
}


def hmmsearch(hmm, faa, tbl_out, evalue, cpu=8):
    subprocess.run([HMMSEARCH, "--tblout", tbl_out, "-E", evalue,
                    "--cpu", str(cpu), hmm, faa],
                   check=True, capture_output=True, text=True)


def write_ids(tbl_out, ids_out):
    ids = []
    for line in open(tbl_out):
        if line.startswith("#"):
            continue
        c = line.split()
        if c:
            ids.append(c[0])
    ids = sorted(set(ids))
    with open(ids_out, "w") as f:
        f.write("\n".join(ids) + ("\n" if ids else ""))
    return len(ids)


def extract_faa(ids_path, validated_path, out_path):
    ids = set()
    for line in open(ids_path):
        s = line.strip()
        if s:
            ids.add(s)
    n = 0
    hdr = None
    buf = []
    with open(validated_path) as fin, open(out_path, "w") as fo:
        for line in fin:
            if line.startswith(">"):
                if hdr is not None and hdr in ids:
                    fo.write(">" + hdr + "\n" + "".join(buf) + "\n")
                    n += 1
                hdr = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if hdr is not None and hdr in ids:
            fo.write(">" + hdr + "\n" + "".join(buf) + "\n")
            n += 1
    return n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", nargs=2, metavar=("IDS", "INPUT"))
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--validate-build")
    parser.add_argument("--families")
    parser.add_argument("--cpu", type=int, default=8)
    args = parser.parse_args()
    if args.extract:
        ids_path, input_path = args.extract
        if not args.output:
            parser.error("--extract requires --output")
        n = extract_faa(ids_path, input_path, args.output)
        print(f"extracted {n} sequences -> {args.output}")
        return
    if args.validate_build:
        families = (args.families or "").split()
        for fam in families:
            for tier in ("tier1", "tier2"):
                ids = os.path.join(args.validate_build, f"{fam}_{tier}.ids")
                faa = os.path.join(args.validate_build, f"{fam}_{tier}.faa")
                if not os.path.isfile(ids) or not os.path.isfile(faa):
                    raise RuntimeError(f"missing tier output: {fam} {tier}")
                expected = sum(1 for line in open(ids, encoding="utf-8") if line.strip())
                observed = sum(1 for line in open(faa, encoding="utf-8") if line.startswith(">"))
                if expected != observed:
                    raise RuntimeError(f"tier count mismatch: {fam} {tier}: {expected}/{observed}")
        return
    if args.cpu < 1:
        parser.error("--cpu must be positive")
    os.makedirs(TIERDIR, exist_ok=True)
    for fam, hmm in CURATED.items():
        faa = os.path.join(SEQDIR, f"{fam}_validated.faa")
        if not os.path.exists(faa):
            print(f"[skip] {fam}: 无 validated.faa")
            continue

        # tier2
        hmmsearch(hmm, faa, os.path.join(TIERDIR, f"{fam}_tier2.tbl"), "1e-10", args.cpu)
        n2_ids = write_ids(os.path.join(TIERDIR, f"{fam}_tier2.tbl"),
                           os.path.join(TIERDIR, f"{fam}_tier2.ids"))
        n2 = extract_faa(os.path.join(TIERDIR, f"{fam}_tier2.ids"), faa,
                         os.path.join(TIERDIR, f"{fam}_tier2.faa"))

        # tier1（对 tier2 子集再筛）
        hmmsearch(hmm, os.path.join(TIERDIR, f"{fam}_tier2.faa"),
                  os.path.join(TIERDIR, f"{fam}_tier1.tbl"), "1e-20", args.cpu)
        n1_ids = write_ids(os.path.join(TIERDIR, f"{fam}_tier1.tbl"),
                           os.path.join(TIERDIR, f"{fam}_tier1.ids"))
        n1 = extract_faa(os.path.join(TIERDIR, f"{fam}_tier1.ids"), faa,
                         os.path.join(TIERDIR, f"{fam}_tier1.faa"))

        # 验证 ids 与序列数一致
        ok = (n1_ids == n1 and n2_ids == n2)
        print(f"{fam}: tier1={n1} tier2={n2} "
              f"({'OK' if ok else 'MISMATCH!'})")
        if not ok:
            print(f"[ERROR] {fam}: ids/序列数不一致（tier1 {n1_ids}/{n1}, tier2 {n2_ids}/{n2}）",
                  file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
