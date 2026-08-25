#!/usr/bin/env python3
"""09i_tree_manifest.py — 系统发育树登记清单（审计用）

扫描 results/trees_tier1/ 下已提交的树文件，对每棵树登记：
  - 工具与树类型（fasttree / iqtree 抽样或全量 / cdhit 去冗余）
  - 实际叶数（解析 Newick 叶标签；标签以字母开头，避免把内部支持值误计为叶）
  - 树文件 SHA-256 与输入 tier1 FASTA 的 SHA-256
  - 抽样名单（results/trees_tier1/{family}.{kind}.leaves.list）

输出: results/trees_tier1/tree_manifest.tsv（UTF-8，制表符分隔）

树类型判定（按文件名）：
  *.fasttree.nwk      → FastTree 抽样树（seed=42, N=1000）
  *.cdhit.treefile    → IQ-TREE2，输入为服务器侧 CD-HIT 去冗余 fasta
  *.treefile          → IQ-TREE2；n_tier ≤ 2000 为全量，否则抽样（seed=42, N=2000）
用法: python 09i_tree_manifest.py [--treedir results/trees_tier1] [--tierdir data/screen/tiers]
"""
import argparse
import glob
import hashlib
import os
import re
import sys

# 叶标签：以字母开头，后接字母/数字/下划线/点/竖线，紧跟 ':'（分支长度）。
# 内部节点支持值（如 FastTree 的 )0.986:）以数字开头，不会误匹配。
LEAF_RE = re.compile(r"([A-Za-z][A-Za-z0-9_.|]*):")


def sha256(path):
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_leaves(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    return sorted(set(LEAF_RE.findall(text)))


def classify(base):
    """返回 (family, tool, kind)。"""
    if base.endswith(".fasttree.nwk"):
        return base[: -len(".fasttree.nwk")], "FastTree", "fasttree"
    if base.endswith(".cdhit.treefile"):
        return base[: -len(".cdhit.treefile")], "IQ-TREE2", "cdhit"
    if base.endswith(".treefile"):
        return base[: -len(".treefile")], "IQ-TREE2", "iqtree"
    if base.endswith(".tree"):
        return base[: -len(".tree")], "IQ-TREE2", "iqtree"
    return base.rsplit(".", 1)[0], "unknown", "unknown"


def classify_status(n_leaves, current_tier_count, recorded_input_sha256, current_input_sha256):
    """Classify tree evidence against the currently registered tier input."""
    if not n_leaves:
        return "EMPTY_OR_PARSE_FAIL"
    if current_input_sha256 is None:
        return "input_not_registered"
    if recorded_input_sha256 and recorded_input_sha256 != current_input_sha256:
        return "stale_input"
    if current_tier_count and n_leaves > current_tier_count:
        return "leaf_count_exceeds_input"
    if not recorded_input_sha256:
        return "unverified_input"
    return "ok"


def input_hashes(recorded_input_sha256, current_input_sha256):
    """Keep historical tree provenance separate from the current tier input."""
    return {
        "recorded_input_sha256": recorded_input_sha256 or "",
        "current_input_sha256": current_input_sha256 or "",
    }


def read_previous_manifest(path):
    """Read prior tree/input provenance before writing a replacement manifest."""
    previous = {}
    if not os.path.exists(path):
        return previous
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        required = {"tree_file", "tree_sha256"}
        if not required.issubset(index):
            return previous
        input_column = "recorded_input_sha256" if "recorded_input_sha256" in index else "input_sha256"
        if input_column not in index:
            return previous
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(index.values()):
                continue
            key = (fields[index["tree_file"]], fields[index["tree_sha256"]])
            previous[key] = fields[index[input_column]]
    return previous


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--treedir", default="results/trees_tier1")
    ap.add_argument("--tierdir", default="data/screen/tiers")
    args = ap.parse_args()
    out = os.path.join(args.treedir, "tree_manifest.tsv")
    previous = read_previous_manifest(out)

    tree_files = sorted(set(
        glob.glob(os.path.join(args.treedir, "*.treefile")) +
        glob.glob(os.path.join(args.treedir, "*.fasttree.nwk")) +
        glob.glob(os.path.join(args.treedir, "*.tree"))
    ))
    if not tree_files:
        print(f"[ERROR] 未在 {args.treedir} 找到树文件", file=sys.stderr)
        sys.exit(1)

    rows = []
    for tf in tree_files:
        base = os.path.basename(tf)
        family, tool, kind = classify(base)

        try:
            leaves = parse_leaves(tf)
        except (OSError, UnicodeError):
            leaves = []
        n_leaves = len(leaves)

        leaf_list = os.path.join(args.treedir, f"{family}.{kind}.leaves.list")
        with open(leaf_list, "w", encoding="utf-8") as f:
            f.write("\n".join(leaves) + ("\n" if leaves else ""))

        # 输入 tier1.faa（cdhit 树的实际输入是服务器侧去冗余 fasta，本地无）
        if kind == "cdhit":
            input_faa = ""
            input_sha256 = None
            n_tier = 0
            note = f"CD-HIT 去冗余后建树（{n_leaves} 叶；输入为服务器侧去冗余 fasta，待登记）"
        else:
            input_faa = os.path.join(args.tierdir, f"{family}_tier1.faa")
            input_sha256 = sha256(input_faa)
            n_tier = 0
            if os.path.exists(input_faa):
                with open(input_faa, encoding="utf-8") as f:
                    n_tier = sum(1 for line in f if line.startswith(">"))
            if kind == "fasttree":
                note = f"FastTree 抽样树（seed=42, N=1000；叶数可因去重/过滤低于 1000）"
            else:  # iqtree
                if n_tier > 2000:
                    note = f"抽样树（seed=42, N=2000；叶数可因去重/过滤低于 2000）"
                else:
                    note = f"全量树（{n_tier} 序列 ≤ N=2000，未抽样）"

        tree_sha = sha256(tf)
        recorded_sha = previous.get((tf, tree_sha), "")
        hashes = input_hashes(recorded_sha, input_sha256)
        status = classify_status(n_leaves, n_tier, recorded_sha or None, input_sha256)
        if status == "stale_input":
            note += "；树文件未随当前 tier1 输入重建，需重新建树"
        rows.append({
            "family": family,
            "kind": kind,
            "tree_file": tf,
            "tool": tool,
            "n_leaves": n_leaves,
            "tree_sha256": tree_sha,
            "input_faa": input_faa,
            "recorded_input_sha256": hashes["recorded_input_sha256"],
            "current_input_sha256": hashes["current_input_sha256"],
            "tier1_seq_count": n_tier,
            "leaf_list": leaf_list,
            "status": status,
            "note": note,
        })

    cols = ["family", "kind", "tree_file", "tool", "n_leaves", "tree_sha256",
            "input_faa", "recorded_input_sha256", "current_input_sha256",
            "tier1_seq_count", "leaf_list", "status", "note"]
    with open(out, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"tree_manifest.tsv -> {out}")
    for r in rows:
        print(f"  {r['family']} [{r['kind']}]: {r['tool']} 叶数={r['n_leaves']} 状态={r['status']}")


if __name__ == "__main__":
    main()
