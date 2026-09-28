#!/usr/bin/env python3
"""08c_tier_rescore.py — 用 curated 金标准 HMM 对命中做三级重评分（tier1/tier2/tier3）

Python 版，替代 bash 版 08c_tier_rescore.sh：cwd 相对路径、无引号陷阱、fail-closed。
  tier1: curated HMM E<1e-20（严格）
  tier2: curated HMM E<1e-10（中等）
  tier3: 现有 validated（宽模型 + 通用验证）

流程（每个核心家族）：
  validated.faa → hmmsearch(tier2, E<1e-10) → tier2.ids + tier2.faa
  tier2.faa    → hmmsearch(tier1, E<1e-20) → tier1.ids + tier1.faa（从 validated.faa 提取）

数据库尺度（Task F4，与 Task 9 的 06_screen.sh 同规）
----------------------------------------------
HMMER 的 E 值满足 `E = Z * P`，Z 是搜索被告知的目标序列总数。若调用不带 `-Z`，
HMMER 默认 `Z = 本次输入 FASTA 的序列数`。本步骤的输入是**已提取的子集**
（`*_validated.faa`、`*_tier2.faa`），一旦使用该默认值，每个子集都会得到**自己的尺度**，
既不可互相比较，也不可与全库 scale 混用。

因此本步骤：
  * 必须显式给出 `--database-size-z`（或环境变量 `PHB_DATABASE_SIZE_Z`），
    缺失/为零/为负/非整数一律 fail-closed，退出码非 0，**没有默认值、没有推导**；
  * 每一次 hmmsearch 都经 `hmmer_command.build_hmmsearch_command(...)` 构造，
    传入**同一个**全库 Z，`--domZ` 永不出现（下游没有任何环节重标 domain E 值）；
  * `--cpu` 默认 1：并发交给外层驱动，单个 HMMER 任务不再独占多线程；
  * 尺度事实由入口写进 `tier_processing_manifest.json`
    （`database_size_Z`、`database_size_basis`、逐片序列数、命令模板）。

用法（服务器 T141，cwd=工作区根）:
  ~/miniconda3/envs/phb_gtdb/bin/python scripts/08c_tier_rescore.py \
      --database-size-z 615969589 --database-size-basis "parent scan-13 shard counts"

  ⚠️ 上面这个 Z 是 **scan-13 全库的实测值**（2026-09-28 对冻结的 100 个分片逐片计数之
  和，见 `runs/20260928_phaded_scan13_z_scale_reconciliation_01/results/shard_record_counts.tsv`）。
  它**取代**了文档中长期沿用的 `~2.92×10⁸` —— 后者是「单分片 2.92M × 100」的外推，
  实测证明**偏低约 2.1 倍**。**不要照抄旧值**：Z 偏低会使 E 值过于宽松（同比例放大），
  正是本步骤的尺度契约要防止的失败模式。换用任何其他父扫描时，Z 必须**实测**其分片计数之和，
  不得外推、不得沿用本文档的历史数字。
"""
import argparse
import importlib.util
import os
import shutil
import subprocess
import sys

HMMSEARCH = shutil.which("hmmsearch") or os.path.expanduser(
    "~/miniconda3/envs/phb_gtdb/bin/hmmsearch")

#: The one shared command builder (``pipeline/scripts/hmmer_command.py``).
#: It sits next to this file, including when the tier entrypoint copies both into
#: ``runs/<run_id>/inputs/``; it is therefore resolved relative to ``__file__``
#: and never through the current working directory.
HMMER_COMMAND_MODULE = "hmmer_command.py"

#: Documented override for the required CLI flag.  It exists so a driver can pass
#: the scale through the environment; it is never a fallback value.
DATABASE_SIZE_ENV_VAR = "PHB_DATABASE_SIZE_Z"
DATABASE_SIZE_BASIS_ENV_VAR = "PHB_DATABASE_SIZE_BASIS"

SEQDIR = "data/screen/family_seqs"
TIERDIR = "data/screen/tiers"
CURATED = {
    "ePhaZ": "data/hmms/ePhaZ.hmm",
    "iPhaZ": "data/hmms/iPhaZ.hmm",
    "OH": "data/hmms/OH.hmm",
    "ArchPhaZ_patatin": "data/hmms/v2/ArchPhaZ_patatin.hmm",
    "ArchPhaZ_hydrolase": "data/hmms/v2/ArchPhaZ_hydrolase.hmm",
}


def validate_database_size_z(value):
    """Return the full-library ``Z`` as a positive ``int``.

    ``bool`` is rejected even though it is an ``int`` in Python: it is never a
    sequence count.  Anything else — empty, zero, negative, decimal, text — is a
    ``ValueError``, because a rescore run that cannot state the full-library
    target count would silently produce per-subset E-values (Task F4).
    """
    if isinstance(value, bool):
        raise ValueError(
            f"--database-size-z must be a positive integer, got {value!r}"
        )
    if isinstance(value, int):
        if value < 1:
            raise ValueError(
                f"--database-size-z must be a positive integer, got {value!r}"
            )
        return value
    if isinstance(value, str):
        text = value.strip()
        if text.isdigit() and int(text) > 0:
            return int(text)
    raise ValueError(
        f"--database-size-z must be a positive integer (no sign, spaces, or "
        f"decimals), got {value!r}"
    )


def load_hmmer_command():
    """Load the sibling ``hmmer_command`` module by path."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), HMMER_COMMAND_MODULE)
    if not os.path.isfile(path):
        raise RuntimeError(
            f"missing {HMMER_COMMAND_MODULE} next to {os.path.abspath(__file__)}: "
            "the shared full-database-Z command builder is required, and a "
            "hand-written hmmsearch call would drop the one -Z scale"
        )
    spec = importlib.util.spec_from_file_location("hmmer_command", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if getattr(module, "DATABASE_SIZE_FLAG", None) != "-Z":
        raise RuntimeError(  # pragma: no cover - guards a broken sibling module
            "hmmer_command does not declare the -Z database-size flag"
        )
    return module


def hmmsearch(hmm, faa, tbl_out, evalue, cpu=1, *, total_targets):
    """Run one curated-HMM search with the shared full-library ``-Z``.

    Every call is built by ``hmmer_command.build_hmmsearch_command`` so the whole
    tier rerun shares exactly one Z (and never ``--domZ``).  ``total_targets`` is
    keyword-only and has no default: there is no path through this function that
    lets HMMER fall back to the subset's own sequence count.

    Returns the argv that was executed, so a caller (or a test) can prove the
    scale that reached HMMER.
    """
    command = load_hmmer_command().build_hmmsearch_command(
        hmm,
        faa,
        total_targets=total_targets,
        evalue=evalue,
        tblout=tbl_out,
        cpu=cpu,
    )
    argv = [HMMSEARCH] + list(command[1:])
    subprocess.run(argv, check=True, capture_output=True, text=True)
    return argv


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


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--extract", nargs=2, metavar=("IDS", "INPUT"))
    parser.add_argument("--input")
    parser.add_argument("--output")
    parser.add_argument("--validate-build")
    parser.add_argument("--families")
    parser.add_argument(
        "--database-size-z",
        default=os.environ.get(DATABASE_SIZE_ENV_VAR) or None,
        help=(
            "REQUIRED for scoring: the COMPLETE library target sequence count Z "
            "passed to every hmmsearch -Z (never this subset's own count)"
        ),
    )
    parser.add_argument(
        "--database-size-basis",
        default=os.environ.get(DATABASE_SIZE_BASIS_ENV_VAR) or None,
        help="where the full-library Z came from (recorded by the entrypoint)",
    )
    parser.add_argument("--cpu", type=int, default=1)
    args = parser.parse_args(argv)
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
    # Full-database scale is fail-closed (Task F4).  This step's inputs are
    # extracted subsets, so without an explicit full-library Z every E-value it
    # reports would be on that subset's own scale.
    if args.database_size_z is None or (
        isinstance(args.database_size_z, str) and not args.database_size_z.strip()
    ):
        print(
            "[ERROR] --database-size-z is required: it must be the COMPLETE library "
            "target sequence count Z (all shards, not this rescore subset).",
            file=sys.stderr,
        )
        print(
            "        Without it HMMER falls back to Z = this subset's sequence "
            "count, producing incomparable per-subset E-values.",
            file=sys.stderr,
        )
        print(
            "        Pass --database-size-z <positive integer> "
            "--database-size-basis <source>.  This step is fail-closed: there is "
            "no default.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    try:
        database_size_z = validate_database_size_z(args.database_size_z)
    except ValueError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
    if args.cpu < 1:
        parser.error("--cpu must be positive")
    builder = load_hmmer_command()
    basis = args.database_size_basis or "unspecified_unverified"
    print(
        f"full-database scale: Z={database_size_z} ({basis}); "
        f"command template: {builder.command_template()}"
    )
    os.makedirs(TIERDIR, exist_ok=True)
    for fam, hmm in CURATED.items():
        faa = os.path.join(SEQDIR, f"{fam}_validated.faa")
        if not os.path.exists(faa):
            print(f"[skip] {fam}: 无 validated.faa")
            continue

        # tier2
        hmmsearch(hmm, faa, os.path.join(TIERDIR, f"{fam}_tier2.tbl"), "1e-10",
                  args.cpu, total_targets=database_size_z)
        n2_ids = write_ids(os.path.join(TIERDIR, f"{fam}_tier2.tbl"),
                           os.path.join(TIERDIR, f"{fam}_tier2.ids"))
        n2 = extract_faa(os.path.join(TIERDIR, f"{fam}_tier2.ids"), faa,
                         os.path.join(TIERDIR, f"{fam}_tier2.faa"))

        # tier1（对 tier2 子集再筛）
        hmmsearch(hmm, os.path.join(TIERDIR, f"{fam}_tier2.faa"),
                  os.path.join(TIERDIR, f"{fam}_tier1.tbl"), "1e-20",
                  args.cpu, total_targets=database_size_z)
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
    try:
        main()
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
