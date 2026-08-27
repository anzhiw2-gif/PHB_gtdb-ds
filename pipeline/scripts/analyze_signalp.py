#!/usr/bin/env python3
"""analyze_signalp.py — 分析 ePhaZ tier1 的 SignalP 结果与 N 端截断线索

回答：ePhaZ 38,692 条里仅 ~42.8% 是分泌型(SP)，其余"无信号肽"里有多少是
N 端截断/预测不完整导致的假阴性？

截断代理信号（无需重跑 Pyrodigal，从序列本身判断）：
  1) 不以 Met 开头 —— Pyrodigal 预测的完整 CDS 起始密码子(ATG/GTG/TTG)都译成 Met，
     若 N 端被 contig 边界或预测起始点截断，则第一个残基非 Met。
  2) 序列内部含 '*'（终止符）—— 说明读码框中断/截断。
  3) 长度明显偏短 —— 完整胞外 ePhaZ 通常 >400 aa。

输入: data/screen/tiers/ePhaZ_tier1.faa
      results/signalp/ePhaZ/prediction_results.txt
输出: results/tables/signalp_ePhaZ_analysis.tsv + 打印摘要
"""
import os
from collections import Counter

TIER = "data/screen/tiers/ePhaZ_tier1.faa"
PRED = "results/signalp/ePhaZ/prediction_results.txt"
OUT = "results/tables/signalp_ePhaZ_analysis.tsv"


def read_faa(path):
    seqs = {}
    hdr = None
    buf = []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None:
                seqs[hdr] = "".join(buf)
            hdr = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if hdr is not None:
        seqs[hdr] = "".join(buf)
    return seqs


def read_pred(path):
    pred = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            pred[c[0]] = c[1]
    return pred


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return 0
    if n % 2:
        return xs[n // 2]
    return (xs[n // 2 - 1] + xs[n // 2]) / 2


def main():
    seqs = read_faa(TIER)
    pred = read_pred(PRED)
    print(f"tier1 序列: {len(seqs)}, SignalP 记录: {len(pred)}")

    # 1) Prediction 类型分布
    print("\n=== SignalP Prediction 类型分布 ===")
    pc = Counter(pred.values())
    for k, v in pc.most_common():
        print(f"  {k}: {v} ({v/len(pred)*100:.1f}%)")

    # 有信号肽 = SP + LIPO + TAT + TATLIPO + PILIN；无 = OTHER
    has_sp = {k for k, v in pred.items() if v != "OTHER"}
    no_sp = {k for k, v in pred.items() if v == "OTHER"}
    print(f"\n有信号肽(SP+LIPO+TAT+TATLIPO+PILIN): {len(has_sp)} ({len(has_sp)/len(pred)*100:.1f}%)")
    print(f"无信号肽(OTHER): {len(no_sp)} ({len(no_sp)/len(pred)*100:.1f}%)")

    # 2) 长度分布对比
    def lens(accs):
        return [len(seqs[a]) for a in accs if a in seqs]

    l_sp = lens(has_sp)
    l_no = lens(no_sp)
    print("\n=== 长度分布（aa）===")
    print(f"  有信号肽组: n={len(l_sp)} min={min(l_sp)} median={median(l_sp):.0f} max={max(l_sp)}")
    print(f"  无信号肽组: n={len(l_no)} min={min(l_no)} median={median(l_no):.0f} max={max(l_no)}")

    # 3) N 端截断代理信号（无信号肽组）
    print("\n=== 无信号肽(OTHER)组 截断代理信号 ===")
    other_seqs = {a: seqs[a] for a in no_sp if a in seqs}
    n_other = len(other_seqs)
    no_met = [a for a, s in other_seqs.items() if not s.startswith("M")]
    internal_star = [a for a, s in other_seqs.items() if "*" in s[:-1]]
    short = [a for a, s in other_seqs.items() if len(s) < 300]
    print(f"  总数: {n_other}")
    print(f"  不以 Met 开头(N端截断嫌疑): {len(no_met)} ({len(no_met)/n_other*100:.1f}%)")
    print(f"  序列内部含 '*'(读码框中断): {len(internal_star)} ({len(internal_star)/n_other*100:.1f}%)")
    print(f"  长度 < 300 aa(偏短): {len(short)} ({len(short)/n_other*100:.1f}%)")
    # 组合：N端截断嫌疑（不以Met 或 内部* 或 短）
    trunc = {a for a in no_met} | {a for a in internal_star} | {a for a in short}
    print(f"  任一截断嫌疑(并集): {len(trunc)} ({len(trunc)/n_other*100:.1f}%)")

    # 4) 对照组：有信号肽组 不以 Met 开头（应极低，作为方法校验）
    sp_seqs = {a: seqs[a] for a in has_sp if a in seqs}
    sp_no_met = [a for a, s in sp_seqs.items() if not s.startswith("M")]
    print(f"\n=== 对照：有信号肽组 不以 Met 开头 ===")
    print(f"  {len(sp_no_met)} / {len(sp_seqs)} ({len(sp_no_met)/max(len(sp_seqs),1)*100:.1f}%)")

    # 写 TSV
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"tier1_total\t{len(seqs)}\n")
        for k, v in pc.most_common():
            f.write(f"pred_{k}\t{v}\n")
        f.write(f"has_signal_peptide\t{len(has_sp)}\n")
        f.write(f"no_signal_peptide\t{len(no_sp)}\n")
        f.write(f"sp_len_median\t{median(l_sp):.1f}\n")
        f.write(f"other_len_median\t{median(l_no):.1f}\n")
        f.write(f"other_no_met\t{len(no_met)}\n")
        f.write(f"other_internal_star\t{len(internal_star)}\n")
        f.write(f"other_len_lt300\t{len(short)}\n")
        f.write(f"other_trunc_union\t{len(trunc)}\n")
        f.write(f"sp_no_met\t{len(sp_no_met)}\n")
    print(f"\n分析结果 -> {OUT}")


if __name__ == "__main__":
    main()
