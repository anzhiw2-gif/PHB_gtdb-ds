#!/usr/bin/env python3
"""filter_long_seqs.py — 流式过滤 FASTA 中超长序列（>MAXLEN aa，HMMER 硬限制 100K）
用法: python filter_long_seqs.py <in.faa> <out.faa> [maxlen]
"""
import sys

MAXLEN = 100000


def filter_fasta(inpath, outpath, maxlen=MAXLEN):
    kept = dropped = 0
    with open(inpath) as fin, open(outpath, "w") as fout:
        name = None
        seq = []
        seqlen = 0

        def flush():
            nonlocal kept, dropped, name, seq, seqlen
            if name is not None:
                if seqlen <= maxlen:
                    fout.write(name + "\n")
                    fout.write("".join(seq))
                    kept += 1
                else:
                    dropped += 1
            name = None
            seq = []
            seqlen = 0

        for line in fin:
            if line.startswith(">"):
                flush()
                name = line.rstrip("\n")
            else:
                seq.append(line)
                seqlen += len(line.rstrip("\n"))
        flush()
    return kept, dropped


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    inpath, outpath = sys.argv[1], sys.argv[2]
    maxlen = int(sys.argv[3]) if len(sys.argv) > 3 else MAXLEN
    kept, dropped = filter_fasta(inpath, outpath, maxlen)
    print(f"{inpath}: kept={kept} dropped={dropped} (> {maxlen} aa)")
