#!/usr/bin/env python3
import argparse, csv
from pathlib import Path

def main():
    p=argparse.ArgumentParser(); p.add_argument('--input',type=Path,required=True); p.add_argument('--output',type=Path,required=True); a=p.parse_args()
    with a.input.open(newline='') as src, a.output.open('w',newline='') as dst:
        r=csv.DictReader(src,delimiter='\t'); w=csv.writer(dst,delimiter='\t'); w.writerow(['genome','family','copies','best_E','best_score'])
        for row in r: w.writerow([row['genome'],row['family'],row['copies'],'',''])
if __name__=='__main__': main()
