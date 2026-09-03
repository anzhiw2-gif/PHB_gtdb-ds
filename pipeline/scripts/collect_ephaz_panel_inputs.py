#!/usr/bin/env python3
"""Capture accession and PMID responses through the local browser CDP proxy."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen


ACCESSION_URLS = {
    "AAB40611.1": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=AAB40611.1&rettype=fasta&retmode=text",
    "O24719": "https://rest.uniprot.org/uniprotkb/O24719.fasta",
    "A0A8W8": "https://rest.uniprot.org/uniprotkb/A0A8W8.fasta",
    "Q9LBN6": "https://rest.uniprot.org/uniprotkb/Q9LBN6.fasta",
    "Q5YEW3": "https://rest.uniprot.org/uniprotkb/Q5YEW3.fasta",
    "Q6UFW4": "https://rest.uniprot.org/uniprotkb/Q6UFW4.fasta",
    "Q84C08": "https://rest.uniprot.org/uniprotkb/Q84C08.fasta",
    "Q71KW6": "https://rest.uniprot.org/uniprotkb/Q71KW6.fasta",
    "P26495": "https://rest.uniprot.org/uniprotkb/P26495.fasta",
    "Q9AGB6": "https://rest.uniprot.org/uniprotkb/Q9AGB6.fasta",
    "J7K890": "https://rest.uniprot.org/uniprotkb/J7K890.fasta",
    "O87189": "https://rest.uniprot.org/uniprotkb/O87189.fasta",
    "Q7WT48": "https://rest.uniprot.org/uniprotkb/Q7WT48.fasta",
    "Q7WT49": "https://rest.uniprot.org/uniprotkb/Q7WT49.fasta",
    "Q5Q138": "https://rest.uniprot.org/uniprotkb/Q5Q138.fasta",
    "Q7X5S3": "https://rest.uniprot.org/uniprotkb/Q7X5S3.fasta",
    "A0A375HYL0": "https://rest.uniprot.org/uniprotkb/A0A375HYL0.fasta",
    "A0A1C9W3H4": "https://rest.uniprot.org/uniprotkb/A0A1C9W3H4.fasta",
    "Q939Q9": "https://rest.uniprot.org/uniprotkb/Q939Q9.fasta",
}

PMID_URLS = {
    "9297825": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=9297825&retmode=xml",
    "17064368": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=17064368&retmode=xml",
    "9177489": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=9177489&retmode=xml",
    "16232882": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=16232882&retmode=xml",
    "15340791": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=15340791&retmode=xml",
    "15995648": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=15995648&retmode=xml",
    "12898135": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12898135&retmode=xml",
    "1989978": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=1989978&retmode=xml",
    "11114905": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=11114905&retmode=xml",
    "15135527": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=15135527&retmode=xml",
    "12813072": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=12813072&retmode=xml",
    "11457823": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=11457823&retmode=xml",
    "18706425": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=18706425&retmode=xml",
    "20516591": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=20516591&retmode=xml",
    "24007310": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=24007310&retmode=xml",
    "28370478": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=28370478&retmode=xml",
}


def cdp_json(path: str, method: str = "GET", body: str | None = None) -> dict:
    request = Request(f"http://localhost:3456/{path}", method=method, data=body.encode() if body else None)
    with urlopen(request, timeout=60) as response:
        payload = response.read().decode("utf-8")
    try:
        return json.loads(payload)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"CDP returned non-JSON for {path}: {payload[:200]!r}") from exc


def capture(url: str, out_path: Path, wait_seconds: float = 2.0) -> None:
    target = cdp_json("new", method="POST", body=url)["targetId"]
    try:
        value = None
        last_error = None
        for attempt in range(8):
            time.sleep(wait_seconds if attempt == 0 else 1.0)
            try:
                candidate = cdp_json(f"eval?target={target}", method="POST", body="document.body.innerText").get("value")
                if isinstance(candidate, str) and candidate.strip():
                    value = candidate
                    break
            except Exception as exc:  # browser page may still be attaching
                last_error = exc
        if value is None:
            raise RuntimeError(f"empty CDP page for {url}; last_error={last_error}")
        out_path.write_text(json.dumps({"url": url, "target": target, "value": value}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    finally:
        try:
            cdp_json(f"close?target={target}")
        except Exception:
            pass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--accessions", nargs="*", default=sorted(ACCESSION_URLS))
    parser.add_argument("--pmids", nargs="*", default=sorted(PMID_URLS))
    args = parser.parse_args()
    args.input_dir.mkdir(parents=True, exist_ok=True)
    for accession in args.accessions:
        capture(ACCESSION_URLS[accession], args.input_dir / f"{accession.replace('.', '_')}_fasta_eval.json")
    for pmid in args.pmids:
        capture(PMID_URLS[pmid], args.input_dir / f"PMID_{pmid}_xml_eval.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
