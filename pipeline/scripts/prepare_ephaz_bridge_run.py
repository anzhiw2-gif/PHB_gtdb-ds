#!/usr/bin/env python3
"""Prepare a new evidence-gated ePhaZ bridge curation run."""
from __future__ import annotations

import csv
import importlib.util
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("run_context", HERE / "run_context.py")
run_context = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(run_context)

FIELDS = ["candidate_id", "accession", "protein_name", "organism", "sequence_source", "sequence_url", "sequence_length", "primary_reference", "pmid", "doi", "experimental_evidence_type", "substrate", "extracellular_evidence", "architecture_evidence", "completeness_status", "core_overlap_status", "iPhaZ_challenge_status", "decision", "decision_reason", "review_date", "notes"]

BRIDGE = [
 {"candidate_id":"ephaz_bridge_001","accession":"Q51871","protein_name":"Poly(3-hydroxybutyrate) depolymerase A","organism":"Paucimonas lemoignei","sequence_source":"UniProtKB","sequence_url":"https://rest.uniprot.org/uniprotkb/Q51871.fasta","sequence_length":"433","primary_reference":"Jendrossek et al., J Bacteriol. 1995;177:596-607","pmid":"7836292","doi":"10.1128/jb.177.3.596-607.1995","experimental_evidence_type":"cloned_gene_and_biochemical_characterization_of_extracellular_PHA_depolymerase","substrate":"PHB","extracellular_evidence":"primary paper reports extracellular localization; UniProt signal peptide 1-26","architecture_evidence":"UniProt IPR010126 Esterase_PHB; alpha_beta_hydrolase; NCBIfam e_dPHAscl_type1","completeness_status":"complete","core_overlap_status":"not_in_current_core","iPhaZ_challenge_status":"not_challenge","decision":"accept_bridge_candidate","decision_reason":"primary experiment, stable accession, full length, secretory architecture, no iPhaZ-like assignment","review_date":"2026-08-30","notes":"Unreviewed UniProt record directly links nucleotide sequence to PMID 7836292.","sequence":"MRNTLKAAFKLGVISAALLAPFATQAATAGPGAWSSQQTWAADSVNGGNLTGFYYWPATQPVHANGKRALVLVLHGCAQTASGDVINNGDNGYNWKAAADQYGAVILAPNATGNVSSQHCWDYSRTSHSRSTGHEYVLLDLINRFKNDPQYEIDPNQVYVTGLSSGGGETIVLGCIAPDVFAGWASNAGPTPGTTTLQIGAVPSGYTATNAKNNCLSLAGSNSSYFSTQIAGVVWGTSDFTVAPGYNPLMMDAMRQIYGGTFTKQASTSVATGGTNTTYKDSSGRVRTHELSVSGMSHAWPAGTGGQNTNYVTSQYVNYPLFVMDYFFTNNSRAGSGGGTTTTTAGGTTTTTAAGTTTTAATTTTTASSTTTTVAATCYTSSNYAHVTAGRAHNSSGYALANGSNQNMGLNNTFYTSTLKQTSPGYYVIGTCP"},
 {"candidate_id":"ephaz_bridge_002","accession":"Q5SLU4","protein_name":"Carboxylic ester hydrolase (TTHA0199)","organism":"Thermus thermophilus HB8","sequence_source":"UniProtKB","sequence_url":"https://rest.uniprot.org/uniprotkb/Q5SLU4.fasta","sequence_length":"515","primary_reference":"Papaneophytou et al., Appl Microbiol Biotechnol. 2009;83:659-668","pmid":"19214501","doi":"10.1007/s00253-008-1842-2","experimental_evidence_type":"purified_enzyme_and_biochemical_characterization; gene_identified_as_TTHA0199","substrate":"PHB","extracellular_evidence":"primary paper purified extracellular enzyme; UniProt signal peptide 1-19","architecture_evidence":"Carboxylesterase type B domain 22-483; alpha_beta_hydrolase; primary paper reports lipase box","completeness_status":"complete","core_overlap_status":"not_in_current_core","iPhaZ_challenge_status":"not_challenge","decision":"accept_bridge_candidate","decision_reason":"paper names TTHA0199 and stable accession maps to HB8 locus; full length and secretory architecture","review_date":"2026-08-30","notes":"Paper catalytic positions use mature-protein numbering after 19-aa signal peptide.","sequence":"MLRRLLPFLALLGGALAQAFWVETPLGRAQGRLEGGAIAFYGLPYAEAERFRAPKPLKAWPPGVGQEAVACPQAPGITAWFGGPIPLEREDCLVLNVYLPAQIPPPGGFPVMVYLHGGGFTSGAGAEPIYRGHRLSEEGVVVVAPNYRLGPLGFLALPALAEEDPKAVGNYGLLDVLEALRFVRDYIRYFGGDPKNVTLFGESAGGMLVCTLLATPEARGLFQKAIVQSGGCGYVRALEEDYAQGEAWAKARGCDPKDLACLRALPLERLLPEEPTLEATGRFLSNPSLFRTGPFKPHLSPFLLPQDPREALREGKAAGTPLIAGANAEEVAFPSLQALLGPGDWEEAERRLLESGLSREKAQALLAHYRKGVPDPKRAWGEVQTDLTLLCPSLKAARLQAPHAPTYAYLFTFRAPGFEGLGAFHGLELAPLFGNLLERPFLPLFLRQEAQEEAEYLGKKMRRYWTSFAKDGEPKGWPRWPLYREGLLLRLDVPLGLLPDLYEERCGALEVLGLL"},
]

def read_fasta(path: Path):
    records, header, chunks = {}, None, []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if header: records[header[1:].split("|", 1)[0]] = (header, "".join(chunks))
            header, chunks = raw, []
        elif raw.strip(): chunks.append(raw.strip())
    if header: records[header[1:].split("|", 1)[0]] = (header, "".join(chunks))
    return records

def write_fasta(path: Path, records):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records: handle.write(f"{header}\n{sequence}\n")

def prepare(root, run_id, core_source, positive_controls, negative_controls):
    run = run_context.create_run_layout(root, run_id)
    inputs, results = run / "inputs", run / "results"
    protocol = """# ePhaZ bridge curation protocol\n\nDate: 2026-08-30\n\nInclusion requires direct primary experimental evidence for a specific extracellular PHA/PHB depolymerase, an unambiguous stable public protein accession, a complete sequence of at least 200 aa, extracellular/secretory architecture evidence, and no iPhaZ/PhaDED-like confounding.\n\nDiscovery used PubMed, UniProtKB, and NCBI Protein with exact enzyme/locus/organism queries. Annotation-only, unresolved sequence identity, short/incomplete, duplicate, iPhaZ-like, nylon hydrolase, and generic lipase/cutinase records are excluded. This curation does not authorize a GTDB full scan.\n"""
    (inputs / "search_protocol.md").write_text(protocol, encoding="utf-8")
    core = read_fasta(Path(core_source)); controls = read_fasta(Path(positive_controls))
    core_ids = ("B2NHN2", "O05527", "P12625", "Q51718")
    challenge_ids = ("A0ABY7N197", "A0A3Q9BTL6", "A0ABY9RN41", "A0A3N5XWX8", "A0ABT5L4W6", "A0ABT7T1G3")
    write_fasta(inputs / "ephaz_curated_core.faa", [core[x] for x in core_ids])
    write_fasta(inputs / "iPhaZ_like_challenge.faa", [controls[x] for x in challenge_ids])
    shutil.copyfile(negative_controls, inputs / "negative_controls.faa")
    selected = [{key: row.get(key, "") for key in FIELDS} for row in BRIDGE]
    for target in (inputs / "literature_candidates.tsv", results / "ephaz_bridge_candidate.tsv"):
        with target.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(selected)
    write_fasta(results / "ephaz_bridge_candidate.faa", [(f">{x['accession']}|ePhaZ_bridge_candidate|{x['organism']}", x["sequence"]) for x in BRIDGE])
    with (results / "ephaz_bridge_rejections.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(("candidate_id", "rejection_class", "reason")); writer.writerows([("PhaZSex2", "experimental_activity_but_no_retrievable_sequence", "PMID 26156240; no unambiguous stable protein accession found"), ("fkbU", "experimental_activity_but_no_retrievable_sequence", "PMID 23951224; no unambiguous stable protein accession found"), ("PhaZBm", "sequence_paper_identity_unresolved", "PMID 17064368; strain-specific stable protein accession unresolved"), ("six_former_broad_controls", "iPhaZ_like_or_cross_family_confounded", "Six controlled records retained only as challenge set")])
    (results / "bridge_curation_summary.md").write_text("# ePhaZ bridge curation summary\n\n- Screened leads: 6\n- Accepted bridge candidates: 2 (Q51871, Q5SLU4)\n- Rejected or held: 4 lead classes\n- The bridge remains provisional and is not merged into ePhaZ_curated_core.\n- The next step is restricted leave-one-core-out calibration, not a GTDB full scan.\n", encoding="utf-8")
    run_context.write_input_contract(run, run_id=run_id, inputs={"search_protocol": inputs / "search_protocol.md", "literature_candidates": inputs / "literature_candidates.tsv", "curated_core": inputs / "ephaz_curated_core.faa", "bridge_candidates": results / "ephaz_bridge_candidate.faa", "iPhaZ_like_challenge": inputs / "iPhaZ_like_challenge.faa", "negative_controls": inputs / "negative_controls.faa"})
    return run
