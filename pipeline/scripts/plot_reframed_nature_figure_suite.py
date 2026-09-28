#!/usr/bin/env python3
"""Reframed Nature-style candidate figures: ePhaZ/iPhaZ main line + archaeal side branch."""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

MAIN_FAMILIES = ("ePhaZ", "iPhaZ")
WIDTH_IN = 183 / 25.4
COLORS = {"ePhaZ":"#0F4D92", "iPhaZ":"#42949E", "secreted":"#4C9F70", "review":"#C99000", "gray":"#777777", "light":"#D9D9D9", "dark":"#262626", "arch":"#7B61A8", "red":"#B64342"}
PALETTE = {"dark": COLORS["dark"], "blue": COLORS["ePhaZ"], "teal": COLORS["iPhaZ"]}
SUBTYPE_LABELS = {
 "ePhaZ_SCL_secreted_evidence":"ePhaZ SCL secreted", "ePhaZ_core_nonsecreted_evidence":"ePhaZ core nonsecreted",
 "ePhaZ_tier2_review_evidence":"ePhaZ tier2 review", "ePhaZ_competition_review_evidence":"ePhaZ competition review",
 "iPhaZ_intracellular_evidence":"iPhaZ intracellular", "iPhaZ_tier2_review_evidence":"iPhaZ tier2 review"}
LAYER_ORDER = ["PhaZh1_like_high","PhaZh1_like_review","archaeal_patatin_exploratory","non_PhaZh1_patatin"]
LAYER_LABELS = {"PhaZh1_like_high":"PhaZh1-like high","PhaZh1_like_review":"PhaZh1-like review","archaeal_patatin_exploratory":"patatin exploratory","non_PhaZh1_patatin":"non-PhaZh1 patatin"}
LAYER_COLORS = {"PhaZh1_like_high":"#0F4D92","PhaZh1_like_review":"#C99000","archaeal_patatin_exploratory":"#D9D9D9","non_PhaZh1_patatin":"#B64342"}

def read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0: raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

def subtype_counts(frame: pd.DataFrame) -> pd.DataFrame:
    counts=frame["pha_subtype_evidence"].value_counts()
    rows=[{"source_category":k,"category":SUBTYPE_LABELS[k],"count":int(counts.get(k,0))} for k in SUBTYPE_LABELS if counts.get(k,0)]
    return pd.DataFrame(rows)

def species_from_taxonomy(taxonomy: str) -> str:
    for token in taxonomy.split(";"):
        if token.startswith("s__") and token[3:].strip(): return token[3:].strip()
    return "unclassified"

def validate_export_bundle(output_dir: Path, stems: list[str]) -> list[str]:
    return [f"{stem}{ext}" for stem in stems for ext in (".svg",".pdf",".tiff",".png") if not (output_dir/f"{stem}{ext}").is_file() or (output_dir/f"{stem}{ext}").stat().st_size==0]

def validate_main_family_scope(frame: pd.DataFrame) -> list[str]:
    return sorted(set(frame.get("family", [])) - set(MAIN_FAMILIES))

def signalp_scope_counts(frame: pd.DataFrame) -> dict[str, int]:
    scoped = frame[(frame.family == "ePhaZ") & (frame.source_layer == "ePhaZ_curated_core")]
    return {"input_count": int(len(scoped)), "excluded_iPhaZ": int((frame.family == "iPhaZ").sum()), "excluded_ePhaZ_tier2": int(((frame.family == "ePhaZ") & (frame.source_layer != "ePhaZ_curated_core")).sum())}

def genome_distribution(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame[frame.family.isin(MAIN_FAMILIES)].drop_duplicates(["genome","family"])
    return frame.groupby(["phylum","family"], as_index=False).size().rename(columns={"size":"genomes"})

def competition_display_groups(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep the small competition outcomes visible without rescaling their counts."""
    ordered = (
        frame.set_index("assignment")
        .reindex(["ePhaZ_like", "iPhaZ_like", "ambiguous"])
        .dropna()
        .reset_index()
    )
    return ordered.iloc[:1].copy(), ordered.iloc[1:].copy()

def workflow_table() -> pd.DataFrame:
    return pd.DataFrame([
      {"step":"GTDB Release 11 R232","unit":"representative genomes","count":199923,"side_branch":False,"tool":"GTDB taxonomy + metadata"},
      {"step":"Protein prediction / sequence binding","unit":"genomes","count":199923,"side_branch":False,"tool":"Pyrodigal + FAA/GFF QC"},
      {"step":"HMM raw hits","unit":"hit rows","count":6769772,"side_branch":False,"tool":"HMMER 3.4"},
      {"step":"Profile search and tiering","unit":"hit rows","count":6740900,"side_branch":False,"tool":"HMMER 3.4"},
      {"step":"ePhaZ/iPhaZ deduplicated candidates","unit":"proteins","count":109087,"side_branch":False,"tool":"family-specific evidence layering"},
      {"step":"Genome-level ePhaZ/iPhaZ representation","unit":"genome×family","count":76141,"side_branch":False,"tool":"taxonomy-bound de-duplication"},
    ])

def resolve_input_tables(run_dir: Path) -> dict[str, Path]:
    names = {"overall":"overall_pha_phb.tsv","archaea":"archaeal_phazh1_full_audit.tsv","domain":"domain_annotation.tsv","calibration":"positive_negative_calibration.tsv","manual":"cross_phylum_manual_review.tsv","taxonomy":"tier1_genome_family.tsv","competition":"no_signal_competition.tsv"}
    paths = {k: run_dir/"inputs"/v for k,v in names.items()}
    missing = [str(v) for v in paths.values() if not v.is_file() or v.stat().st_size == 0]
    if missing: raise FileNotFoundError("missing inputs: "+", ".join(missing))
    return paths

def style(ax, grid="x"):
    ax.grid(axis=grid, color="#E8E8E8", lw=.6, zorder=0); ax.tick_params(length=2.5,width=.7,pad=2)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
def label(ax, text): ax.text(.01,.99,text,transform=ax.transAxes,va="top",fontweight="bold",fontsize=8)
def bar_label(value: int, ymax: float):
    text=f"{int(value):,}"
    return (text,"inside","white") if value >= ymax*0.45 else (text,"outside",COLORS["dark"])
def save(fig, out: Path, stem: str):
    out.mkdir(parents=True, exist_ok=True)
    for ext,dpi in (("svg",None),("pdf",None),("tiff",600),("png",300)):
        kw={"facecolor":"white","bbox_inches":"tight"};
        if dpi: kw["dpi"]=dpi
        fig.savefig(out/f"{stem}.{ext}", **kw)
    plt.close(fig)
def write(frame,path): path.parent.mkdir(parents=True,exist_ok=True); frame.to_csv(path,sep="\t",index=False,lineterminator="\n")

def file_record(path: Path) -> dict[str, str | int]:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

def write_input_contract(run: Path, inputs: dict[str, Path]) -> None:
    contract = {
        "run_id": run.name,
        "purpose": "Nature-style visualization only: ePhaZ/iPhaZ main line with a detached archaeal PhaZh1-like supplementary branch.",
        "status": "completed_candidate_only",
        "input_snapshot": {name: file_record(path) for name, path in sorted(inputs.items())},
        "renderer": file_record(Path(__file__).resolve()),
        "environment": {"python_executable": sys.executable, "python_version": sys.version.split()[0], "pandas_version": pd.__version__, "matplotlib_version": matplotlib.__version__},
        "upstream_assets": {
            "GTDB_taxonomy": {"version": "Release 11 R232", "path": "pending", "bytes": "pending", "sha256": "pending"},
            "GTDB_metadata": {"version": "Release 11 R232", "path": "pending", "bytes": "pending", "sha256": "pending"},
            "GTDB_tree": {"version": "Release 11 R232", "path": "pending", "bytes": "pending", "sha256": "pending"},
            "HMM_profiles": {"version": "pending", "path": "pending", "bytes": "pending", "sha256": "pending"},
        },
        "scope": {
            "main_families": list(MAIN_FAMILIES),
            "signalp_scope": "ePhaZ curated-core tier1 only (38,692 proteins)",
            "competition_scope": "ePhaZ curated-core OTHER only (16,836 proteins)",
            "archaeal_side_branch": "PhaZh1-like domain/motif/neighborhood and control evidence; excluded from main family totals.",
        },
        "phenotype_boundary": "All displayed HMM, SignalP, domain, motif, neighborhood and taxonomy evidence denotes candidate homology or functional potential, not validated PHB degradation phenotype.",
        "non_modification_statement": "This visualization run did not execute a formal scan or modify a registry, deploy directory, server state, or historical run.",
    }
    (run / "input_contract.json").write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")

def build_data(run: Path):
    p=resolve_input_tables(run); overall=read_tsv(p["overall"]); tax=read_tsv(p["taxonomy"]); arch=read_tsv(p["archaea"]); dom=read_tsv(p["domain"]); cal=read_tsv(p["calibration"]); manual=read_tsv(p["manual"]); comp=read_tsv(p["competition"])
    if validate_main_family_scope(overall): raise ValueError("overall table contains non-main family")
    main_tax=tax[tax.family.isin(MAIN_FAMILIES)].drop_duplicates(["genome","family"])
    subtype=overall.groupby(["family","pha_subtype_evidence"],as_index=False).size().rename(columns={"size":"proteins"})
    fam=overall.groupby("family",as_index=False).agg(proteins=("accession","nunique"),genomes=("genome","nunique"))
    sig=overall[(overall.family=="ePhaZ")&(overall.source_layer=="ePhaZ_curated_core")].groupby("signal_type",as_index=False).size().rename(columns={"size":"proteins"})
    ph=genome_distribution(main_tax)
    gen=main_tax.assign(genus=main_tax.gtdb_taxonomy.str.extract(r";g__([^;]+)")[0].fillna("unclassified"),species=main_tax.gtdb_taxonomy.str.extract(r";s__([^;]+)")[0].fillna("unclassified"))
    genus=gen.groupby("genus",as_index=False).genome.nunique().rename(columns={"genome":"genomes"}).sort_values("genomes",ascending=False).head(15)
    co=gen.groupby("genome").family.nunique().value_counts().rename_axis("family_count").reset_index(name="genomes").sort_values("family_count")
    layers=arch.groupby("layer",as_index=False).size().rename(columns={"size":"records"}); layers["layer"]=pd.Categorical(layers.layer,LAYER_ORDER,ordered=True); layers=layers.sort_values("layer")
    motif=dom.groupby(["layer","nterm_gxsxg_motif"],as_index=False).size().rename(columns={"size":"records"})
    comp_counts=comp.assignment.value_counts().rename_axis("assignment").reset_index(name="proteins") if "assignment" in comp else pd.DataFrame(columns=["assignment","proteins"])
    return {"overall":overall,"taxonomy":main_tax,"subtype":subtype,"family":fam,"signal":sig,"phylum":ph,"genus":genus,"co":co,"arch":arch,"layers":layers,"motif":motif,"domain":dom,"cal":cal,"manual":manual,"competition":comp_counts,"workflow":workflow_table()}

def figure1(d,out):
    fig=plt.figure(figsize=(WIDTH_IN,4.2)); gs=fig.add_gridspec(1,2,width_ratios=[1.05,1.2],wspace=.28); ax=fig.add_subplot(gs[0]); fu=fig.add_subplot(gs[1]); ax.axis("off")
    steps=[("1","GTDB Release 11 R232","199,923 representative genomes","GTDB taxonomy + metadata"),("2","Protein prediction and QC","FAA/GFF sequence binding","Pyrodigal"),("3","HMM raw hits","6,769,772 hit rows","HMMER 3.4"),("4","HMMER profile search","HMMER 3.4; family-specific thresholds","HMMER 3.4"),("5","ePhaZ / iPhaZ evidence layering","109,087 deduplicated proteins","SignalP6 only on ePhaZ core")]
    y=np.linspace(.86,.18,len(steps))
    for i,(n,title,note,tool) in enumerate(steps):
        box=FancyBboxPatch((.08,y[i]-.065),.84,.13,boxstyle="round,pad=.012,rounding_size=.018",fc="#EAF2F8" if i<4 else "#D8ECEA",ec="#666",lw=.8); ax.add_patch(box); ax.text(.14,y[i]+.03,n,ha="center",va="center",fontsize=8,fontweight="bold"); ax.text(.52,y[i]+.028,title,ha="center",va="center",fontsize=7,fontweight="bold"); ax.text(.52,y[i]-.003,note,ha="center",va="center",fontsize=5.5); ax.text(.52,y[i]-.033,tool,ha="center",va="center",fontsize=5,color="#555")
        if i<len(steps)-1: ax.add_patch(FancyArrowPatch((.5,y[i]-.08),(.5,y[i+1]+.08),arrowstyle="-|>",mutation_scale=9,lw=.8,color="#777"))
    ax.text(.5,.04,"Archaeal PhaZh1-like audit is a detached supplementary branch (2,307 records).",ha="center",fontsize=5.8,color=COLORS["arch"]); label(ax,"a"); ax.set_title("Main GTDB-to-candidate workflow",fontsize=8,fontweight="bold",pad=10)
    w=d["workflow"]; names=["GTDB\nR232","HMM raw\nhits","HMMER\naccepted","e/i\nproteins","genome×\nfamily"]; vals=[199923,6769772,6740900,109087,76141]; y2=np.arange(len(vals)); fu.barh(y2,vals,color=["#C4D3E2","#8FAECC","#6F95BC",COLORS["ePhaZ"],COLORS["iPhaZ"]],zorder=2); fu.set_yticks(y2, names,fontsize=6); fu.invert_yaxis(); fu.set_xlabel("Count (log10 scale)",fontsize=7); fu.set_xscale("log"); fu.set_title("Main-line data funnel",fontsize=8,fontweight="bold"); style(fu); label(fu,"b")
    for yi,v in zip(y2,vals): fu.text(v*1.08,yi,f"{v:,}",va="center",fontsize=6)
    fig.text(.5,.01,"Main funnel uses ePhaZ/iPhaZ only; counts are not phenotype confirmations. Archaeal records are excluded from these totals.",ha="center",fontsize=6,color="#666"); fig.subplots_adjust(top=.84,bottom=.14,left=.06,right=.98); save(fig,out,"figure_1_workflow_and_funnel")

def figure2(d,out):
    fig,ax=plt.subplots(figsize=(WIDTH_IN,3.8)); rows=[("ePhaZ · total",70682,"#0F4D92"),("  curated secreted",21856,"#4C9F70"),("  curated nonsecreted",16818,"#6FA6C8"),("  tier2 review",31990,"#C99000"),("  competition review",18,"#B64342"),("iPhaZ · total",38405,"#42949E"),("  tier1",32926,"#42949E"),("  tier2 review",5479,"#C99000")]; y=np.arange(len(rows)); ax.barh(y,[r[1] for r in rows],color=[r[2] for r in rows],zorder=2); ax.set_yticks(y,[r[0] for r in rows],fontsize=6.5); ax.invert_yaxis(); ax.set_xlabel("Deduplicated candidate proteins"); ax.set_title("Two-family hierarchy and evidence layers (n = 109,087)",fontsize=9,fontweight="bold"); style(ax); label(ax,"a");
    for yi,(_,v,_) in zip(y,rows): ax.text(v*1.01,yi,f"{v:,}",va="center",fontsize=6)
    fig.text(.5,.035,"ePhaZ and iPhaZ are separate families; SignalP-derived layers are not symmetric evidence for iPhaZ.",ha="center",fontsize=6,color="#666"); fig.subplots_adjust(left=.24,right=.96,top=.85,bottom=.25); save(fig,out,"figure_2_family_hierarchy")

def figure3(d,out):
    sig=d["signal"].set_index("signal_type").reindex(["SP","LIPO","TAT","TATLIPO","OTHER"]).fillna(0); fig,ax=plt.subplots(figsize=(WIDTH_IN,3.3)); vals=sig.proteins.astype(int); y=np.arange(len(vals)); colors=[COLORS["secreted"],"#6FA6C8","#8C6BB1","#B07AA1",COLORS["light"]]; ax.barh(y,vals,color=colors,zorder=2); ax.set_yticks(y,sig.index); ax.invert_yaxis(); ax.set_xlabel("ePhaZ curated-core tier1 proteins"); ax.set_title("SignalP6 localization evidence is ePhaZ-only (input n = 38,692)",fontsize=9,fontweight="bold"); style(ax); label(ax,"a");
    for yi,v in zip(y,vals): ax.text(v+250,yi,f"{v:,}",va="center",fontsize=6)
    fig.text(.5,.035,"SignalP6 was not run as a family-wide layer for iPhaZ; localization does not establish substrate specificity.",ha="center",fontsize=6,color="#666"); fig.subplots_adjust(left=.18,right=.96,top=.85,bottom=.25); save(fig,out,"figure_3_ephaz_signalp_scope")

def figure4(d,out):
    ph=d["phylum"]; piv=ph.pivot(index="phylum",columns="family",values="genomes").fillna(0); piv["total"]=piv.sum(axis=1); piv=piv.sort_values("total").tail(18); fig,ax=plt.subplots(figsize=(WIDTH_IN,4.4)); y=np.arange(len(piv)); ax.barh(y,piv.get("ePhaZ",0),color=COLORS["ePhaZ"],label="ePhaZ"); ax.barh(y,piv.get("iPhaZ",0),left=piv.get("ePhaZ",0),color=COLORS["iPhaZ"],label="iPhaZ"); ax.set_yticks(y,piv.index,fontsize=6); ax.set_xlabel("Unique tier1 genomes (family-specific)"); ax.set_title("GTDB phylum distribution of the two main families",fontsize=9,fontweight="bold"); ax.legend(fontsize=6); style(ax); label(ax,"a"); fig.text(.5,.01,"Counts are unique genome×family assignments from the taxonomy-bound tier1 table; families may co-occur in one genome.",ha="center",fontsize=6,color="#666"); fig.subplots_adjust(left=.18,right=.97,top=.85,bottom=.15); save(fig,out,"figure_4_main_family_phylum_distribution")

def figure5(d,out):
    fig,axs=plt.subplots(1,2,figsize=(WIDTH_IN,3.8),gridspec_kw={"width_ratios":[1.15,.85]}); g=d["genus"].sort_values("genomes"); y=np.arange(len(g)); axs[0].barh(y,g.genomes,color=COLORS["ePhaZ"]); axs[0].set_yticks(y,g.genus,fontsize=5.8); axs[0].set_xlabel("Unique genomes"); axs[0].set_title("Leading genera",fontsize=8,fontweight="bold"); style(axs[0]); label(axs[0],"a")
    co=d["co"]; axs[1].bar(co.family_count,co.genomes,color=COLORS["iPhaZ"],width=.65); axs[1].set_xlabel("Families detected per genome"); axs[1].set_ylabel("Genomes"); axs[1].set_title("ePhaZ/iPhaZ co-occurrence",fontsize=8,fontweight="bold"); style(axs[1]); label(axs[1],"b"); fig.suptitle("Genome representation of the two main families",fontsize=9,fontweight="bold"); fig.text(.5,.01,"Taxonomy and co-occurrence are computational census summaries, not ecological or phenotype claims.",ha="center",fontsize=6,color="#666"); fig.subplots_adjust(top=.84,bottom=.17,left=.14,right=.98,wspace=.38); save(fig,out,"figure_5_main_family_genome_representation")

def figure6(d,out):
    c=d["competition"]; majority, minority=competition_display_groups(c); cols={"ePhaZ_like":COLORS["ePhaZ"],"iPhaZ_like":COLORS["iPhaZ"],"ambiguous":COLORS["review"]}
    fig=plt.figure(figsize=(WIDTH_IN,3.7)); outer=fig.add_gridspec(1,2,width_ratios=[.9,1.1],wspace=.42); left=outer[0].subgridspec(2,1,height_ratios=[1,.9],hspace=.65); major_ax=fig.add_subplot(left[0]); minor_ax=fig.add_subplot(left[1]); boundary_ax=fig.add_subplot(outer[1])
    major_ax.barh(majority.assignment,majority.proteins,color=[cols[x] for x in majority.assignment],zorder=2); major_ax.set_xlabel("Proteins",fontsize=6); major_ax.set_title("ePhaZ OTHER competition audit",fontsize=8,fontweight="bold"); major_ax.set_xlim(0,int(majority.proteins.max()*1.18)); style(major_ax); label(major_ax,"a")
    for yi, value in enumerate(majority.proteins): major_ax.text(int(value)*1.01,yi,f"{int(value):,}",va="center",fontsize=6)
    minor_ax.barh(minority.assignment,minority.proteins,color=[cols[x] for x in minority.assignment],zorder=2); minor_ax.set_xlabel("Proteins (minority assignments)",fontsize=6); minor_ax.set_xlim(0,max(22,int(minority.proteins.max()*1.25))); minor_ax.invert_yaxis(); style(minor_ax)
    for yi, value in enumerate(minority.proteins): minor_ax.text(int(value)+.5,yi,str(int(value)),va="center",fontsize=6)
    counts=pd.DataFrame({"class":["ePhaZ curated core","iPhaZ tier1","All ePhaZ/iPhaZ"],"proteins":[38692,32926,109087]}); boundary_ax.barh(np.arange(3),counts.proteins,color=[COLORS["ePhaZ"],COLORS["iPhaZ"],COLORS["gray"]]); boundary_ax.set_yticks(np.arange(3),counts['class'],fontsize=6); boundary_ax.set_xlabel("Candidate proteins"); boundary_ax.set_title("Conservative PHB boundary",fontsize=8,fontweight="bold"); style(boundary_ax); label(boundary_ax,"b"); fig.suptitle("Competition evidence is a review layer, not a family reclassification",fontsize=9,fontweight="bold"); fig.text(.5,.025,"Input to competition audit: ePhaZ curated-core OTHER = 16,836; results remain candidate-only and PHB-unresolved.",ha="center",fontsize=6,color="#666"); fig.subplots_adjust(top=.82,bottom=.18,left=.14,right=.98); save(fig,out,"figure_6_competition_and_phb_boundary")

def supp1(d,out):
    lay=d["layers"]; fig,axs=plt.subplots(1,2,figsize=(WIDTH_IN,3.4),gridspec_kw={"width_ratios":[1,.95]}); y=np.arange(len(lay)); axs[0].bar(y,lay.records,color=[LAYER_COLORS[x] for x in lay.layer]); axs[0].set_xticks(y,[LAYER_LABELS[x].replace(" ","\n") for x in lay.layer],fontsize=5.8); axs[0].set_ylabel("Archaeal records"); axs[0].set_title("Independent archaeal side branch (n = 2,307)",fontsize=8,fontweight="bold"); style(axs[0]); label(axs[0],"a")
    piv=d["motif"].pivot(index="layer",columns="nterm_gxsxg_motif",values="records").reindex(["PhaZh1_like_high","PhaZh1_like_review"]).fillna(0); yy=np.arange(len(piv)); axs[1].barh(yy,piv.get("present",0),color="#4C9F70",label="GxSxG present"); axs[1].barh(yy,piv.get("absent",0),left=piv.get("present",0),color=COLORS["light"],label="absent"); axs[1].set_yticks(yy,[LAYER_LABELS[x] for x in piv.index]); axs[1].set_xlabel("Domain-annotated records"); axs[1].set_title("Motif evidence in priority strata",fontsize=8,fontweight="bold"); axs[1].legend(fontsize=5.5); style(axs[1]); label(axs[1],"b"); fig.suptitle("Supplementary Figure 1 | Archaeal PhaZh1-like audit",fontsize=9,fontweight="bold"); fig.text(.5,.01,"This branch is excluded from the bacterial ePhaZ/iPhaZ main conclusion; all labels remain candidate-only.",ha="center",fontsize=6,color="#666"); fig.subplots_adjust(top=.82,bottom=.2,left=.11,right=.98,wspace=.35); save(fig,out,"supplementary_figure_1_archaeal_side_branch")

def supp2(d,out):
    cal=d["cal"].copy(); fig,axs=plt.subplots(1,2,figsize=(WIDTH_IN,3.4),gridspec_kw={"width_ratios":[1,1]}); groups=["positive_control","negative_control"]; colors=["#4C9F70",COLORS["red"]]
    for i,g in enumerate(groups):
        sub=cal[cal.panel==g]; axs[0].scatter(np.full(len(sub),i),sub.bitscore.astype(float),s=28,c=colors[i],edgecolor="white",lw=.4,label=g.replace("_"," ")); axs[0].text(i,0,f"n={len(sub)}",ha="center",va="bottom",fontsize=6)
    axs[0].set_xticks([0,1],["positive","negative"]); axs[0].set_ylabel("ArchPhaZ_patatin bitscore"); axs[0].set_title("Control calibration",fontsize=8,fontweight="bold"); style(axs[0]); label(axs[0],"a")
    m=d["arch"].copy(); marks=m.nearby_markers.fillna("").str.split(";").explode(); counts=marks[marks!=""].value_counts(); axs[1].barh(np.arange(len(counts)),counts.values,color=COLORS["arch"]); axs[1].set_yticks(np.arange(len(counts)),counts.index,fontsize=6); axs[1].set_xlabel("Records with marker"); axs[1].set_title("Archaeal neighborhood markers",fontsize=8,fontweight="bold"); style(axs[1]); label(axs[1],"b"); fig.suptitle("Supplementary Figure 2 | Domain controls and neighborhood context",fontsize=9,fontweight="bold"); fig.text(.5,.01,"Positive/negative controls and local markers assess evidence quality; they do not validate PHB degradation or pathway membership.",ha="center",fontsize=6,color="#666"); fig.subplots_adjust(top=.82,bottom=.2,left=.12,right=.98,wspace=.36); save(fig,out,"supplementary_figure_2_archaeal_controls_neighborhood")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--run-dir",type=Path,required=True); args=ap.parse_args(); run=args.run_dir; out=run/"results"/"figures"; src=run/"results"/"source_data"; inputs=resolve_input_tables(run); write_input_contract(run,inputs); d=build_data(run)
    for name,frame in {"figure_1_workflow.tsv":d["workflow"],"figure_2_family_hierarchy.tsv":d["subtype"],"figure_3_signalp.tsv":d["signal"],"figure_4_phylum.tsv":d["phylum"],"figure_5_genus.tsv":d["genus"],"figure_5_cooccurrence.tsv":d["co"],"figure_6_competition.tsv":d["competition"],"supplementary_1_archaeal_layers.tsv":d["layers"],"supplementary_1_archaeal_motif.tsv":d["motif"],"supplementary_2_calibration.tsv":d["cal"],"supplementary_2_neighborhood.tsv":d["arch"][["layer","nearby_markers"]]}.items(): write(frame,src/name)
    for fn in (figure1,figure2,figure3,figure4,figure5,figure6,supp1,supp2): fn(d,out)
    stems=["figure_1_workflow_and_funnel","figure_2_family_hierarchy","figure_3_ephaz_signalp_scope","figure_4_main_family_phylum_distribution","figure_5_main_family_genome_representation","figure_6_competition_and_phb_boundary","supplementary_figure_1_archaeal_side_branch","supplementary_figure_2_archaeal_controls_neighborhood"]
    manifest={"run_id":run.name,"status":"completed_candidate_only","backend":"Python matplotlib","main_families":list(MAIN_FAMILIES),"signalp_scope":"ePhaZ curated-core tier1 only (38,692)","archaeal_side_branch":True,"figures":stems,"source_data":sorted(x.name for x in src.glob("*.tsv")),"phenotype_boundary":"All HMM, SignalP, domain, motif, neighborhood and taxonomy evidence remains candidate-only; archaeal branch excluded from main family conclusion."}
    manifest["exports"]=[{"file":x.name,"bytes":x.stat().st_size,"sha256":hashlib.sha256(x.read_bytes()).hexdigest()} for x in sorted(out.iterdir()) if x.is_file()]; (run/"results"/"figure_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
if __name__=="__main__": main()
