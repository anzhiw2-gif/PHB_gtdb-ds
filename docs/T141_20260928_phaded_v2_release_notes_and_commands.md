# v2 release: state, the one open choice, and the exact commands

Prepared by the agent; **nothing here has been executed**, because the GitHub CLI is not authenticated and no Zenodo token exists in this environment (`gh auth status` reports no logged-in hosts; `GH_TOKEN`, `GITHUB_TOKEN`, `GH_ENTERPRISE_TOKEN`, `ZENODO_TOKEN`, `ZENODO_ACCESS_TOKEN` are all unset; none of `~/.config/gh/hosts.yml`, `%APPDATA%\GitHub CLI\hosts.yml`, `~/.zenodo`, `~/.config/zenodo` exists).

This document exists so that the remaining work is **one command**, not one investigation.

## Current state (measured)

| item | value |
|---|---|
| annotated tag | `phaded-evidence-model-v2-20260928` → tag object `8694888` → commit **`fc944fe`** |
| pushed to origin | **yes** (`refs/tags/phaded-evidence-model-v2-20260928` and its peeled `fc944fe` are both on the remote) |
| `origin/main` | **`3b1f72e`** |
| commits since the tag | **23** |
| working tree | clean; frozen `runs/`, `results/`, `deploy/` untouched |
| gate at `3b1f72e` | 1,889 tests OK (skipped 1); compileall 0; `git diff --check` 0; `test_public_repo_safety` 7 OK |

The tag's own annotation is a complete, accurate release-candidate note for the state at `fc944fe`, and it should not be rewritten.

## The one choice that needs you

The 23 commits since the tag are not cosmetic — they are the entire structural evidence line plus three corrections. Two honest options:

**Option A — release the existing tag as-is.** The tag already exists and carries good notes; `gh release create` on it publishes exactly what `fc944fe` contains. Simple, and the tag's meaning stays "the v2 machinery, as frozen on 2026-09-28". But the published release would be **23 commits behind `main`**, and none of the structural findings would appear in it.

**Option B — tag the current state and release that.** Create e.g. `phaded-evidence-model-v2.1-20260928` at `3b1f72e` and release it with the notes below. The release then matches what a reader gets from `main` today, and the structural line is included.

I recommend **B**, because several of the 23 commits are *corrections to claims the v2 tag itself makes* — most sharply, the tag says the reference-layer gate "detects 100% of experimentally documented extracellular references but only 46.8% of inherited-label ones" without the follow-on that the gate's stringency is therefore unbacked, and it does not carry the finding that sequence-level competition cannot separate nPHAMCL from its measured confounders. Releasing only the tag would publish a picture the project has since refined. But this is a versioning decision, so it is yours.

## Release notes for Option B (paste into `gh release create`)

> **PhaDED evidence-model v2.1 (2026-09-28)** — the structural evidence line, and three corrections to v2.
>
> Builds on `phaded-evidence-model-v2-20260928` (`fc944fe`). The frozen 2026-09-21 results and the v1 comparison baseline are **unchanged**; `runs/`, `results/` and `deploy/` are untouched.
>
> **The structural stage was never blocked — that was my error, now corrected.** v2 recorded it as blocked on the documented hard constraint that single-sequence structure prediction is unusable for comparison (pLDDT 36.1). The frozen evidence had only ever tested `--msa-mode single_sequence`; ColabFold's default MSA mode was never tried. Measured with MSA mode: **pLDDT 92.6–93.4**, comparable models, five sequences in under an hour on one idle GPU. The constraint is real but scoped to single-sequence mode.
>
> **Structural competition, first result (n=2 pilot):** both probe candidates are structurally closer to the **measured confounder** than to the nPHAMCL reference — including the candidate that won its sequence comparison by **19.3 bits**, the largest anchor margin in the whole 987, which is still nearer Q88N36 by 0.081 TM (0.838 vs 0.757). **The sequence-level `nphamcl_like_supported` label carries no structural warranty.** This is a third independent instrument agreeing with the frozen `hfam_4` verdict, and it adds that the nearest neighbour is a *measured* confounder rather than an arbitrary protein.
>
> **F7 step C, executed on the real 985 MB table:** the corrected E-value scale removes all sequence-hit support from **831** of the 109,087 universe accessions (0.76%; 830 consequential), all from the `tier2` review layer — none from curated/tier1. 24,191 newly rejected rows fall inside the universe (853 consequential + 23,338 redundant); 53 accessions appear in no hit row, reproducing the frozen `109,034/109,087`.
>
> **P5: the deferred pool's tail is not its body.** Scored against the same 11-member panel as the nPHAMCL set, the tail below `1e-30` has a median margin of **−31.50 bits** against the body's **−2.70**, and every tail member hits the anchor at `E<1e-5` while only 92.4% of the body does. The body holds 99.6% of the pool, so a proportional sample would have measured a coin flip and reported it as the answer. The tail's margins reach **−538.9**, far beyond the nPHAMCL set's −19.3.
>
> **P3: the Cys anchor maps to four explicit states** — `pattern_present` 28,271 (94.3%), `substitution` 776, `truncation` 558, `uncertain` 369 — with the caveat that `uncertain` is not a negative, because `C183S` retains activity. Closing a coverage gap inside the task changed the answer: 1,492 pool-external rows had no sequence in the candidate union, and only after adding the pool-external source did the `truncation` state appear at all.
>
> **Two plan premises refuted by measurement.** (1) F14's follow-up claim that adding the frozen geometry table would close a 668-row gap removes **zero** rows: the table's verdicts are already reflected row for row, and **82.3% of the universe has geometry that could not determine a catalytic-domain type** — an intrinsic bound, not a missing input. (2) The gRodon 66 and identity-gate premises were refuted in v2 and remain so.
>
> **The 616-row disposition question is now pinned rather than vague.** Of the 613 assigned `remote_homolog_candidate`, only **140** carry a discovery-layer claim; **472** have no discovery row at all and no profile call, so that label asserts evidence they lack. The six-item vocabulary has **no item for "no admissible homology call"** — a spec gap the 616 exposed. A seventh item is recommended; the catalogue is deliberately not rewritten, because choosing between stretching item 4 and adding item 7 is the operator's call.
>
> **Repository hygiene, measured.** The public-repo identity gate was a complete no-op before v2's fix, so the history was audited: **93 violating blobs across 24 commits, 0 in `HEAD`**, oldest at the initial commit. **No credential is present** — the one apparent PEM hit was my auditor's false positive, a test fixture in the safety module that the working-tree gate deliberately skips. The working tree is clean; the history is not; rewriting it is destructive and has **not** been done.
>
> Candidate-only boundary: nothing here proves a PHB/PHA degradation phenotype. The 987 nPHAMCL candidates remain `function_unresolved`.
>
> Gate: 1,889 tests OK (skipped 1, environment-only); compileall 0; `git diff --check` 0; `test_public_repo_safety` 7 OK; frozen trees untouched.

## Exact commands

### Option A — release the existing tag

```bash
gh release create phaded-evidence-model-v2-20260928 \
  --title "PhaDED evidence-model v2 (2026-09-28)" \
  --notes-file docs/T141_20260928_phaded_v2_release_notes_and_commands.md \
  --verify-tag
```
(Point `--notes-file` at the tag's own annotation instead if you prefer: `git tag -l --format='%(contents)' phaded-evidence-model-v2-20260928 > /tmp/notes.md`.)

### Option B — tag the current state and release it

```bash
git tag -a phaded-evidence-model-v2.1-20260928 3b1f72e \
  -m "PhaDED evidence-model v2.1 (2026-09-28): structural evidence line + three corrections"
git push origin phaded-evidence-model-v2.1-20260928

gh release create phaded-evidence-model-v2.1-20260928 \
  --title "PhaDED evidence-model v2.1 (2026-09-28)" \
  --notes-file <the Option B notes above, saved to a file> \
  --verify-tag
```

### Zenodo DOI

Zenodo's GitHub integration is enabled **per repository, in the Zenodo web UI**, and it mints the DOI on release — there is no API call that substitutes for it. If the integration is not already on: log in to Zenodo → Settings → GitHub → flip the toggle for this repository, then publish the release. If you would rather deposit manually, Zenodo's API needs a personal token (`ZENODO_TOKEN`), which does not exist in this environment.

Once a DOI exists, record it in `docs/T141_20260917_project_handoff.md` §12 and in the README; **do not fabricate a DOI placeholder** in the meantime.

## What remains after that

Nothing in F1–F17 except the two items this document cannot reach: the credentials above, and the structural survey currently running on GPU 1 (8 of 30 done at the time of writing), whose results will be reported when it completes.
