# T141 Formal Frozen Scan 11: Failed Run

- Run ID: `20260831_formal_frozen_scan_11`
- Deploy: `deploy/20260831_formal_frozen_preflight_10/`
- Status: failed and terminated; retained as failure evidence

HMMER 3.4 aborted on targets longer than 100,000 aa (`Target sequence length > 100K, over comparison pipeline limit`) in ePhaZ jobs. Five failed tasks were recorded before termination. The run was not accepted as a formal result and was not resumed or overwritten.

The failure is a tool input-limit issue, not a biological negative. A replacement run uses a new run ID and records excluded overlength accessions separately.
