#!/usr/bin/env python3
"""Audit discovery-layer / with-lipase deferred-layer usage across ``pipeline/scripts``.

Governance rules under audit
----------------------------

``R1_recall_only_only``
    The PhaDED discovery layer (``model_layer="discovery_hmm_uncalibrated"``) is for
    **recall only**.  It must never filter, delete or demote a candidate.  (Its measured
    behaviour on the 563 x1-confounder hits —443 hits = 78.69 % —is why it has no
    discriminating power.)

``R2_no_family_call``
    The discovery layer must never produce a family/superfamily/subtype judgement, and
    every one of its outputs must carry the discovery-layer label.

``R3_registry_exclusion``
    The discovery layer must not be written into ``pipeline/config/formal_scan_models.tsv``.

``R4_deferred_preserved``
    The with-lipase (``DED_hfam_2``) pool-external hits (1,206,655; archive
    ``runs/20260919_phaded_with_lipase_deferred_archive_01``) must be preserved as a
    deferred structural-validation layer and must never be deleted.

``R5_deferred_out_of_counts``
    The deferred layer must never enter a candidate or high-confidence count.

``R6_reference_query_only_unscoreable``
    A ``reference_query_only`` profile has no HMM at all and must never be scored as if
    it had one.

``R7_layer_status_independent``
    ``model_layer`` and ``functional_calibration_status`` are independent fields: neither
    may be derived from the other.

Method (and why it is not grep)
-------------------------------

Every ``pipeline/scripts/*.py`` is parsed with :mod:`ast`; findings are attached to real
code paths only.  String literals that sit inside a docstring, and every comment, are
invisible to the AST, so a *commented-out* ``if discovery_score > 0.5: drop(...)`` can
never be reported as a violation —:func:`scan_text_signals` records such prose matches
separately (they are reported as ``unreviewed`` prose, never as violations).  A
documented text scan is used for two things the AST cannot see: config/path literals
(``formal_scan_models.tsv``, the with-lipase deferred archive) and ``*.sh`` scripts,
which have no AST at all.

Values are tracked with a small taint analysis over four provenance kinds —
``layer`` (something that identifies the discovery layer), ``score`` (a discovery-layer
score/rank/E-value/family field), ``family`` (a discovery-derived family value),
``status`` (a functional-calibration status).  A tainted value used in a decision
(``if`` / comprehension ``if``) is classified by the *effect* of the guarded block:

* drops/keeps rows        -> ``thresholding`` (R1) or ``family_call`` (R2) when it assigns
* assigns a family field  -> ``family_call`` (R2)
* counts rows             -> ``counting`` (R1/R5)
* ranks/scores/labels only-> ``recall_only`` (compliant)
* undecidable effect      -> ``unknown`` / ``unreviewed`` (never silently ``false``)

The audit never executes, imports or modifies a scanned script, and never writes into
``runs/``, ``results/`` or ``deploy/``.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Vocabulary
# ---------------------------------------------------------------------------

REQUIRED_COLUMNS = (
    "script", "line", "evidence_kind", "signal", "use_class", "violation",
    "rule_cited", "note",
)

USE_CLASSES = (
    "recall_only",       # scores/ranks/labels only, never filters
    "thresholding",      # includes/excludes without a validation gate
    "family_call",       # writes a family/superfamily/subtype/status judgement
    "counting",          # the governance layer text mentions it in a count
    "registry_write",    # references pipeline/config/formal_scan_models.tsv
    "deferred_layer_touch",  # reads/writes the with-lipase deferred archive
    "unknown",           # nothing decidable from static source
)

VIOLATION_VALUES = ("true", "false", "unreviewed")

RULES = {
    "R1_recall_only_only":
        "AGENTS.md: \u53d1\u73b0\u5c42\u53ea\u7528\u4e8e\u53ec\u56de\uff1b"
        "\u4e0d\u5f97\u7528\u4e8e\u7b5b\u9009\u3001\u5220\u9664\u6216\u964d\u7ea7"
        "\u4efb\u4f55\u5019\u9009 (discovery layer is recall-only; it must never filter, "
        "delete or demote a candidate)",
    "R2_no_family_call":
        "AGENTS.md: \u53d1\u73b0\u5c42\u4e0d\u5f97\u4ea7\u751f family \u5224\u5b9a "
        "(the discovery layer must never produce a family/superfamily/subtype call; every "
        "output carries its layer label)",
    "R3_registry_exclusion":
        "AGENTS.md: \u53d1\u73b0\u5c42\u4e0d\u5f97\u5199\u5165 "
        "pipeline/config/formal_scan_models.tsv",
    "R4_deferred_preserved":
        "AGENTS.md: with-lipase (DED_hfam_2) \u6c60\u5916\u547d\u4e2d\u5fc5\u987b"
        "\u4fdd\u7559\u4e3a deferred \u7ed3\u6784\u9a8c\u8bc1\u5c42\uff0c\u7981\u5220",
    "R5_deferred_out_of_counts":
        "AGENTS.md: deferred \u5c42\u4e0d\u8fdb\u5165\u4efb\u4f55\u5019\u9009/"
        "\u9ad8\u53ef\u4fe1\u5ea6\u8ba1\u6570",
    "R6_reference_query_only_unscoreable":
        "AGENTS.md: reference_query_only \u65e0 HMM\uff0c\u4e0d\u5f97\u88ab\u5f53\u4f5c"
        "\u6709 HMM \u6253\u5206",
    "R7_layer_status_independent":
        "AGENTS.md: model_layer \u4e0e functional_calibration_status \u72ec\u7acb\uff0c"
        "\u4e0d\u5f97\u4e92\u76f8\u63a8\u5bfc",
    "NONE":
        "informational: no governance rule is broken by this row",
}

#: ``use_class`` for R7 findings is ``family_call``: manufacturing a verdict from
#: ``model_layer`` is the same failure shape as manufacturing a family call from the
#: discovery layer, and the task fixes this vocabulary.  The precise construct is always
#: named in ``signal`` and ``note``.

AUDIT_LIMITS = (
    "No execution: no scanned script is imported, executed or traced; every judgement "
    "comes from source text and the AST.",
    "Dynamic dispatch is invisible: getattr/setattr, importlib, string-built module "
    "names, dict-driven callables, plugin registries and subprocess argv assembled at "
    "runtime cannot be followed, so a violation reached only through dynamic dispatch "
    "will be reported as recall_only or unknown rather than as a violation.",
    "Data-dependent branching is invisible: whether a table read from a path holds "
    "discovery-layer rows, whether a manifest row is discovery-layer, and which branch "
    "actually runs on real data are runtime facts.  A construct whose discovery "
    "provenance arrives only through data (for example a value read from a table whose "
    "path does not name the layer) is reported unknown/unreviewed, never as false.",
    "Non-Python entry points are text-scanned only: *.sh files have no AST here, so "
    "their rows use the documented text heuristics and any judgement beyond a literal "
    "path reference is unreviewed.",
    "Shell/Python-interop and code generated at runtime (heredocs, eval, exec, "
    "generated modules) are not analysed.",
    "Absence of a violation row is not a proof of compliance: it means the audit could "
    "not see the misuse, not that the misuse does not exist.",
)


# ---------------------------------------------------------------------------
# Tokens
# ---------------------------------------------------------------------------

MODEL_LAYERS = (
    "reference_query_only",
    "discovery_hmm_uncalibrated",
    "sequence_family_hmm_validated",
    "calibrated_candidate_model",
)
FUNCTIONAL_CALIBRATION_STATUSES = (
    "calibrated_candidate_model",
    "not_function_calibrated",
    "blocked_contradictory_evidence",
)

#: Exact discovery-layer identifiers.  The bare word ``discovery`` is deliberately NOT a
#: token: ``ePhaZ_broad_discovery`` is an ePhaZ model in the formal-scan registry and has
#: nothing to do with the PhaDED discovery layer.
DISCOVERY_LAYER_LITERALS = (
    "discovery_hmm_uncalibrated",
    "cys_discovery_uncalibrated",
)
DISCOVERY_LAYER_TOKEN_RE = re.compile(
    "|".join(re.escape(token) for token in DISCOVERY_LAYER_LITERALS)
)

#: Identifiers/keys that *may* carry discovery provenance.  ``disc`` is included only with
#: a word boundary (``disc_e``), never as a prefix of ``discard``/``discover``.
DISCOVERY_NAME_RE = re.compile(r"discovery|uncalibrated|(^|_)disc(_|$)", re.I)
#: Names that contain ``discovery`` but have nothing to do with the PhaDED discovery
#: layer: ``ePhaZ_broad_discovery`` is an ePhaZ family model in the formal-scan registry.
NOT_DISCOVERY_LAYER_RE = re.compile(r"broad_discovery|ephaZ_broad", re.I)
SCORE_HINT_RE = re.compile(
    r"score|evalue|e_value|bits|rank|pvalue|p_value|threshold|thresh|count"
    r"|(^|_)hits?(_|$)|(^|_)cov(erage)?(_|$)",
    re.I,
)
FAMILY_HINT_RE = re.compile(
    r"family|superfamily|subtype|clan|confiden|clade|(^|_)calls?(_|$)", re.I
)
CALIBRATION_HINT_RE = re.compile(r"calibration_status|functional_calibration|calibrated", re.I)
LAYER_FIELD_RE = re.compile(r"^(model_)?layer$|_layer$", re.I)


def _discovery_name(text: str) -> bool:
    """Whether a name/key/literal refers to the PhaDED discovery layer."""
    return bool(DISCOVERY_NAME_RE.search(text)) and not NOT_DISCOVERY_LAYER_RE.search(text)

REGISTRY_TOKEN_RE = re.compile(r"formal_scan_models\.tsv|FORMAL_SCAN_REGISTRY", re.I)
DEFERRED_ARCHIVE_TOKEN_RE = re.compile(
    r"DED_hfam_2"
    r"|pool_external_with_lipase_deferred"
    r"|with_lipase_deferred"
    r"|deferred_structural_validation"
    r"|20260919_phaded_with_lipase_deferred_archive_01",
    re.I,
)

FAMILY_TARGET_RE = re.compile(
    r"family|superfamily|subtype|clan|(^|_)calls?(_|$)|verdict|disposition|confiden", re.I
)
#: A target that *is* a family call rather than a family-shaped local: only these turn a
#: discovery-derived value into a proven family call at assignment level.
STRONG_CALL_TARGET_RE = re.compile(
    r"(^|_)(family_call|superfamily_call|subtype_call|call|family|superfamily|subtype)s?$"
    r"|sequence_family$|primary_disposition$",
    re.I,
)
#: A module-level constant that merely *names* a rule/reason/column is not a call.
CONSTANT_TARGET_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")
COUNT_TARGET_RE = re.compile(r"count|_n$|(^|_)num|total", re.I)
CANDIDATE_COUNT_RE = re.compile(
    r"candidates?|high_confidence|high-confiden|(^|_)core(_|$)|promoted|primary",
    re.I,
)
SCORE_TARGET_RE = re.compile(r"score|evalue|bits|rank", re.I)
DROP_CALL_RE = re.compile(
    r"^(drop|exclude|reject|remove|discard|filter_out|delete|purge|demote|skip|blacklist|prune)",
    re.I,
)
DELETE_CALL_RE = re.compile(r"(unlink|rmtree|remove|rmdir|truncate|shred)", re.I)
SCORING_CALL_RE = re.compile(r"hmmsearch|hmmscan|hmmbuild|pyhmmer|tblout|build_hmmsearch_command", re.I)
ASSIGN_CALL_RE = re.compile(r"assign|classify|resolve|call_", re.I)
DROP_TARGET_RE = re.compile(r"excluded|dropped|rejected|removed|discarded|skipped|pruned|demoted", re.I)
TRAINING_HINT_RE = re.compile(r"training|train_|eligible|discovery_set", re.I)
WEAK_VALUES = frozenset(
    {"", "low", "unresolved", "pending", "unknown", "none", "false", "not_tested",
     "not_function_calibrated", "candidate_only", "remote_homolog_candidate",
     "low_candidate_only", "reference_only_nearest", "moderate_candidate_only"}
)

SHELL_DELETE_RE = re.compile(r"(^|[;&|]\s*)(rm|shred|truncate)\s")

COMPARISON_RE = re.compile(r"<=|>=|==|!=|<|>")
THRESHOLD_RE = re.compile(r"[<>]=?\s*[-+]?[\d.]|==\s*[\"']?(?:discovery|uncalibrated|reference_only|trained)")

PROSE_VERBS_RE = re.compile(
    r"filter|exclude|delete|drop|demote|discard|remove|never enters|not enter|recall[- ]only|"
    r"只用于召回|不得|禁删|不进入|不筛选|降级",
    re.I,
)
PROSE_FAMILY_RE = re.compile(r"family call|superfamily|subtype|family assignment|family 判定", re.I)

SELF_EXCLUSION_REASON = (
    "the audit excludes its own file: it contains the governance token tables and regexes "
    "this scan looks for, so self-scanning would report the detector, not a use"
)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

_REPO_ROOT = Path(__file__).resolve().parents[2]


class AuditError(Exception):
    """Usage/validation error: the audit refuses to run rather than guess."""


def _display_path(path: object) -> str:
    """Return a repo-relative POSIX path, or the caller's label verbatim."""
    text = str(path)
    candidate = Path(text)
    root = _REPO_ROOT.resolve()
    try:
        resolved = candidate.resolve() if candidate.is_absolute() else (root / candidate).resolve()
        return resolved.relative_to(root).as_posix()
    except (ValueError, OSError):
        return text.replace("\\", "/")


def _finding(script, line, evidence_kind, signal, use_class, violation, rule_cited, note):
    return {
        "script": script,
        "line": int(line),
        "evidence_kind": evidence_kind,
        "signal": signal,
        "use_class": use_class,
        "violation": violation,
        "rule_cited": rule_cited,
        "note": note,
    }


def _sort_findings(findings):
    return sorted(
        findings,
        key=lambda row: (row["script"], row["line"], row["signal"], row["use_class"]),
    )


def _dedupe(findings):
    seen = set()
    out = []
    for row in findings:
        key = (row["script"], row["line"], row["evidence_kind"], row["signal"], row["use_class"])
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def _call_name(func) -> str:
    if isinstance(func, ast.Name):
        return func.id
    if isinstance(func, ast.Attribute):
        return func.attr
    return ""


def _target_names(node) -> list:
    if isinstance(node, ast.Name):
        return [node.id]
    if isinstance(node, ast.Attribute):
        return [node.attr]
    if isinstance(node, (ast.Tuple, ast.List)):
        names = []
        for element in node.elts:
            names.extend(_target_names(element))
        return names
    if isinstance(node, ast.Subscript):
        return []
    return []


def _assign_pair(node):
    """Return ``(target_names, value_node)`` for an assignment statement."""
    if isinstance(node, ast.Assign):
        names = []
        for target in node.targets:
            names.extend(_target_names(target))
        return names, node.value
    if isinstance(node, ast.AnnAssign):
        return _target_names(node.target), node.value
    if isinstance(node, ast.AugAssign):
        return _target_names(node.target), node.value
    return [], None


def _string_constants(node) -> list:
    out = []
    if node is None:
        return out
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return [node.value]
    for child in ast.iter_child_nodes(node):
        out.extend(_string_constants(child))
    return out


def _dict_string_keys(node) -> list:
    return [
        key.value for key in node.keys
        if isinstance(key, ast.Constant) and isinstance(key.value, str)
    ]


def _source_names(node) -> str:
    """Every identifier mentioned in a test, as one space-joined string."""
    return " ".join(
        child.id for child in ast.walk(node) if isinstance(child, ast.Name)
    )


def _is_empty_constant(node) -> bool:
    if isinstance(node, ast.Constant):
        if node.value is None or node.value is False:
            return True
        if isinstance(node.value, str) and not node.value.strip():
            return True
    if isinstance(node, (ast.Dict, ast.List, ast.Tuple, ast.Set)) and not len(node.elts if hasattr(node, "elts") else node.keys):
        return True
    if isinstance(node, ast.List) and not node.elts:
        return True
    return False


def _is_weak_value(node) -> bool:
    if _is_empty_constant(node):
        return True
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value.strip().lower() in WEAK_VALUES
    return False


def _docstring_lines(tree) -> set:
    lines = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = getattr(node, "body", [])
            if not body:
                continue
            first = body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) \
                    and isinstance(first.value.value, str):
                end = getattr(first, "end_lineno", first.lineno)
                lines.update(range(first.lineno, end + 1))
    return lines


# ---------------------------------------------------------------------------
# Taint
# ---------------------------------------------------------------------------

def _name_kinds(name: str) -> set:
    if not name:
        return set()
    lowered = name.lower()
    kinds = set()
    if CALIBRATION_HINT_RE.search(lowered):
        kinds.add("status")
        return kinds
    if LAYER_FIELD_RE.search(lowered) and lowered not in {"layer", "layers", "model_layer"}:
        # `*_layer` names are layer-typed but not necessarily the discovery layer; the
        # `layer`/`model_layer` field itself is only a *carrier* and is not tainted.
        kinds.add("field")
    if lowered in {"layer", "model_layer", "model_layers"}:
        kinds.add("field")
    if _discovery_name(lowered):
        if SCORE_HINT_RE.search(lowered):
            kinds.add("score")
        elif FAMILY_HINT_RE.search(lowered):
            kinds.add("score")
            kinds.add("family")
        else:
            kinds.add("layer")
    return kinds


def _literal_kinds(value: str) -> set:
    lowered = value.lower()
    kinds = set()
    if DISCOVERY_LAYER_TOKEN_RE.search(lowered):
        kinds.add("layer")
    elif _discovery_name(lowered) and " " not in lowered.strip():
        # A single-token discovery identifier (a table name, field name or dotted path).
        # Prose that merely contains the word "discovery" is not provenance.
        kinds.add("score")
    if lowered in FUNCTIONAL_CALIBRATION_STATUSES:
        kinds.add("status")
    if lowered in MODEL_LAYERS:
        kinds.add("field")
    return kinds


def _key_kinds(key: str) -> set:
    lowered = key.strip().lower()
    kinds = set()
    if lowered in FUNCTIONAL_CALIBRATION_STATUSES or CALIBRATION_HINT_RE.search(lowered):
        kinds.add("status")
    if lowered in {"model_layer", "layer"}:
        kinds.add("field")
    if _discovery_name(lowered):
        if FAMILY_HINT_RE.search(lowered):
            kinds.add("score")
            kinds.add("family")
        elif SCORE_HINT_RE.search(lowered):
            kinds.add("score")
        else:
            kinds.add("layer")
    return kinds


def _kinds_of(node, env) -> frozenset:
    """Provenance kinds flowing into an expression."""
    kinds = set()
    if node is None:
        return frozenset()
    if isinstance(node, ast.Name):
        kinds |= env.get(node.id, set())
    elif isinstance(node, ast.Attribute):
        kinds |= _kinds_of(node.value, env)
        kinds |= _name_kinds(node.attr)
    elif isinstance(node, ast.Subscript):
        kinds |= _kinds_of(node.value, env)
        for key in _string_constants(node.slice):
            kinds |= _key_kinds(key)
    elif isinstance(node, ast.Constant):
        if isinstance(node.value, str):
            kinds |= _literal_kinds(node.value)
    elif isinstance(node, ast.Call):
        kinds |= _name_kinds(_call_name(node.func))
        kinds |= _kinds_of(node.func, env)
        for arg in node.args:
            kinds |= _kinds_of(arg, env)
        for keyword in node.keywords:
            kinds |= _kinds_of(keyword.value, env)
    elif isinstance(node, ast.Lambda):
        kinds |= _kinds_of(node.body, env)
    elif isinstance(node, ast.Compare):
        kinds |= _kinds_of(node.left, env)
        for comparator in node.comparators:
            kinds |= _kinds_of(comparator, env)
    elif isinstance(node, ast.BoolOp):
        for value in node.values:
            kinds |= _kinds_of(value, env)
    elif isinstance(node, ast.UnaryOp):
        kinds |= _kinds_of(node.operand, env)
    elif isinstance(node, ast.BinOp):
        kinds |= _kinds_of(node.left, env) | _kinds_of(node.right, env)
    elif isinstance(node, ast.IfExp):
        kinds |= _kinds_of(node.test, env)
        kinds |= _kinds_of(node.body, env)
        kinds |= _kinds_of(node.orelse, env)
    elif isinstance(node, ast.Dict):
        # Container literals do not propagate element taint onto the container name: a
        # module-level table that happens to *hold* a discovery-named member (an ePhaZ
        # family map, a scenario table) must not taint every later use of that table.
        # What a container *emits* is classified field-wise by ``visit_Dict``.
        for key in node.keys:
            kinds |= _kinds_of(key, env)
    elif isinstance(node, (ast.Tuple, ast.List, ast.Set)):
        pass
    elif isinstance(node, ast.JoinedStr):
        for value in node.values:
            kinds |= _kinds_of(value, env)
    elif isinstance(node, ast.FormattedValue):
        kinds |= _kinds_of(node.value, env)
    else:
        for child in ast.iter_child_nodes(node):
            kinds |= _kinds_of(child, env)
    return frozenset(kinds)


def _index_parents(tree) -> dict:
    parents: dict = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parents[child] = node
    return parents


def _enclosing_function(parents, node):
    current = parents.get(node)
    while current is not None:
        if isinstance(current, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
            return current
        current = parents.get(current)
    return None


def _arg_names(func) -> list:
    args = getattr(func, "args", None)
    if args is None:
        return []
    names = []
    for group in ("posonlyargs", "args", "kwonlyargs"):
        names.extend(argument.arg for argument in getattr(args, group, []) or [])
    return names


def _propagate_env(base: dict, assignments, depth: int = 3) -> dict:
    env = {name: set(kinds) for name, kinds in base.items()}
    for _ in range(depth):
        for node in assignments:
            names, value = _assign_pair(node)
            if not names or value is None:
                continue
            kinds = set(_kinds_of(value, env))
            for name in names:
                kinds |= _name_kinds(name)
            if not kinds:
                continue
            for name in names:
                env.setdefault(name, set()).update(kinds)
    return env


def _build_scoped_taint_env(tree):
    """Return ``(module_env, {function_node: env})``.

    Taint is scoped: a discovery-derived value assigned inside one function never
    taints a same-named local of another function.  File-global name propagation is what
    turns a 1,700-line script into hundreds of false "discovery decision" rows.
    """
    parents = _index_parents(tree)
    module_assignments = []
    function_assignments: dict = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            function_assignments.setdefault(node, [])
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            continue
        owner = _enclosing_function(parents, node)
        if owner is None:
            module_assignments.append(node)
        elif owner in function_assignments:
            function_assignments[owner].append(node)
    module_env = _propagate_env({}, module_assignments)
    function_envs = {}
    for func, assignments in function_assignments.items():
        seed = {name: set(_name_kinds(name)) for name in _arg_names(func)}
        base = {name: set(kinds) for name, kinds in module_env.items()}
        for name, kinds in seed.items():
            base.setdefault(name, set()).update(kinds)
        function_envs[func] = _propagate_env(base, assignments)
    return module_env, function_envs, parents


# ---------------------------------------------------------------------------
# Effects of a guarded block
# ---------------------------------------------------------------------------

def _effects(nodes) -> set:
    effects = set()
    collect_lines = []
    skip_lines = []
    for statement in nodes:
        for node in ast.walk(statement):
            if isinstance(node, (ast.Continue, ast.Break, ast.Delete)):
                effects.add("drop")
                skip_lines.append(node.lineno)
            elif isinstance(node, ast.Raise):
                effects.add("raise")
            elif isinstance(node, ast.Return):
                effects.add("return")
                if node.value is None or _is_empty_constant(node.value):
                    effects.add("suppress")
                elif _is_weak_value(node.value):
                    effects.add("weaken")
                else:
                    effects.add("return_value")
            elif isinstance(node, ast.Call):
                name = _call_name(node.func)
                if name and DROP_CALL_RE.search(name):
                    effects.add("drop")
                if name and DELETE_CALL_RE.search(name):
                    effects.add("delete")
                if name and SCORING_CALL_RE.search(name):
                    effects.add("score")
                if name and ASSIGN_CALL_RE.search(name):
                    effects.add("assign")
                if name in {"len", "sum", "Counter"}:
                    effects.add("count")
                if name in {"append", "extend", "add", "update", "insert", "setdefault"}:
                    effects.add("collect")
                    collect_lines.append(node.lineno)
            elif isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for target_node in targets:
                    for target in _target_names(target_node):
                        if DROP_TARGET_RE.search(target):
                            effects.add("drop")
                        if FAMILY_TARGET_RE.search(target) and not COUNT_TARGET_RE.search(target):
                            effects.add("assign")
                        if COUNT_TARGET_RE.search(target):
                            effects.add("count")
                            if CANDIDATE_COUNT_RE.search(target):
                                effects.add("candidate_count")
                        if SCORE_TARGET_RE.search(target):
                            effects.add("score")
                if _is_empty_constant(getattr(node, "value", None)):
                    effects.add("suppress")
                elif _is_weak_value(getattr(node, "value", None)):
                    effects.add("weaken")
            elif isinstance(node, ast.Dict):
                keys = _dict_string_keys(node)
                for key in keys:
                    if FAMILY_TARGET_RE.search(key):
                        effects.add("assign")
                    if COUNT_TARGET_RE.search(key):
                        effects.add("count")
                        if CANDIDATE_COUNT_RE.search(key):
                            effects.add("candidate_count")
                    if SCORE_TARGET_RE.search(key):
                        effects.add("score")
                    if LAYER_FIELD_RE.search(key) or key == "model_layer":
                        effects.add("label")
    if collect_lines and skip_lines and min(collect_lines) < min(skip_lines):
        # `if <tainted>: emit a labelled row ... ; continue` retains the row: the block
        # labels candidates rather than filtering them out.
        effects.add("emit_then_skip")
    return effects


# ---------------------------------------------------------------------------
# AST findings
# ---------------------------------------------------------------------------

def _has_discovery_label(nodes, env) -> bool:
    """Whether an emitted row/assignment is explicitly labelled as the discovery layer."""
    for statement in nodes:
        for node in ast.walk(statement):
            if isinstance(node, ast.Dict):
                for key, value in zip(node.keys, node.values):
                    if not (isinstance(key, ast.Constant) and key.value == "model_layer"):
                        continue
                    if any(DISCOVERY_LAYER_TOKEN_RE.search(item) for item in _string_constants(value)):
                        return True
                    if "layer" in _kinds_of(value, env):
                        return True
            names, value = _assign_pair(node)
            if "model_layer" in names and value is not None:
                if any(DISCOVERY_LAYER_TOKEN_RE.search(item) for item in _string_constants(value)):
                    return True
                if "layer" in _kinds_of(value, env):
                    return True
    return False


def _decision_result(kinds, effects, labelled: bool, context: str):
    """Map (provenance kinds, block effects) to ``(use_class, violation, rule, note)``."""
    data_tainted = bool(kinds & {"score", "family"})
    layer_gate = bool(kinds & {"layer"}) and not data_tainted

    if data_tainted:
        if "training_selection" in effects:
            return ("recall_only", "false", "R1_recall_only_only",
                    "selects discovery-layer *training* rows (AGENTS.md explicitly permits "
                    "training a discovery HMM from DED family sequences); no candidate is "
                    "filtered, deleted or demoted")
        if "drop" in effects or "delete" in effects or "select" in effects:
            if "score_threshold" not in effects:
                return ("unknown", "unreviewed", "R1_recall_only_only",
                        "a value that can carry discovery-layer provenance gates which rows "
                        "are kept or dropped, but no discovery score is compared against a "
                        "threshold here: the audit cannot tell whether this filters "
                        "candidates, labels them or merely cites them")
            if "emit_then_skip" in effects:
                return ("thresholding", "unreviewed", "R1_recall_only_only",
                        "the block emits a labelled row and then skips the remaining work "
                        "(so the row is retained, not filtered), but whether the collection "
                        "it feeds is a labelled recall/hit table or a candidate set cannot "
                        "be decided from this file alone")
            if "candidate_set" in effects:
                return ("thresholding", "true", "R1_recall_only_only",
                        "a discovery-layer score/family value decides which candidates are "
                        "kept: the filtered collection is named as a candidate/"
                        "high-confidence set, and R1 forbids the discovery layer from "
                        "filtering, deleting or demoting candidates")
            return ("thresholding", "unreviewed", "R1_recall_only_only",
                    "a discovery-layer score decides which rows are kept or dropped, but "
                    "the collection is not named as a candidate/high-confidence set: this is "
                    "either the permitted recall/hit-table use or a forbidden candidate "
                    "filter, and static source cannot tell which -- a human must follow the "
                    "collection to its consumer")
        if "assign" in effects:
            if labelled:
                return ("family_call", "unreviewed", "R2_no_family_call",
                        "a discovery-layer value selects the emitted family/superfamily/"
                        "subtype value; the row is explicitly labelled with the discovery "
                        "layer, so whether a labelled discovery-layer assignment is still a "
                        "'family call' under R2 needs a human decision")
            return ("family_call", "true", "R2_no_family_call",
                    "a discovery-layer score/family value is written into a family/"
                    "superfamily/subtype/confidence field; R2 forbids the discovery layer "
                    "from producing a family call")
        if "candidate_count" in effects:
            return ("counting", "true", "R1_recall_only_only",
                    "a discovery-layer score/family value feeds a candidate or "
                    "high-confidence count; R1 forbids the discovery layer from being used "
                    "to decide the candidate set")
        if "count" in effects:
            return ("counting", "false", "R1_recall_only_only",
                    "a discovery-layer value feeds a count that is not named as a "
                    "candidate/high-confidence count")
        if "collect" in effects:
            if "score_threshold" not in effects:
                return ("unknown", "unreviewed", "R1_recall_only_only",
                        "the block collects or annotates rows selected by a value that can "
                        "carry discovery-layer provenance, without comparing a discovery "
                        "score against a threshold; the audit cannot tell whether this is a "
                        "recall/annotation list or a candidate filter")
            return ("thresholding", "unreviewed", "R1_recall_only_only",
                    "rows are collected through a discovery-layer score threshold "
                    "(inclusion without a validation gate); if the collected list is a "
                    "candidate/high-confidence set this breaks R1, if it is a recall/hit "
                    "list it is the permitted recall use -- the audit cannot tell which")
        if effects & {"score", "label", "rank", "weaken"}:
            return ("recall_only", "false", "R1_recall_only_only",
                    "the discovery-layer value is only scored, ranked, labelled or "
                    "downgraded; no candidate is included, excluded or deleted")
        return ("unknown", "unreviewed", "R1_recall_only_only",
                "a discovery-layer value drives a decision whose effect this audit cannot "
                "classify; a human must read the block")

    if layer_gate:
        if "drop" in effects or "delete" in effects:
            if "candidate_set" in effects:
                return ("thresholding", "true", "R1_recall_only_only",
                        "rows are dropped by a discovery-layer *label* test while building a "
                        "candidate/high-confidence set; R1 forbids the discovery layer from "
                        "deciding the candidate set")
            return ("thresholding", "unreviewed", "R1_recall_only_only",
                    "rows are dropped/deleted by a discovery-layer *label* test; R1 forbids "
                    "the discovery layer from filtering or deleting candidates, but a label "
                    "gate can also be boundary enforcement, so the effect needs review")
        if "raise" in effects:
            return ("recall_only", "false", "R2_no_family_call",
                    "fails closed on a row that violates the discovery-layer boundary "
                    "(enforcement of R1/R2), it does not filter on a score")
        if "suppress" in effects or "weaken" in effects:
            return ("recall_only", "false", "R2_no_family_call",
                    "a discovery-layer row is emitted with an empty/weakest value: the "
                    "discovery layer makes no family call and cannot support a stronger "
                    "claim")
        if effects & {"score", "label", "rank"} or ("select" in effects and "collect" in effects):
            return ("recall_only", "false", "R1_recall_only_only",
                    "the discovery-layer label only routes scoring/labelling/recall "
                    "collection; no candidate is filtered, deleted or demoted")
        if "collect" in effects or "select" in effects:
            return ("recall_only", "false", "R1_recall_only_only",
                    "rows are routed into (or collected from) a discovery-layer labelled "
                    "recall set; R1 permits recall, so this is compliant as long as the set "
                    "is not consumed as a candidate set")
        return ("unknown", "unreviewed", "R1_recall_only_only",
                "a discovery-layer label drives a decision whose effect this audit cannot "
                "classify; a human must read the block")

    if not (kinds & {"score", "family", "layer"}):
        # Only the model-layer / calibration field is compared: this is a layer-aware
        # decision, not a discovery-layer use.  Whether the discovery layer flows through
        # one of its branches is a runtime manifest fact.
        return ("unknown", "unreviewed", "NONE",
                "the block tests the model-layer / functional-calibration field "
                f"({'/'.join(sorted(kinds))}) but does not name the discovery layer; whether "
                "the discovery layer reaches this branch (and therefore whether any rule is "
                "engaged) depends on runtime manifest values, so no rule is cited")

    return ("unknown", "unreviewed", "R1_recall_only_only",
            f"a discovery/deferred-derived value drives a decision ({context}) that this "
            "audit cannot classify; the provenance arrives only through data")


def _test_compares_threshold(test, env) -> bool:
    """Whether the test compares a discovery score against an ordering threshold.

    ``if discovery_score > 0.5`` is a threshold; ``if discovery_claim:`` or
    ``if x not in discovery_members`` is a flag/membership test and must not be reported
    as score thresholding.
    """
    for node in ast.walk(test):
        if not isinstance(node, ast.Compare):
            continue
        if not any(isinstance(op, (ast.Lt, ast.LtE, ast.Gt, ast.GtE)) for op in node.ops):
            continue
        operands = [node.left, *node.comparators]
        if any(_kinds_of(operand, env) & {"score", "family"} for operand in operands):
            return True
    return False


def _comparison_is_exclusion(test) -> bool:
    for node in ast.walk(test):
        if isinstance(node, ast.Compare):
            for operator in node.ops:
                if isinstance(operator, (ast.NotIn, ast.NotEq)):
                    return True
    return False


def _comparison_is_inclusion(test) -> bool:
    for node in ast.walk(test):
        if isinstance(node, ast.Compare):
            for operator in node.ops:
                if isinstance(operator, (ast.In, ast.Eq)):
                    return True
    return False


REFERENCE_ONLY_LITERALS = ("reference_query_only", "reference_only")


def _test_uses_reference_only_literal(test) -> bool:
    return any(
        isinstance(node, ast.Constant) and node.value in REFERENCE_ONLY_LITERALS
        for node in ast.walk(test)
    )


def _test_uses_calibration_literal(test) -> bool:
    """A *pure* calibration-status literal in the test.

    ``calibrated_candidate_model`` is both a model layer and a functional-calibration
    status, so a test on it is not by itself evidence of a status-derived layer.
    """
    return any(
        isinstance(node, ast.Constant)
        and node.value in FUNCTIONAL_CALIBRATION_STATUSES
        and node.value not in MODEL_LAYERS
        for node in ast.walk(test)
    )


def _test_uses_layer_literal(test) -> bool:
    return any(
        isinstance(node, ast.Constant) and node.value in MODEL_LAYERS
        for node in ast.walk(test)
    )


def _test_uses_model_status(test) -> bool:
    if any(
        isinstance(node, ast.Constant) and node.value == "model_status"
        for node in ast.walk(test)
    ):
        return True
    return any(
        (isinstance(node, ast.Attribute) and node.attr == "model_status")
        for node in ast.walk(test)
    )


def _blocks_emit_layer(bodies) -> bool:
    """Whether a guarded block actually produces a model layer (or a layer field)."""
    for body in bodies:
        for statement in body or []:
            for node in ast.walk(statement):
                if isinstance(node, ast.Constant) and node.value in MODEL_LAYERS:
                    return True
                if isinstance(node, ast.Dict):
                    if "model_layer" in _dict_string_keys(node):
                        return True
                targets, value = _assign_pair(node)
                if any(LAYER_FIELD_RE.search(name) for name in targets):
                    return True
    return False


def _score_scoring_call(nodes) -> bool:
    """Whether a block reaches an HMM scoring action (call, binary name or tblout path)."""
    for statement in nodes:
        for node in ast.walk(statement):
            if isinstance(node, ast.Call) and SCORING_CALL_RE.search(_call_name(node.func) or ""):
                return True
            if isinstance(node, ast.Name) and SCORING_CALL_RE.search(node.id or ""):
                return True
            if isinstance(node, ast.Attribute) and SCORING_CALL_RE.search(node.attr or ""):
                return True
            if isinstance(node, ast.Constant) and isinstance(node.value, str) \
                    and SCORING_CALL_RE.search(node.value):
                return True
    return False


def _enclosing(function_stack, node):
    for func in reversed(function_stack):
        if func.lineno <= node.lineno <= getattr(func, "end_lineno", node.lineno):
            return func
    return None


class _AstAudit(ast.NodeVisitor):
    def __init__(self, tree, path):
        self.tree = tree
        self.path = path
        self.env, self.function_envs, self.parents = _build_scoped_taint_env(tree)
        self.docstring_lines = _docstring_lines(tree)
        self.findings = []
        self.functions = []
        self.decision_lines = set()
        self.covered_compares = set()

    def env_at(self, node) -> dict:
        owner = _enclosing_function(self.parents, node)
        if owner is not None and owner in self.function_envs:
            return self.function_envs[owner]
        return self.env

    # -- emission ---------------------------------------------------------
    def emit(self, line, signal, use_class, violation, rule, note, kind="ast"):
        self.findings.append(
            _finding(self.path, line, kind, signal, use_class, violation, rule, note)
        )
        if kind == "ast" and violation != "false":
            self.decision_lines.add(int(line))

    # -- visitors ---------------------------------------------------------
    def visit_FunctionDef(self, node):
        self.functions.append(node)
        self.generic_visit(node)
        self.functions.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_If(self, node):
        self._handle_decision(node.test, [node.body], [node.orelse], node)
        self.generic_visit(node)

    def visit_While(self, node):
        self._handle_decision(node.test, [node.body], [node.orelse], node)
        self.generic_visit(node)

    def visit_Compare(self, node):
        """Classify a reference_query_only comparison that is not a decision test."""
        if id(node) in self.covered_compares:
            return
        if not _test_uses_reference_only_literal(node):
            return
        enclosing = _enclosing(self.functions, node)
        body = [enclosing] if enclosing is not None else []
        self._emit_reference_only(node, node, [body] if body else [[]], [[]])

    def visit_ListComp(self, node):
        self._handle_comprehension(node)
        self.generic_visit(node)

    visit_SetComp = visit_ListComp
    visit_GeneratorExp = visit_ListComp
    visit_DictComp = visit_ListComp

    def visit_Call(self, node):
        name = _call_name(node.func)
        if name in {"sorted", "max", "min"} or (name == "sort"):
            for keyword in node.keywords:
                if keyword.arg in {"key", "reverse"} and _kinds_of(keyword.value, self.env_at(node)):
                    kinds = _kinds_of(keyword.value, self.env_at(node))
                    self.emit(
                        node.lineno, "discovery_score_sort_or_rank", "recall_only", "false",
                        "R1_recall_only_only",
                        f"discovery-layer value ({'/'.join(sorted(kinds))}) is only used to "
                        f"order candidates for recall/review; ranking is permitted under R1",
                    )
        self.generic_visit(node)

    def visit_Dict(self, node):
        env = self.env_at(node)
        kind_by_key = {
            key.value: value for key, value in zip(node.keys, node.values)
            if isinstance(key, ast.Constant) and isinstance(key.value, str)
        }
        # A dict whose values are all plain string constants is a naming/rename map
        # (``{"superfamily_claim": "discovery_superfamily_claim"}``), not a data
        # emission: no row field is filled from a value there.
        if kind_by_key and all(
            isinstance(value, ast.Constant) and isinstance(value.value, str)
            for value in kind_by_key.values()
        ):
            self.generic_visit(node)
            return
        if "model_layer" in kind_by_key:
            value = kind_by_key["model_layer"]
            value_kinds = _kinds_of(value, env)
            literals = _string_constants(value)
            if "layer" in value_kinds or any(
                DISCOVERY_LAYER_TOKEN_RE.search(literal) for literal in literals
            ):
                self.emit(
                    value.lineno, "discovery_layer_label_emitted", "recall_only", "false",
                    "R2_no_family_call",
                    "the emitted row carries an explicit discovery-layer model_layer label "
                    "(mandated provenance label, not a family call)",
                )
        for key, value in kind_by_key.items():
            value_kinds = _kinds_of(value, env)
            data_tainted = bool(value_kinds & {"score", "family"})
            if _discovery_name(key) and FAMILY_HINT_RE.search(key):
                # The field is *named* as discovery-derived: that is a labelled provenance
                # field, not a family call.  Whether a labelled discovery-layer hit field
                # still counts as a "family call" under R2 is a human decision.
                self.emit(
                    value.lineno, "discovery_family_field_emitted", "family_call",
                    "unreviewed", "R2_no_family_call",
                    f"row field {key!r} names a discovery-derived family/confidence output; "
                    "the field name labels its provenance, so whether this is a 'family call' "
                    "under R2 or a permitted labelled recall field needs a human decision",
                )
            elif FAMILY_HINT_RE.search(key) and data_tainted:
                labelled = _has_discovery_label([node], env)
                self.emit(
                    value.lineno, "discovery_family_field_emitted", "family_call",
                    "unreviewed" if labelled else "true", "R2_no_family_call",
                    (
                        f"row field {key!r} is filled from a discovery-layer value inside a "
                        "row that also carries an explicit discovery-layer label; whether a "
                        "labelled discovery-layer assignment is still a 'family call' under "
                        "R2 needs a human decision"
                    ) if labelled else (
                        f"row field {key!r} (a family/superfamily/subtype/confidence call "
                        "field) is filled from a discovery-layer value; R2 forbids the "
                        "discovery layer from producing a family call"
                    ),
                )
        self.generic_visit(node)

    def visit_Assign(self, node):
        self._handle_assignment(node.targets, node.value, node)
        self.generic_visit(node)

    def visit_AnnAssign(self, node):
        self._handle_assignment([node.target], node.value, node)
        self.generic_visit(node)

    # -- helpers ----------------------------------------------------------
    def _context_names(self, node, extra=()) -> set:
        """Names that describe *what* a decision builds: function, assignment targets."""
        names = {name for name in extra if name}
        owner = _enclosing_function(self.parents, node)
        if owner is not None:
            names.add(getattr(owner, "name", "") or "")
        statement = node
        while statement is not None and not isinstance(statement, ast.stmt):
            statement = self.parents.get(statement)
        if statement is not None:
            targets, _value = _assign_pair(statement)
            names.update(targets)
            if isinstance(statement, ast.Return):
                owner_name = getattr(owner, "name", "") if owner is not None else ""
                names.add(owner_name)
        return {name for name in names if name}

    def _context_effects(self, node, test, extra_names=()) -> set:
        """Extra effect flags that need names, not just the block body."""
        effects = set()
        haystack = " ".join(self._context_names(node, extra_names))
        if CANDIDATE_COUNT_RE.search(haystack):
            effects.add("candidate_set")
        if TRAINING_HINT_RE.search(haystack) or TRAINING_HINT_RE.search(
            " ".join(_string_constants(test)) + " " + _source_names(test)
        ):
            effects.add("training_selection")
        if self._function_has_scoring(node):
            effects.add("scoring_context")
        return effects

    def _handle_decision(self, test, bodies, other_bodies, node):
        env = self.env_at(node)
        kinds = _kinds_of(test, env)
        self.covered_compares.add(id(test))
        for child in ast.walk(test):
            if isinstance(child, ast.Compare):
                self.covered_compares.add(id(child))

        body_effects = _effects(bodies[0]) if bodies else set()
        else_effects = _effects(other_bodies[0]) if other_bodies else set()
        effects = body_effects | else_effects
        effects |= self._context_effects(node, test)
        if _test_compares_threshold(test, env):
            effects.add("score_threshold")

        if _test_uses_reference_only_literal(test):
            self._emit_reference_only(node, test, bodies, other_bodies)
            return

        if _test_uses_model_status(test):
            self._emit_model_status(node, test, effects)
            return

        if (_test_uses_layer_literal(test) and (kinds & {"status"})) or (
            _test_uses_calibration_literal(test) and (kinds & {"field", "layer"})
        ):
            if _blocks_emit_layer(bodies):
                self._emit_layer_status(test, bodies, node, direction="reverse")
                return

        if not kinds:
            return

        labelled = _has_discovery_label(bodies[0], env) if bodies else False
        use_class, violation, rule, note = _decision_result(
            kinds, effects, labelled, f"lineno {node.lineno}",
        )
        signal = "discovery_decision"
        if use_class == "thresholding":
            signal = "discovery_score_filter"
        elif use_class == "family_call":
            signal = "discovery_family_assignment"
        elif use_class == "counting":
            signal = "discovery_count"
        elif use_class == "recall_only":
            signal = "discovery_layer_gate_compliant"
        self.emit(node.lineno, signal, use_class, violation, rule, note)

        if _test_uses_layer_literal(test) and (kinds & {"field"}) and not (kinds & {"score"}):
            self._emit_layer_status(test, bodies, node, direction="forward")

    def _handle_comprehension(self, node):
        target_names = set()
        for generator in node.generators:
            target_names.update(_target_names(generator.target))
        for generator in node.generators:
            for condition in generator.ifs:
                self.covered_compares.add(id(condition))
                for child in ast.walk(condition):
                    if isinstance(child, ast.Compare):
                        self.covered_compares.add(id(child))
                kinds = _kinds_of(condition, self.env_at(node))
                element = getattr(node, "key", None) or node.elt
                if isinstance(element, ast.Name) and element.id in target_names:
                    effects = {"select"}
                else:
                    effects = _effects([ast.Expr(value=element)]) if not isinstance(element, ast.Name) else set()
                effects |= self._context_effects(node, condition)
                if _test_compares_threshold(condition, self.env_at(node)):
                    effects.add("score_threshold")
                if _test_uses_reference_only_literal(condition):
                    self._emit_reference_only(node, condition, [[ast.Pass()]], [[]])
                    continue
                if _test_uses_model_status(condition):
                    self._emit_model_status(node, condition, effects | {"select"})
                    continue
                if not kinds:
                    continue
                use_class, violation, rule, note = _decision_result(
                    kinds, effects, False, f"comprehension at lineno {node.lineno}",
                )
                signal = "discovery_score_filter" if use_class == "thresholding" else (
                    "discovery_family_assignment" if use_class == "family_call" else (
                        "discovery_count" if use_class == "counting" else (
                            "discovery_layer_gate_compliant" if use_class == "recall_only"
                            else "discovery_decision"
                        )
                    )
                )
                self.emit(node.lineno, signal, use_class, violation, rule, note)

    def _handle_assignment(self, targets, value, node):
        if value is None:
            return
        names = []
        for target in targets:
            names.extend(_target_names(target))
        env = self.env_at(node)
        value_kinds = _kinds_of(value, env)
        # A container literal is not a scalar discovery value: what it *emits* is
        # classified by visit_Dict (field-level), so a dict/list literal assigned to a
        # family-named variable must not be read as a family call.
        container_literal = isinstance(value, (ast.Dict, ast.List, ast.Set, ast.Tuple))
        discovery_data = bool(value_kinds & {"score", "family"})
        if not discovery_data:
            discovery_data = any(
                DISCOVERY_LAYER_TOKEN_RE.search(constant) for constant in _string_constants(value)
            )

        # R7: functional_calibration_status derived from a layer test.
        status_targets = [name for name in names if CALIBRATION_HINT_RE.search(name)]
        if status_targets and any(
            isinstance(constant, str) and constant in FUNCTIONAL_CALIBRATION_STATUSES
            for constant in _string_constants(value)
        ):
            enclosing = _enclosing(self.functions, node)
            if enclosing is not None and self._function_tests_layer(enclosing):
                self.emit(
                    node.lineno, "calibration_status_derived_from_layer", "family_call",
                    "true", "R7_layer_status_independent",
                    f"{status_targets[0]} is assigned inside a test on the model layer: R7 "
                    "forbids deriving functional_calibration_status from model_layer",
                )

        # R7 reverse: model_layer derived from a functional calibration status.
        if any(LAYER_FIELD_RE.search(name) or name == "model_layer" for name in names):
            if any(constant in MODEL_LAYERS for constant in _string_constants(value)):
                enclosing = _enclosing(self.functions, node)
                if enclosing is not None and self._function_tests_calibration(enclosing):
                    self.emit(
                        node.lineno, "model_layer_derived_from_calibration_status",
                        "family_call", "unreviewed", "R7_layer_status_independent",
                        "model_layer is returned inside a test on "
                        "functional_calibration_status: R7 forbids deriving either field "
                        "from the other; a human must decide whether resolving a stated "
                        "calibrated status into its layer is a sanctioned normalisation",
                    )

        if not discovery_data or container_literal:
            return
        # A module-level SCREAMING_CASE constant that merely *names* a rule, reason or
        # column (REASON_DISCOVERY_FAMILY_CALL = "discovery_family_call_present") is not a
        # family call.
        constant_label = (
            len(names) == 1
            and CONSTANT_TARGET_RE.match(names[0]) is not None
            and all(isinstance(constant, str) for constant in _string_constants(value))
            and bool(_string_constants(value))
        )
        for name in names:
            if constant_label:
                break
            if STRONG_CALL_TARGET_RE.search(name):
                labelled = _has_discovery_label([node], env)
                self.emit(
                    node.lineno, "discovery_family_assignment", "family_call",
                    "unreviewed" if labelled else "true", "R2_no_family_call",
                    (
                        f"{name!r} is bound from a discovery-layer value inside a row that "
                        "carries an explicit discovery-layer label; whether that is still a "
                        "'family call' under R2 needs a human decision"
                    ) if labelled else (
                        f"{name!r} is bound from a discovery-layer value "
                        f"({'/'.join(sorted(value_kinds))}); R2 forbids the discovery layer "
                        "from producing a family call"
                    ),
                )
            elif FAMILY_TARGET_RE.search(name):
                self.emit(
                    node.lineno, "discovery_family_assignment", "family_call", "unreviewed",
                    "R2_no_family_call",
                    f"{name!r} is a family/superfamily/confidence-shaped local bound from a "
                    f"discovery-layer value ({'/'.join(sorted(value_kinds))}); whether that "
                    "value becomes a family call is not visible here, so a human must read "
                    "the data flow",
                )
            if name == "model_layer" and "layer" in value_kinds:
                self.emit(
                    node.lineno, "discovery_layer_label_emitted", "recall_only", "false",
                    "R2_no_family_call",
                    "model_layer is set to the discovery-layer label (mandated provenance)",
                )
            if COUNT_TARGET_RE.search(name) and discovery_data:
                if CANDIDATE_COUNT_RE.search(name):
                    self.emit(
                        node.lineno, "discovery_count", "counting", "true",
                        "R1_recall_only_only",
                        f"{name!r} is a candidate/high-confidence count fed by a "
                        "discovery-layer value; R1 forbids the discovery layer from deciding "
                        "the candidate set",
                    )
                else:
                    self.emit(
                        node.lineno, "discovery_count", "counting", "false",
                        "R1_recall_only_only",
                        f"{name!r} counts discovery-layer rows; it is not named as a "
                        "candidate/high-confidence count",
                    )
            if SCORE_TARGET_RE.search(name) and discovery_data:
                self.emit(
                    node.lineno, "discovery_score_recorded", "recall_only", "false",
                    "R1_recall_only_only",
                    f"{name!r} records a discovery-layer score for recall/review only",
                )

    def _function_tests_layer(self, func) -> bool:
        env = self.function_envs.get(func, self.env)
        for node in ast.walk(func):
            if isinstance(node, (ast.If, ast.While)) and _test_uses_layer_literal(node.test):
                return True
            if isinstance(node, ast.If) and _kinds_of(node.test, env) & {"field", "layer"}:
                return True
        return False

    def _function_has_scoring(self, node) -> bool:
        """Whether a scoring call (hmmsearch/tblout/pyhmmer) is reachable in this function.

        A conservative ``True`` at module level: without an enclosing function the audit
        cannot bound the reachable code.
        """
        owner = _enclosing_function(self.parents, node)
        if owner is None:
            return True
        return _score_scoring_call([owner])

    def _function_tests_calibration(self, func) -> bool:
        env = self.function_envs.get(func, self.env)
        for node in ast.walk(func):
            if isinstance(node, ast.If):
                if _test_uses_calibration_literal(node.test):
                    return True
                if _kinds_of(node.test, env) & {"status"}:
                    return True
        return False

    def _emit_reference_only(self, node, test, bodies, other_bodies):
        body_effects = _effects(bodies[0]) if bodies else set()
        other_effects = _effects(other_bodies[0]) if other_bodies else set()
        scoring = _score_scoring_call((bodies[0] if bodies else []) + (other_bodies[0] if other_bodies else []))
        exclusion = _comparison_is_exclusion(test)
        if exclusion:
            self.emit(
                node.lineno, "reference_query_only_model_layer", "recall_only", "false",
                "R6_reference_query_only_unscoreable",
                "reference_query_only is excluded from the discriminating set: it has no "
                "HMM at all, so excluding it is the R6-compliant behaviour",
            )
            return
        if scoring or {"score"} & (body_effects | other_effects):
            self.emit(
                node.lineno, "reference_query_only_model_layer", "thresholding", "true",
                "R6_reference_query_only_unscoreable",
                "a reference_query_only profile is selected for scoring; it has no HMM, so "
                "scoring it fabricates a score that does not exist (R6)",
            )
            return
        if _comparison_is_inclusion(test):
            labeling_effects = {"collect", "label", "suppress", "weaken", "raise", "count", "return"}
            if body_effects and body_effects <= labeling_effects and not self._function_has_scoring(node):
                self.emit(
                    node.lineno, "reference_query_only_model_layer", "recall_only", "false",
                    "R6_reference_query_only_unscoreable",
                    "the block records/refuses a reference_query_only profile (label, reason "
                    "or count) and no scoring call (hmmsearch/tblout/pyhmmer) is reachable in "
                    "this function: R6 is not engaged",
                )
                return
            self.emit(
                node.lineno, "reference_query_only_model_layer", "unknown", "unreviewed",
                "R6_reference_query_only_unscoreable",
                "reference_query_only is selected by an inclusion test whose effect this "
                "audit cannot classify; no scoring call is visible in this function, but "
                "whether a caller scores the selected profiles (fabricating a score for a "
                "profile that has no HMM) needs a human to follow the value",
            )

    def _emit_model_status(self, node, test, effects):
        if not _comparison_is_inclusion(test):
            return
        note_tail = (
            "the legacy model_status column reads 'trained' for the discovery layer too "
            "(pipeline/scripts/build_phaded_profiles.py:407-408 writes model_status='trained' "
            "together with model_layer='discovery_hmm_uncalibrated'), so a decision keyed on "
            "model_status cannot separate the discovery layer from a discriminating model; "
            "the criterion must read model_layer"
        )
        if effects & {"assign"}:
            self.emit(node.lineno, "model_status_used_as_layer_criterion", "family_call",
                      "true", "R2_no_family_call",
                      "a family/subtype/confidence value is produced from a model_status "
                      "criterion; " + note_tail)
        elif effects & {"candidate_set", "candidate_count", "scoring_context"}:
            self.emit(node.lineno, "model_status_used_as_layer_criterion", "thresholding",
                      "true", "R1_recall_only_only",
                      "profiles are selected for HMM scoring (or into a candidate/"
                      "high-confidence set) from a model_status criterion; " + note_tail)
        elif effects & {"drop", "select", "score"}:
            self.emit(node.lineno, "model_status_used_as_layer_criterion", "thresholding",
                      "unreviewed", "R1_recall_only_only",
                      "profiles/rows are partitioned on a model_status criterion, but the "
                      "collection is not named as a candidate/high-confidence set (report "
                      "partition or recall table); " + note_tail)
        elif "count" in effects:
            self.emit(node.lineno, "model_status_used_as_layer_criterion", "counting",
                      "false", "R1_recall_only_only",
                      "the legacy model_status column is only counted/reported here; " + note_tail)
        else:
            self.emit(node.lineno, "model_status_used_as_layer_criterion", "unknown",
                      "unreviewed", "R1_recall_only_only",
                      "a model_status criterion drives a block whose effect this audit "
                      "cannot classify; " + note_tail)

    def _emit_layer_status(self, test, bodies, node, direction):
        literals = [
            constant for constant in _string_constants(test)
            if constant in MODEL_LAYERS or constant in FUNCTIONAL_CALIBRATION_STATUSES
        ]
        if direction == "reverse":
            self.emit(
                node.lineno, "model_layer_derived_from_calibration_status", "family_call",
                "unreviewed", "R7_layer_status_independent",
                "a model-layer literal is returned inside a test on "
                f"{'/'.join(sorted(set(literals)))}: R7 forbids deriving either field from "
                "the other; a human must decide whether resolving a stated calibrated "
                "status into its layer is a sanctioned normalisation",
            )
        else:
            self.emit(
                node.lineno, "layer_gate_on_model_layer", "recall_only", "false",
                "R1_recall_only_only",
                "a test on the model_layer value routes the block; it does not score or "
                "filter candidates on a discovery-layer score",
            )


# ---------------------------------------------------------------------------
# Text scan (config/path literals and *.sh)
# ---------------------------------------------------------------------------


def scan_text_signals(source_text: str, path="<memory>") -> list:
    """Raw text matches for governance-relevant literals (no classification).

    This is the documented text heuristic.  It deliberately over-matches: a *comment*
    or docstring that merely mentions a discovery-score comparison produces a row here,
    which is exactly why violations are decided by the AST in :func:`classify_usage`.
    """
    script = _display_path(path)
    lines = source_text.splitlines()
    rows = []
    for index, line in enumerate(lines, start=1):
        if REGISTRY_TOKEN_RE.search(line):
            rows.append(_finding(
                script, index, "text", "formal_scan_models_reference", "registry_write",
                "unreviewed", "R3_registry_exclusion",
                "text match on the frozen formal-scan registry path",
            ))
        if DEFERRED_ARCHIVE_TOKEN_RE.search(line):
            rows.append(_finding(
                script, index, "text", "with_lipase_deferred_archive_reference",
                "deferred_layer_touch", "unreviewed", "R4_deferred_preserved",
                "text match on the with-lipase (DED_hfam_2) deferred archive",
            ))
        if DISCOVERY_NAME_RE.search(line) and COMPARISON_RE.search(line) \
                and THRESHOLD_RE.search(line):
            rows.append(_finding(
                script, index, "text", "discovery_score_comparison", "unknown",
                "unreviewed", "R1_recall_only_only",
                "text match on a discovery-layer comparison (may be prose)",
            ))
    return rows


def _line_is_prose(source_text: str, tree, line: int) -> bool:
    lines = source_text.splitlines()
    if line - 1 >= len(lines):
        return False
    text = lines[line - 1]
    stripped = text.strip()
    if stripped.startswith("#"):
        return True
    if tree is None:
        return False
    return line in _docstring_lines(tree)


def _token_bound_names(tree, token_re) -> set:
    """Names bound to an expression that contains the literal governance path."""
    names = set()
    if tree is None:
        return names
    for node in ast.walk(tree):
        targets, value = _assign_pair(node)
        if not targets or value is None:
            continue
        if any(token_re.search(item) for item in _string_constants(value)):
            names.update(targets)
    return names


def _targets_at_line(tree, line: int) -> list:
    """Assignment targets of the statement containing *line* (no token check)."""
    if tree is None:
        return []
    best = None
    for node in ast.walk(tree):
        if not isinstance(node, ast.stmt):
            continue
        start = getattr(node, "lineno", None)
        end = getattr(node, "end_lineno", start)
        if start is None or not (start <= line <= end):
            continue
        span = end - start
        if best is None or span < best[0]:
            best = (span, node)
    if best is None:
        return []
    targets, _value = _assign_pair(best[1])
    return targets


def _names_bound_at_line(tree, token_re, line: int) -> set:
    """Names bound on *this* line to an expression holding the governance path.

    Scoping matters: one writing name anywhere in a file must not turn every comment
    that mentions the deferred archive into a write.
    """
    if tree is None:
        return set()
    best = None
    for node in ast.walk(tree):
        if not isinstance(node, ast.stmt):
            continue
        start = getattr(node, "lineno", None)
        end = getattr(node, "end_lineno", start)
        if start is None or not (start <= line <= end):
            continue
        span = end - start
        if best is None or span < best[0]:
            best = (span, node)
    if best is None:
        return set()
    targets, value = _assign_pair(best[1])
    if value is None or not targets:
        return set()
    if any(token_re.search(item) for item in _string_constants(value)):
        return set(targets)
    return set()


def _name_write_mode(name: str, source_text: str) -> str:
    """``write``/``delete``/``read`` for a name bound to a governance path."""
    escaped = re.escape(name)
    patterns = (
        (rf"\b{escaped}\s*\.\s*unlink\(|shutil\.rmtree\(\s*{escaped}\b|"
         rf"os\.remove\(\s*{escaped}\b|os\.unlink\(\s*{escaped}\b", "delete"),
        (rf"\b{escaped}\s*\.\s*open\([^)]*[\"'](?:w|a|x)[\"']", "write"),
        (rf"\b{escaped}\s*\.\s*(?:write_text|write_bytes|touch|mkdir)\(", "write"),
        (rf"\b{escaped}\s*\.\s*(?:replace|rename)\(", "write"),
    )
    for pattern, mode in patterns:
        if re.search(pattern, source_text):
            return mode
    return "read"


def _classify_path_reference(source_text, tree, row, token_re) -> tuple:
    """Decide read-only vs write for one registry/deferred reference row.

    Two mechanisms, in order: (1) Python —the name the governance path is bound to is
    opened/written/deleted somewhere in the file; (2) shell/text —the path literal is the
    destination of a redirection or of ``cp``/``mv`` on the same line.  A governance path
    that is only hashed, read or copied *into a run* is read-only.
    """
    line = row["line"]
    lines = source_text.splitlines()
    text = lines[line - 1] if line - 1 < len(lines) else ""

    if tree is not None:
        bound = _names_bound_at_line(tree, token_re, line)
        modes = {_name_write_mode(name, source_text) for name in bound}
        if "delete" in modes:
            return "delete", False
        if "write" in modes:
            return "write", False

    match = token_re.search(text)
    if match is not None:
        before = text[: match.start()]
        after = text[match.end():]
        if re.search(
            r"(?:open\([^)]*[\"']?(?:w|a|x)|write_text\(|write_bytes\(|\.unlink\(|"
            r"shutil\.rmtree\(|os\.remove\(|>>?\s*[\"']?)\s*$",
            before,
        ):
            return "write", False
        if re.search(r"(^|[;&|]\s*)(cp|mv|install|tee)\s", before):
            # `cp "$REGISTRY" "$RUN_ROOT/inputs/formal_scan_models.tsv"`
            destination_is_config = "config" in before or "config" in after
            return ("write" if destination_is_config else "copy_into_run"), not destination_is_config
        if SHELL_DELETE_RE.search(text):
            return "delete", False
    elif SHELL_DELETE_RE.search(text):
        return "delete", False
    return "read", True


DELETE_MARKER_RE = re.compile(r"\.unlink\(|shutil\.rmtree\(|os\.remove\(|\.rmdir\(|os\.unlink\(")


def classify_text_usage(source_text: str, path="<memory>", tree=None) -> list:
    """Classify the documented text matches into audit rows."""
    script = _display_path(path)
    findings = []
    raw_rows = scan_text_signals(source_text, path)
    ast_decision_lines = set()
    if tree is not None:
        ast_decision_lines = {
            node.lineno for node in ast.walk(tree)
            if isinstance(node, (ast.If, ast.ListComp, ast.SetComp, ast.GeneratorExp,
                                 ast.DictComp))
        }
    for row in raw_rows:
        if row["signal"] == "formal_scan_models_reference":
            mode, read_only = _classify_path_reference(source_text, tree, row, REGISTRY_TOKEN_RE)
            if mode == "delete":
                findings.append(_finding(
                    script, row["line"], "text", "formal_scan_models_reference",
                    "registry_write", "true", "R3_registry_exclusion",
                    "the frozen formal-scan registry is deleted/renamed by this file; the "
                    "registry is a frozen governance artefact (read_only=false)",
                ))
            elif mode == "write":
                findings.append(_finding(
                    script, row["line"], "text", "formal_scan_models_reference",
                    "registry_write", "true", "R3_registry_exclusion",
                    "this file opens the formal-scan registry for writing; the discovery "
                    "layer must never be written into pipeline/config/formal_scan_models.tsv "
                    "(read_only=false)",
                ))
            elif mode == "copy_into_run":
                findings.append(_finding(
                    script, row["line"], "text", "formal_scan_models_reference",
                    "registry_write", "false", "R3_registry_exclusion",
                    "copies the frozen registry into a run's inputs/ for provenance "
                    "(read_only=true: the config itself is only read)",
                ))
            else:
                findings.append(_finding(
                    script, row["line"], "text", "formal_scan_models_reference",
                    "registry_write", "false", "R3_registry_exclusion",
                    "read-only reference to the frozen formal-scan registry; the audit's "
                    "own registry-absence check depends on this file staying frozen "
                    "(read_only=true)",
                ))
            continue
        if row["signal"] == "with_lipase_deferred_archive_reference":
            mode, read_only = _classify_path_reference(
                source_text, tree, row, DEFERRED_ARCHIVE_TOKEN_RE
            )
            source_lines = source_text.splitlines()
            line_text = source_lines[row["line"] - 1] if row["line"] - 1 < len(source_lines) else ""
            count_names = sorted(
                name for name in _targets_at_line(tree, row["line"])
                if COUNT_TARGET_RE.search(name)
            )
            if not count_names and _targets_at_line(tree, row["line"]) and re.search(
                r"(?<![\w.])[0-9]{4,}(?![\w.])", line_text
            ):
                # `REVIEW_WITH_LIPASE_DEFERRED = 1206655`: the layer's own size recorded as
                # a review constant.  R5 is about where that number is *used*; the audit can
                # only confirm that the constant itself is not a candidate count.
                count_names = sorted(_targets_at_line(tree, row["line"]))
            if count_names and mode != "delete":
                findings.append(_finding(
                    script, row["line"], "text", "with_lipase_deferred_archive_reference",
                    "counting", "false", "R5_deferred_out_of_counts",
                    f"{count_names[0]!r} records the deferred layer's own size; R5 requires "
                    "this number to stay outside every candidate/high-confidence count "
                    "(read_only=true). Whether the constant is consumed by a candidate count "
                    "elsewhere is not visible here",
                ))
                continue
            if mode == "delete":
                findings.append(_finding(
                    script, row["line"], "text", "with_lipase_deferred_archive_reference",
                    "deferred_layer_touch", "true", "R4_deferred_preserved",
                    "the with-lipase deferred archive is deleted/renamed; AGENTS.md requires "
                    "the deferred structural-validation layer to be preserved "
                    "(read_only=false)",
                ))
            elif mode == "write":
                findings.append(_finding(
                    script, row["line"], "text", "with_lipase_deferred_archive_reference",
                    "deferred_layer_touch", "false", "R4_deferred_preserved",
                    "writes a with-lipase deferred artefact; rebuilding/preserving the "
                    "deferred layer is required by AGENTS.md, it is not a candidate output "
                    "(read_only=false)",
                ))
            else:
                findings.append(_finding(
                    script, row["line"], "text", "with_lipase_deferred_archive_reference",
                    "deferred_layer_touch", "false", "R4_deferred_preserved",
                    "read-only touch of the with-lipase deferred archive; the layer is "
                    "preserved and must stay outside candidate counts (read_only=true)",
                ))
            continue
        # discovery_score_comparison
        if _line_is_prose(source_text, tree, row["line"]):
            findings.append(_finding(
                script, row["line"], "text", "discovery_score_comparison_in_prose",
                "unknown", "unreviewed", "R1_recall_only_only",
                "text-only match inside a comment/docstring; not an executable decision, "
                "recorded so the audit does not silently pass over it",
            ))
        elif tree is None:
            findings.append(_finding(
                script, row["line"], "text", "discovery_score_comparison", "unknown",
                "unreviewed", "R1_recall_only_only",
                "shell line comparing a discovery-layer value; *.sh has no AST in this "
                "audit, so the effect is unreviewed",
            ))
        elif row["line"] not in ast_decision_lines:
            findings.append(_finding(
                script, row["line"], "text", "discovery_score_comparison_in_prose",
                "unknown", "unreviewed", "R1_recall_only_only",
                "text-only match that is not an AST decision node; a human must read it",
            ))
        # else: the AST already classified the same line
    return findings


# ---------------------------------------------------------------------------
# Prose claims
# ---------------------------------------------------------------------------

def _prose_lines(source_text: str, tree) -> list:
    lines = source_text.splitlines()
    out = []
    doc_lines = _docstring_lines(tree) if tree is not None else set()
    for index, text in enumerate(lines, start=1):
        stripped = text.strip()
        if stripped.startswith("#") and stripped:
            out.append((index, stripped))
        elif index in doc_lines and stripped:
            out.append((index, stripped))
    return out


def _prose_findings(source_text: str, tree, path) -> list:
    script = _display_path(path)
    findings = []
    emitted = set()
    for line, text in _prose_lines(source_text, tree):
        if not DEFERRED_ARCHIVE_TOKEN_RE.search(text) and not _discovery_name(text):
            continue
        rule = None
        note = ""
        if (DEFERRED_ARCHIVE_TOKEN_RE.search(text) or _discovery_name(text)) \
                and PROSE_VERBS_RE.search(text):
            rule = "R1_recall_only_only"
            note = ("prose-only governance claim (docstring/comment) about filtering, "
                    "deleting, demoting or counting; static AST analysis cannot confirm or "
                    "refute behaviour that is only asserted in prose")
        elif _discovery_name(text) and PROSE_FAMILY_RE.search(text):
            rule = "R2_no_family_call"
            note = ("prose-only claim (docstring/comment) about a discovery-layer family/"
                    "superfamily assignment; the AST cannot see the discovery provenance of "
                    "a value that arrives through data, so a human must read the block")
        if rule is None:
            continue
        key = ("prose", rule)
        if key in emitted:
            continue
        emitted.add(key)
        findings.append(_finding(
            script, line, "ast" if line in (_docstring_lines(tree) if tree else set()) else "text",
            "prose_only_governance_claim", "unknown", "unreviewed", rule, note,
        ))
    return findings


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

def classify_usage(source_text: str, path="<memory>") -> list:
    """Classify one Python source text; pure function, no I/O, no execution."""
    script = _display_path(path)
    findings = []
    tree = None
    try:
        tree = ast.parse(source_text)
    except SyntaxError as error:
        findings.append(_finding(
            script, getattr(error, "lineno", 1) or 1, "ast", "unparseable_source",
            "unknown", "unreviewed", "NONE",
            f"source could not be parsed ({error.msg}); nothing about this file's code "
            "paths could be reviewed (its text references are still scanned)",
        ))
    else:
        visitor = _AstAudit(tree, script)
        visitor.visit(tree)
        findings.extend(visitor.findings)
    findings.extend(classify_text_usage(source_text, path, tree))
    findings.extend(_prose_findings(source_text, tree, path))
    return _sort_findings(_dedupe(findings))


def classify_shell_usage(source_text: str, path="<memory>") -> list:
    """Classify one shell script through the documented text heuristic (no AST)."""
    return _sort_findings(_dedupe(classify_text_usage(source_text, path, None)))


def classify_file(path) -> list:
    target = Path(path)
    source_text = target.read_text(encoding="utf-8", errors="replace")
    if target.suffix == ".py":
        return classify_usage(source_text, target)
    if target.suffix == ".sh":
        return classify_shell_usage(source_text, target)
    raise AuditError(f"unsupported file type for the audit: {target}")


def classify_repository(scripts_dir, exclude=()) -> list:
    """Scan every ``*.py`` (AST) and ``*.sh`` (text) directly under *scripts_dir*."""
    directory = Path(scripts_dir)
    if not directory.is_dir():
        raise AuditError(f"--scripts-dir is not a directory: {directory}")
    excluded = {str(name) for name in exclude}
    findings = []
    for pattern in ("*.py", "*.sh"):
        for path in sorted(directory.glob(pattern)):
            if path.name in excluded:
                continue
            findings.extend(classify_file(path))
    return _sort_findings(_dedupe(findings))


def check_registry_absence(registry_path) -> dict:
    """Read the frozen registry and report whether the discovery layer is listed.

    The check matches the discovery-layer *identifiers*
    (``discovery_hmm_uncalibrated`` / ``cys_discovery_uncalibrated`` / ``uncalibrated``),
    never the bare word ``discovery``: ``ePhaZ_broad_discovery`` is an ePhaZ model that
    legitimately lives in this registry and has nothing to do with the PhaDED discovery
    layer.
    """
    path = Path(registry_path)
    if not path.is_file():
        raise AuditError(f"--registry is not a file: {path}")
    raw = path.read_bytes()
    text = raw.decode("utf-8-sig", errors="replace")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    columns = list(rows[0].keys()) if rows else []
    matched = []
    for row in rows:
        values = [(row.get(column) or "") for column in columns]
        if any(DISCOVERY_LAYER_TOKEN_RE.search(value) for value in values) or any(
            value.strip() == "uncalibrated" for value in values
        ):
            matched.append(
                row.get("model") or row.get(columns[0]) if columns else "<unknown>"
            )
    return {
        "path": _display_path(path),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "size_bytes": len(raw),
        "columns": columns,
        "row_count": len(rows),
        "models": [row.get("model", "") for row in rows],
        "discovery_layer_models": matched,
        "discovery_layer_absent": not matched,
        "note": (
            "the discovery layer (model_layer='discovery_hmm_uncalibrated') is absent from "
            "the frozen formal-scan registry"
            if not matched else
            "the discovery layer appears in the frozen formal-scan registry; R3 is broken"
        ),
    }


def registry_findings(report: dict) -> list:
    findings = []
    for index, model in enumerate(report.get("discovery_layer_models", []), start=2):
        findings.append(_finding(
            report["path"], index, "text", "formal_scan_models_reference", "registry_write",
            "true", "R3_registry_exclusion",
            f"the frozen formal-scan registry lists the discovery layer as {model!r}; "
            "AGENTS.md forbids the discovery layer from entering this registry",
        ))
    return findings


def summarize(findings, scripts_dir, scanned: dict, excluded, registry_report=None) -> dict:
    use_class_counts: dict = {}
    violation_counts: dict = {}
    for row in findings:
        use_class_counts[row["use_class"]] = use_class_counts.get(row["use_class"], 0) + 1
        violation_counts[row["violation"]] = violation_counts.get(row["violation"], 0) + 1
    files_with_violations = sorted(
        {row["script"] for row in findings if row["violation"] == "true"}
    )
    unreviewed = [row for row in findings if row["violation"] == "unreviewed"]
    return {
        "schema_version": 1,
        "status": "read_only_static_audit",
        "scripts_dir": _display_path(scripts_dir),
        "scripts_scanned": scanned.get("python", 0) + scanned.get("shell", 0),
        "scripts_scanned_python": scanned.get("python", 0),
        "scripts_scanned_shell": scanned.get("shell", 0),
        "excluded": sorted(str(name) for name in excluded),
        "findings": len(findings),
        "violations": violation_counts.get("true", 0),
        "unreviewed": violation_counts.get("unreviewed", 0),
        "compliant": violation_counts.get("false", 0),
        "use_class_counts": dict(sorted(use_class_counts.items())),
        "violation_counts": dict(sorted(violation_counts.items())),
        "files_with_violations": files_with_violations,
        "unreviewed_prose_claims": sum(
            1 for row in unreviewed if row["signal"] == "prose_only_governance_claim"
        ),
        "registry": registry_report,
        "limits": list(AUDIT_LIMITS),
        "rules": dict(RULES),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="audit_discovery_layer_usage.py",
        description=(
            "Read-only static audit of discovery-layer (discovery_hmm_uncalibrated) and "
            "with-lipase deferred-layer (DED_hfam_2) downstream usage in pipeline/scripts. "
            "Writes discovery_layer_usage_audit.tsv + discovery_layer_usage_audit_summary.json "
            "into a new/empty --out-dir; exits 1 when a violation is found, 2 on a usage or "
            "validation error. It never writes into runs/, results/ or deploy/ and never "
            "modifies a scanned script."
        ),
        epilog=(
            "use_class vocabulary: " + ", ".join(USE_CLASSES) + ". "
            "violation: true | false | unreviewed (a construct the audit could not decide is "
            "never silently false). Rules: " + ", ".join(RULES) + "."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--scripts-dir", type=Path, default=Path("pipeline/scripts"),
                        help="directory holding pipeline scripts (default: pipeline/scripts)")
    parser.add_argument("--registry", type=Path, default=None,
                        help="frozen formal-scan registry to check for the discovery layer "
                             "(read-only; e.g. pipeline/config/formal_scan_models.tsv)")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="output directory; must not exist or must be empty")
    parser.add_argument("--exclude", action="append", default=None,
                        help="script filename to exclude (repeatable). Default: this audit's "
                             "own file, because it contains the token tables it searches for")
    parser.add_argument("--print-table", action="store_true",
                        help="also print the TSV rows to stdout")
    return parser


def _validated_directory(path: Path, flag: str) -> Path:
    if not path.exists():
        raise AuditError(f"{flag} does not exist: {path}")
    if not path.is_dir():
        raise AuditError(f"{flag} is not a directory: {path}")
    return path


def _validated_file(path: Path, flag: str) -> Path:
    if not path.exists():
        raise AuditError(f"{flag} does not exist: {path}")
    if not path.is_file():
        raise AuditError(f"{flag} is not a file: {path}")
    return path


def _validated_out_dir(path: Path) -> Path:
    if path.exists():
        if not path.is_dir():
            raise AuditError(f"--out-dir is not a directory: {path}")
        if any(path.iterdir()):
            raise AuditError(
                f"--out-dir is a non-empty directory: {path}; refusing to write the audit "
                "into an existing non-empty output directory (frozen evidence must not be "
                "overwritten)"
            )
    return path


def write_outputs(findings, summary: dict, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    table_path = out_dir / "discovery_layer_usage_audit.tsv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(REQUIRED_COLUMNS), delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(findings)
    summary_path = out_dir / "discovery_layer_usage_audit_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return table_path, summary_path


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        scripts_dir = _validated_directory(args.scripts_dir, "--scripts-dir")
        registry = _validated_file(args.registry, "--registry") if args.registry else None
        out_dir = _validated_out_dir(args.out_dir)
    except AuditError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if args.exclude is None:
        excluded = {Path(__file__).name}
    else:
        excluded = {name for name in args.exclude if name}

    try:
        findings = classify_repository(scripts_dir, exclude=excluded)
    except AuditError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    scanned = {
        "python": sum(
            1 for path in scripts_dir.glob("*.py") if path.name not in excluded
        ),
        "shell": sum(
            1 for path in scripts_dir.glob("*.sh") if path.name not in excluded
        ),
    }

    registry_report = None
    if registry is not None:
        try:
            registry_report = check_registry_absence(registry)
        except AuditError as error:
            print(f"error: {error}", file=sys.stderr)
            return 2
        findings = _sort_findings(_dedupe(findings + registry_findings(registry_report)))

    summary = summarize(findings, scripts_dir, scanned, excluded, registry_report)

    try:
        table_path, summary_path = write_outputs(findings, summary, out_dir)
    except OSError as error:
        print(f"error: could not write the audit outputs: {error}", file=sys.stderr)
        return 2

    if args.print_table:
        print("\t".join(REQUIRED_COLUMNS))
        for row in findings:
            print("\t".join(str(row[column]) for column in REQUIRED_COLUMNS))
    print(json.dumps({
        "audit_table": _display_path(table_path),
        "summary": _display_path(summary_path),
        "scripts_scanned": summary["scripts_scanned"],
        "scripts_scanned_python": summary["scripts_scanned_python"],
        "scripts_scanned_shell": summary["scripts_scanned_shell"],
        "findings": summary["findings"],
        "violations": summary["violations"],
        "unreviewed": summary["unreviewed"],
        "compliant": summary["compliant"],
        "use_class_counts": summary["use_class_counts"],
        "files_with_violations": summary["files_with_violations"],
        "registry_discovery_layer_absent": (
            None if registry_report is None else registry_report["discovery_layer_absent"]
        ),
        "registry_sha256": None if registry_report is None else registry_report["sha256"],
    }, indent=2, ensure_ascii=False))
    return 1 if summary["violations"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
