#!/usr/bin/env python3
"""Dynamic server thread governance for PHB-GTDB launches.

Project rule (AGENTS.md): a server task may use at most
``min(40, total_cores - one_minute_load - 10)`` threads, measured immediately
before launch from ``nproc`` and ``/proc/loadavg``.  The reserve keeps roughly
ten cores on the host, and the hard cap keeps the shared server from being
oversubscribed no matter how idle it looks.

Rationale for replacing fixed defaults: a hard-coded ``70`` silently encodes
the assumption "the host has 80 cores and is completely idle".  Both halves of
that assumption are false in general, so the number must be *measured* and the
measurement must be recorded next to the run.

This module never launches anything.  It only computes the limit, records the
measurement, and exits non-zero when no capacity remains.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_RESERVE = 10
DEFAULT_HARD_CAP = 40
BASIS = "min(hard_cap, total_cores - load_1m - reserve)"
MEASUREMENT_SOURCE = "nproc+/proc/loadavg"
DEFAULT_CPUINFO = "/proc/cpuinfo"
DEFAULT_LOADAVG = "/proc/loadavg"


def compute_thread_limit(
    total_cores: int,
    load_1m: float,
    reserve: int = DEFAULT_RESERVE,
    hard_cap: int = DEFAULT_HARD_CAP,
) -> int:
    """Return the auditable per-task thread limit.

    Raises ``ValueError`` for a nonsensical measurement and ``RuntimeError``
    when the reserve leaves no capacity at all -- the caller must abort rather
    than overload the host.
    """
    if total_cores < 1:
        raise ValueError(f"invalid total_cores: {total_cores}")
    if load_1m < 0:
        raise ValueError(f"invalid load_1m: {load_1m}")
    if reserve < 0:
        raise ValueError(f"invalid reserve: {reserve}")
    if hard_cap < 1:
        raise ValueError(f"invalid hard_cap: {hard_cap}")
    available = int(total_cores - load_1m - reserve)
    if available < 1:
        raise RuntimeError(
            "insufficient free capacity after reserve: "
            f"total_cores={total_cores}, load_1m={load_1m}, reserve={reserve}"
        )
    return min(hard_cap, available)


def resource_record(
    *,
    total_cores: int,
    load_1m: float,
    reserve: int,
    hard_cap: int,
    threads: int,
    source: str = MEASUREMENT_SOURCE,
    measured_at_utc: str | None = None,
) -> dict:
    """Build the JSON record that binds a launch to its measurement."""
    return {
        "total_cores": total_cores,
        "load_1m": load_1m,
        "reserve": reserve,
        "hard_cap": hard_cap,
        "threads": threads,
        "basis": BASIS,
        "source": source,
        "measured_at_utc": measured_at_utc or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def parse_loadavg(text: str) -> float:
    """Read the one-minute load average from ``/proc/loadavg`` content."""
    fields = text.split()
    if len(fields) < 3:
        raise ValueError(f"malformed /proc/loadavg content: {text!r}")
    try:
        value = float(fields[0])
    except ValueError as exc:
        raise ValueError(f"malformed /proc/loadavg one-minute load: {fields[0]!r}") from exc
    if value < 0:
        raise ValueError(f"negative one-minute load: {value}")
    return value


def count_cores(cpuinfo_text: str) -> int:
    """Count logical processors from ``/proc/cpuinfo`` content (``nproc`` parity)."""
    processors = sum(1 for line in cpuinfo_text.splitlines() if line.startswith("processor"))
    if processors < 1:
        raise ValueError("no 'processor' entries found in /proc/cpuinfo")
    return processors


def read_measurements(
    cpuinfo_path: str | Path = DEFAULT_CPUINFO,
    loadavg_path: str | Path = DEFAULT_LOADAVG,
) -> tuple[int, float]:
    """Read total cores and one-minute load, failing loudly if unmeasurable."""
    try:
        cpuinfo_text = Path(cpuinfo_path).read_text(encoding="utf-8", errors="replace")
        loadavg_text = Path(loadavg_path).read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise RuntimeError(
            f"cannot measure server capacity: {cpuinfo_path} / {loadavg_path}: {exc}"
        ) from exc
    try:
        return count_cores(cpuinfo_text), parse_loadavg(loadavg_text)
    except ValueError as exc:
        raise RuntimeError(f"cannot measure server capacity: {exc}") from exc


def measure(
    *,
    total_cores: int | None = None,
    load_1m: float | None = None,
    reserve: int = DEFAULT_RESERVE,
    hard_cap: int = DEFAULT_HARD_CAP,
    requested: int | None = None,
    cpuinfo_path: str | Path = DEFAULT_CPUINFO,
    loadavg_path: str | Path = DEFAULT_LOADAVG,
) -> dict:
    """Measure (unless injected), compute the limit, and validate a request."""
    if total_cores is None or load_1m is None:
        measured_cores, measured_load = read_measurements(cpuinfo_path, loadavg_path)
        total_cores = measured_cores if total_cores is None else total_cores
        load_1m = measured_load if load_1m is None else load_1m
    threads = compute_thread_limit(total_cores, load_1m, reserve, hard_cap)
    if requested is not None:
        if requested < 1:
            raise ValueError(f"invalid requested thread count: {requested}")
        if requested > threads:
            raise RuntimeError(
                f"requested {requested} threads exceeds the measured limit of {threads} "
                f"(total_cores={total_cores}, load_1m={load_1m}, reserve={reserve}, hard_cap={hard_cap})"
            )
    return resource_record(
        total_cores=total_cores, load_1m=load_1m, reserve=reserve,
        hard_cap=hard_cap, threads=threads,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--total-cores", type=int, default=None, help="override nproc (testing only)")
    parser.add_argument("--load-1m", type=float, default=None, help="override /proc/loadavg (testing only)")
    parser.add_argument("--reserve", type=int, default=DEFAULT_RESERVE)
    parser.add_argument("--hard-cap", type=int, default=DEFAULT_HARD_CAP)
    parser.add_argument("--requested", type=int, default=None,
                        help="a proposed thread count; rejected if it exceeds the measured limit")
    parser.add_argument("--out", default=None, help="also write the JSON record to this path")
    parser.add_argument("--cpuinfo", default=DEFAULT_CPUINFO)
    parser.add_argument("--loadavg", default=DEFAULT_LOADAVG)
    args = parser.parse_args(argv)
    try:
        record = measure(
            total_cores=args.total_cores, load_1m=args.load_1m,
            reserve=args.reserve, hard_cap=args.hard_cap, requested=args.requested,
            cpuinfo_path=args.cpuinfo, loadavg_path=args.loadavg,
        )
    except (ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    payload = json.dumps(record, sort_keys=True)
    if args.out:
        out = Path(args.out)
        if out.exists():
            print(f"[ERROR] refusing to overwrite existing resource record: {out}", file=sys.stderr)
            return 2
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(payload + "\n", encoding="utf-8", newline="\n")
    print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
