#!/usr/bin/env python3
import argparse
import csv
import glob
import math
import os
import re
import statistics


def parse_elapsed_seconds(raw):
    s = raw.strip()
    if not s:
        return math.nan

    # GNU time may emit d-hh:mm:ss or hh:mm:ss or mm:ss
    days = 0
    if "-" in s:
        day_part, rest = s.split("-", 1)
        if day_part.isdigit():
            days = int(day_part)
            s = rest

    parts = s.split(":")
    try:
        if len(parts) == 3:
            h, m, sec = parts
            return days * 86400 + int(h) * 3600 + int(m) * 60 + float(sec)
        if len(parts) == 2:
            m, sec = parts
            return days * 86400 + int(m) * 60 + float(sec)
        return days * 86400 + float(s)
    except ValueError:
        return math.nan


def extract_first(pattern, text, flags=0):
    m = re.search(pattern, text, flags)
    return m.group(1).strip() if m else None


def extract_sample(text, fallback_name):
    # If logged by your script
    from_processing = extract_first(r"Processing sample:\s*(.+)", text)
    if from_processing:
        return os.path.basename(from_processing)

    # If GNU time has command with -i <sample>
    from_cmd = extract_first(r"Command being timed:\s*\".*?\s-i\s+([^\s\"]+)", text)
    if from_cmd:
        return os.path.basename(from_cmd)

    return fallback_name


def parse_file(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        txt = fh.read()

    sample = extract_sample(txt, os.path.basename(path))

    elapsed_raw = extract_first(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(.+)", txt)
    elapsed_s = parse_elapsed_seconds(elapsed_raw) if elapsed_raw else math.nan

    max_rss_raw = extract_first(r"Maximum resident set size \(kbytes\):\s*(\d+)", txt)
    max_rss_kb = float(max_rss_raw) if max_rss_raw else math.nan
    max_rss_gib = (max_rss_kb / 1024 / 1024) if not math.isnan(max_rss_kb) else math.nan

    user_raw = extract_first(r"User time \(seconds\):\s*([0-9]*\.?[0-9]+)", txt)
    system_raw = extract_first(r"System time \(seconds\):\s*([0-9]*\.?[0-9]+)", txt)
    user_s = float(user_raw) if user_raw else math.nan
    system_s = float(system_raw) if system_raw else math.nan

    cpu_time_s = user_s + system_s if not (math.isnan(user_s) or math.isnan(system_s)) else math.nan
    cpu_equiv = cpu_time_s / elapsed_s if not (math.isnan(cpu_time_s) or math.isnan(elapsed_s) or elapsed_s == 0) else math.nan

    return {
        "sample": sample,
        "elapsed_s": elapsed_s,
        "max_rss_kb": max_rss_kb,
        "max_rss_gib": max_rss_gib,
        "user_s": user_s,
        "system_s": system_s,
        "cpu_time_s": cpu_time_s,
        "cpu_equiv": cpu_equiv,
    }


def mean_ignore_nan(values):
    valid = [v for v in values if not math.isnan(v)]
    return sum(valid) / len(valid) if valid else math.nan


def percentile_ignore_nan(values, p):
    valid = sorted(v for v in values if not math.isnan(v))
    if not valid:
        return math.nan
    if len(valid) == 1:
        return valid[0]

    # Linear interpolation between closest ranks.
    idx = (len(valid) - 1) * (p / 100.0)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return valid[lo]
    frac = idx - lo
    return valid[lo] * (1.0 - frac) + valid[hi] * frac


def main():
    parser = argparse.ArgumentParser(description="Parse sc_Allelome.PRO2.snRNAseq*.err and compute runtime/memory metrics.")
    parser.add_argument(
        "--pattern",
        default="out_error/sc_Allelome.PRO2.snRNAseq*.err",
        help="Glob pattern for .err files (default: out_error/sc_Allelome.PRO2.snRNAseq*.err)",
    )
    parser.add_argument(
        "--out",
        default="sc_Allelome.PRO2.metrics.csv",
        help="Output CSV path (default: sc_Allelome.PRO2.metrics.csv)",
    )
    args = parser.parse_args()

    files = sorted(glob.glob(args.pattern))
    if not files:
        raise SystemExit(f"No files matched pattern: {args.pattern}")

    rows = [parse_file(p) for p in files]

    fieldnames = [
        "sample",
        "elapsed_s",
        "max_rss_kb",
        "max_rss_gib",
        "user_s",
        "system_s",
        "cpu_time_s",
        "cpu_equiv",
    ]

    with open(args.out, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Wrote {len(rows)} rows to {args.out}")
    print("Averages (ignoring missing values):")
    for k in fieldnames[1:]:
        values = [r[k] for r in rows]
        avg = mean_ignore_nan(values)
        med_vals = [v for v in values if not math.isnan(v)]
        med = statistics.median(med_vals) if med_vals else math.nan
        p95 = percentile_ignore_nan(values, 95)
        print(f"  {k}: {avg}")
        print(f"  {k}_median: {med}")
        print(f"  {k}_p95: {p95}")


if __name__ == "__main__":
    main()
