#!/usr/bin/env python3
"""Combine BEAST chain logs, remove burnin, and compute ESS.

Usage:
    mamba run --live-stream -n bdbv python bin/postprocess_logs.py \
        --logs run1.log run2.log \
        --burnin 0.1 \
        --outdir output/ \
        --prefix analysis
"""

import argparse
import math
import sys
from pathlib import Path


def parse_log(path):
    """Parse a BEAST log file. Returns (header_list, rows_as_list_of_floats)."""
    headers = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if headers is None:
                headers = line.split('\t')
                continue
            vals = line.split('\t')
            try:
                rows.append([float(v) for v in vals])
            except ValueError:
                continue
    return headers, rows


def compute_ess(values):
    """Compute effective sample size using autocorrelation."""
    n = len(values)
    if n < 10:
        return float(n)

    mean = sum(values) / n
    centered = [v - mean for v in values]

    var = sum(c * c for c in centered) / n
    if var == 0:
        return float(n)

    max_lag = min(n - 1, 2000)
    act = 0.0
    for lag in range(1, max_lag + 1):
        autocorr = sum(centered[i] * centered[i + lag] for i in range(n - lag)) / (n * var)
        if autocorr < 0.05:
            break
        act += autocorr

    act = 1.0 + 2.0 * act
    if act < 1.0:
        act = 1.0

    return n / act


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--logs', nargs='+', required=True,
                        help='One or more BEAST .log files to combine')
    parser.add_argument('--burnin', type=float, default=0.1,
                        help='Fraction of samples to discard as burnin (default: 0.1)')
    parser.add_argument('--outdir', required=True, type=Path,
                        help='Output directory')
    parser.add_argument('--prefix', required=True,
                        help='Output file prefix')
    parser.add_argument('--min-ess', type=float, default=200,
                        help='Minimum acceptable ESS for key parameters (default: 200)')
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    all_headers = None
    all_rows = []

    for log_path in args.logs:
        print(f"Reading {log_path}")
        headers, rows = parse_log(log_path)
        if headers is None or not rows:
            print(f"  WARNING: no data in {log_path}, skipping", file=sys.stderr)
            continue

        if all_headers is None:
            all_headers = headers
        elif headers != all_headers:
            print(f"  WARNING: headers in {log_path} differ from first file", file=sys.stderr)
            print(f"  Expected: {all_headers[:5]}...", file=sys.stderr)
            print(f"  Got:      {headers[:5]}...", file=sys.stderr)
            sys.exit(1)

        n_burnin = int(len(rows) * args.burnin)
        kept = rows[n_burnin:]
        print(f"  {len(rows)} samples, removed {n_burnin} burnin, kept {len(kept)}")
        all_rows.extend(kept)

    if not all_rows:
        print("ERROR: no data after burnin removal", file=sys.stderr)
        sys.exit(1)

    print(f"\nCombined: {len(all_rows)} samples from {len(args.logs)} chain(s)")

    combined_path = args.outdir / f'{args.prefix}.combined.log'
    with open(combined_path, 'w') as f:
        f.write('\t'.join(all_headers) + '\n')
        for row in all_rows:
            f.write('\t'.join(f'{v:g}' for v in row) + '\n')
    print(f"Wrote {combined_path}")

    KEY_PARAMS = {'posterior', 'likelihood', 'prior', 'treelength', 'clock.rate',
                   'growthrate', 'doublingtime', 'logistic.t50'}
    KEY_PREFIXES = ('treemodel.rootheight', 'skygrid.', 'skyline.')

    def is_key_param(name):
        low = name.lower()
        if low in KEY_PARAMS:
            return True
        return any(low.startswith(p) for p in KEY_PREFIXES)

    ess_results = []
    for i, param in enumerate(all_headers):
        if param.lower() == 'state':
            continue
        values = [row[i] for row in all_rows]
        ess = compute_ess(values)
        ess_results.append((param, ess, min(values), max(values), is_key_param(param)))

    ess_path = args.outdir / f'{args.prefix}.ess_report.txt'
    n_low = 0
    n_low_key = 0
    low_key_params = []
    with open(ess_path, 'w') as f:
        f.write(f'ESS Report\n')
        f.write(f'Chains: {len(args.logs)}\n')
        f.write(f'Combined samples: {len(all_rows)}\n')
        f.write(f'Burnin: {args.burnin:.0%}\n')
        f.write(f'Min ESS threshold: {args.min_ess}\n')
        f.write(f'\n{"Parameter":<50} {"ESS":>10} {"Min":>15} {"Max":>15}  Key?\n')
        f.write('-' * 98 + '\n')
        for param, ess, lo, hi, is_key in ess_results:
            flag = ''
            if ess < args.min_ess:
                n_low += 1
                if is_key:
                    flag = ' *** KEY'
                    n_low_key += 1
                    low_key_params.append((param, ess))
                else:
                    flag = ' *'
            key_marker = '  Y' if is_key else ''
            f.write(f'{param:<50} {ess:>10.1f} {lo:>15.6g} {hi:>15.6g}{key_marker}{flag}\n')

        f.write(f'\n{n_low} parameter(s) below ESS {args.min_ess}')
        if n_low_key:
            f.write(f' ({n_low_key} key)')
        f.write('\n')
    print(f"Wrote {ess_path}")

    print(f"\nESS summary:")
    for param, ess, lo, hi, is_key in ess_results:
        if ess < args.min_ess:
            tag = " [KEY]" if is_key else ""
            print(f"  WARNING: {param} ESS = {ess:.1f} (< {args.min_ess}){tag}")
    if n_low == 0:
        print(f"  All parameters have ESS >= {args.min_ess}")
    else:
        print(f"  {n_low} parameter(s) with ESS < {args.min_ess}")

    if n_low_key > 0:
        print(f"\nERROR: {n_low_key} key parameter(s) have ESS < {args.min_ess}:")
        for param, ess in low_key_params:
            print(f"  {param}: {ess:.1f}")
        print(f"\nRun more chains or increase chain length, then re-run postprocessing.")
        print(f"See {ess_path} for full report.")
        sys.exit(2)


if __name__ == '__main__':
    main()
