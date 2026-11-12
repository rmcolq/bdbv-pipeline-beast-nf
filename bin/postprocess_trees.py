#!/usr/bin/env python3
"""Remove burnin from BEAST .trees files, resample, and generate MCC tree.

Usage:
    mamba run --live-stream -n bdbv python bin/postprocess_trees.py \
        --trees run1.trees run2.trees \
        --burnin 0.1 \
        --resample 1000 \
        --outdir output/ \
        --prefix analysis
"""

import argparse
import math
import subprocess
import sys
from pathlib import Path


def parse_nexus_trees(path, burnin_fraction):
    """Parse a NEXUS .trees file. Returns (preamble_lines, tree_lines, n_total, n_removed).

    preamble_lines: everything before the first 'tree ' line (NEXUS header, translate block)
    tree_lines: post-burnin tree lines
    """
    preamble = []
    trees = []
    in_trees_block = False
    past_translate = False

    with open(path) as f:
        for line in f:
            stripped = line.strip()

            if stripped.lower() == 'begin trees;':
                in_trees_block = True
                preamble.append(line)
                continue

            if not in_trees_block:
                preamble.append(line)
                continue

            if stripped.lower().startswith('tree '):
                past_translate = True
                trees.append(line)
            elif not past_translate:
                preamble.append(line)
            elif stripped.lower() in ('end;', 'end'):
                pass
            else:
                preamble.append(line)

    n_total = len(trees)
    n_burnin = int(n_total * burnin_fraction)
    kept = trees[n_burnin:]

    return preamble, kept, n_total, n_burnin


def resample_trees(tree_lines, target):
    """Evenly resample tree_lines to target count."""
    n = len(tree_lines)
    if n <= target:
        return tree_lines

    step = n / target
    indices = [int(i * step) for i in range(target)]
    return [tree_lines[i] for i in indices]


def write_nexus_trees(path, preamble_lines, tree_lines):
    """Write a NEXUS trees file."""
    with open(path, 'w') as f:
        for line in preamble_lines:
            f.write(line)
        for line in tree_lines:
            f.write(line)
        f.write('End;\n')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--trees', nargs='+', required=True,
                        help='One or more BEAST .trees files')
    parser.add_argument('--burnin', type=float, default=0.1,
                        help='Fraction of trees to discard as burnin (default: 0.1)')
    parser.add_argument('--resample', type=int, default=1000,
                        help='Target number of trees after resampling (default: 1000)')
    parser.add_argument('--outdir', required=True, type=Path,
                        help='Output directory')
    parser.add_argument('--prefix', required=True,
                        help='Output file prefix')
    parser.add_argument('--treeannotator', default='treeannotator',
                        help='Path to TreeAnnotator executable (default: treeannotator)')
    parser.add_argument('--skip-mcc', action='store_true',
                        help='Skip MCC tree generation')
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    all_preamble = None
    all_trees = []

    for trees_path in args.trees:
        print(f"Reading {trees_path}")
        preamble, kept, n_total, n_burnin = parse_nexus_trees(trees_path, args.burnin)

        if all_preamble is None:
            all_preamble = preamble

        print(f"  {n_total} trees, removed {n_burnin} burnin, kept {len(kept)}")
        all_trees.extend(kept)

    if not all_trees:
        print("ERROR: no trees after burnin removal", file=sys.stderr)
        sys.exit(1)

    print(f"\nCombined: {len(all_trees)} trees from {len(args.trees)} chain(s)")

    burnin_path = args.outdir / f'{args.prefix}.burninRemoved.trees'
    write_nexus_trees(burnin_path, all_preamble, all_trees)
    print(f"Wrote {burnin_path} ({len(all_trees)} trees)")

    resampled = resample_trees(all_trees, args.resample)
    resampled_path = args.outdir / f'{args.prefix}.burninRemoved.resampled.trees'
    write_nexus_trees(resampled_path, all_preamble, resampled)
    print(f"Wrote {resampled_path} ({len(resampled)} trees)")

    if not args.skip_mcc:
        mcc_path = args.outdir / f'{args.prefix}.mcc.tree'
        print(f"\nRunning TreeAnnotator...")
        cmd = [
            args.treeannotator,
            '-burnin', '0',
            '-heights', 'median',
            str(burnin_path),
            str(mcc_path),
        ]
        print(f"$ {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"WARNING: TreeAnnotator failed (exit {result.returncode})", file=sys.stderr)
            if result.stderr:
                print(result.stderr, file=sys.stderr)
        else:
            print(f"Wrote {mcc_path}")


if __name__ == '__main__':
    main()
