#!/usr/bin/env python3
"""BEAST phylodynamics pipeline: XML generation → postprocessing → figures → report.

Usage:
    mamba run --live-stream -n bdbv python bin/run.py \
        --alignment /path/to/final.aln \
        --outdir output/20261106 \
        --prefix BDBV2026.drc \
        --template-sg templates/GTR_SG_MAT.template.xml \
        --anchor-date 2025-11-08

    # Postprocess existing BEAST output:
    mamba run --live-stream -n bdbv python bin/run.py \
        --beast-dir /path/to/beast/output \
        --outdir output/20261106 \
        --prefix BDBV2026.drc \
        --skip-to postprocess_logs

    # Rerun figures with updated epi data:
    mamba run --live-stream -n bdbv python bin/run.py \
        --beast-dir /path/to/beast/output \
        --linelist /path/to/linelist.xlsx \
        --case-data /path/to/status_aggregated.csv \
        --outdir output/20261106 \
        --prefix BDBV2026.drc \
        --skip-to generate_figures
"""

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent
PYTHON = sys.executable

STEPS = [
    'generate_xmls',
    'postprocess_logs',
    'postprocess_trees',
    'generate_figures',
    'generate_report',
]


def run(cmd, description=None, allow_exit_2=False):
    if description:
        print(f"\n{'='*60}")
        print(f"  {description}")
        print(f"{'='*60}")
    print(f"$ {cmd}\n")
    result = subprocess.run(cmd, shell=True)
    if result.returncode == 2 and allow_exit_2:
        return 2
    if result.returncode != 0:
        print(f"\nERROR: command failed with return code {result.returncode}")
        sys.exit(1)
    return 0


def count_sequences(fasta_path):
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def should_run(step, skip_to, stop_after):
    """Check whether a step should run given --skip-to and --stop-after."""
    idx = STEPS.index(step)
    if skip_to:
        skip_idx = STEPS.index(skip_to)
        if idx < skip_idx:
            return False
    if stop_after:
        stop_idx = STEPS.index(stop_after)
        if idx > stop_idx:
            return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--alignment', default=None,
                        help='Final alignment from pipeline 1 (required for step 1)')
    parser.add_argument('--beast-dir', default=None,
                        help='Directory with BEAST output logs/trees (required for steps 2-3)')
    parser.add_argument('--outdir', required=True, type=Path,
                        help='Output directory')
    parser.add_argument('--prefix', required=True,
                        help='Output file prefix')

    parser.add_argument('--template-sg', default=None, help='SkyGrid template XML')
    parser.add_argument('--template-egc', default=None, help='EGC template XML')
    parser.add_argument('--template-lgc', default=None, help='LGC template XML')
    parser.add_argument('--template-mat', default=None, help='MAT SkyGrid template XML')
    parser.add_argument('--anchor-date', default='2025-11-08', help='SkyGrid anchor date')
    parser.add_argument('--chain-length', type=int, default=500_000_000,
                        help='MCMC chain length (default: 500M)')
    parser.add_argument('--log-every', type=int, default=50_000,
                        help='Log sampling interval (default: 50000)')

    parser.add_argument('--burnin', type=float, default=0.1,
                        help='Burnin fraction (default: 0.1)')
    parser.add_argument('--resample', type=int, default=1000,
                        help='Target tree count after resampling (default: 1000)')
    parser.add_argument('--treeannotator', default='treeannotator',
                        help='Path to TreeAnnotator executable')
    parser.add_argument('--min-ess', type=float, default=200,
                        help='Minimum acceptable ESS for key parameters (default: 200)')
    parser.add_argument('--force-low-ess', action='store_true',
                        help='Continue past ESS check even if key parameters are below threshold')

    parser.add_argument('--linelist', default=None, help='Genomic linelist Excel')
    parser.add_argument('--sitrep', default=None, help='Sitrep CSV')
    parser.add_argument('--case-data', default=None, help='status_aggregated CSV')

    parser.add_argument('--skip-to', choices=STEPS, default=None,
                        help='Skip to this step')
    parser.add_argument('--stop-after', choices=STEPS, default=None,
                        help='Stop after this step')

    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Pipeline part 2: BEAST phylodynamics")
    print(f"Output: {args.outdir}")
    print(f"Prefix: {args.prefix}")
    if args.skip_to:
        print(f"Skipping to: {args.skip_to}")
    if args.stop_after:
        print(f"Stopping after: {args.stop_after}")
    print()

    # ------------------------------------------------------------------
    # Step 1: Generate BEAST XMLs
    # ------------------------------------------------------------------
    if should_run('generate_xmls', args.skip_to, args.stop_after):
        print(f"\n── Step 1: Generate BEAST XMLs ── [{datetime.now().strftime('%H:%M:%S')}]")

        templates = []
        if args.template_sg:
            templates.append(('--template-sg', args.template_sg))
        if args.template_egc:
            templates.append(('--template-egc', args.template_egc))
        if args.template_lgc:
            templates.append(('--template-lgc', args.template_lgc))
        if args.template_mat:
            templates.append(('--template-mat', args.template_mat))

        if not templates:
            print("  No templates specified, skipping XML generation")
        elif not args.alignment:
            print("  ERROR: --alignment required for XML generation", file=sys.stderr)
            sys.exit(1)
        else:
            n = count_sequences(args.alignment)
            beast_prefix = f'{args.prefix}.n{n}'

            xml_script = SCRIPTS_DIR / 'generate_beast_xml.py'
            cmd = (
                f'{PYTHON} {xml_script}'
                f' -a {args.alignment}'
                f' -o {args.outdir}'
                f' -p {beast_prefix}'
                f' --chain-length {args.chain_length}'
                f' --log-every {args.log_every}'
            )
            for flag, path in templates:
                cmd += f' {flag} {path}'
                if flag in ('--template-sg', '--template-mat'):
                    cmd += f' --anchor-date {args.anchor_date}'
            run(cmd, f'Generating BEAST XMLs (n={n})')

    # ------------------------------------------------------------------
    # Step 2: Postprocess logs
    # ------------------------------------------------------------------
    if should_run('postprocess_logs', args.skip_to, args.stop_after):
        print(f"\n── Step 2: Postprocess logs ── [{datetime.now().strftime('%H:%M:%S')}]")

        beast_dir = Path(args.beast_dir) if args.beast_dir else args.outdir
        log_files = sorted(beast_dir.glob('*.log'))
        log_files = [f for f in log_files if not f.name.endswith('.rates.log')
                     and not f.name.endswith('.ess_report.txt')
                     and not f.name.endswith('.combined.log')]

        if not log_files:
            print("  No .log files found, skipping")
        else:
            log_script = SCRIPTS_DIR / 'postprocess_logs.py'
            logs_arg = ' '.join(str(f) for f in log_files)
            rc = run(
                f'{PYTHON} {log_script}'
                f' --logs {logs_arg}'
                f' --burnin {args.burnin}'
                f' --outdir {args.outdir}'
                f' --prefix {args.prefix}'
                f' --min-ess {args.min_ess}',
                f'Postprocessing {len(log_files)} log file(s)',
                allow_exit_2=True
            )
            if rc == 2:
                if args.force_low_ess:
                    print("\n  --force-low-ess: continuing despite low ESS")
                else:
                    print("\n  Stopping: key parameters have low ESS.")
                    print("  Options:")
                    print("    1. Run more BEAST chains, then re-run postprocessing")
                    print("    2. Pass --force-low-ess to continue anyway")
                    print(f"    3. Re-enter at step 3: --skip-to postprocess_trees --force-low-ess")
                    sys.exit(2)

    # ------------------------------------------------------------------
    # Step 3: Postprocess trees
    # ------------------------------------------------------------------
    if should_run('postprocess_trees', args.skip_to, args.stop_after):
        print(f"\n── Step 3: Postprocess trees ── [{datetime.now().strftime('%H:%M:%S')}]")

        beast_dir = Path(args.beast_dir) if args.beast_dir else args.outdir
        tree_files = sorted(beast_dir.glob('*.trees'))
        tree_files = [f for f in tree_files if 'burninRemoved' not in f.name
                      and 'resampled' not in f.name
                      and 'history' not in f.name
                      and 'default' not in f.name
                      and 'empirical' not in f.name]

        if not tree_files:
            print("  No .trees files found, skipping")
        else:
            trees_script = SCRIPTS_DIR / 'postprocess_trees.py'
            trees_arg = ' '.join(str(f) for f in tree_files)
            cmd = (
                f'{PYTHON} {trees_script}'
                f' --trees {trees_arg}'
                f' --burnin {args.burnin}'
                f' --resample {args.resample}'
                f' --outdir {args.outdir}'
                f' --prefix {args.prefix}'
                f' --treeannotator {args.treeannotator}'
            )
            run(cmd, f'Postprocessing {len(tree_files)} tree file(s)')

    # ------------------------------------------------------------------
    # Step 4: Generate figures (TODO)
    # ------------------------------------------------------------------
    if should_run('generate_figures', args.skip_to, args.stop_after):
        print(f"\n── Step 4: Generate figures ── [{datetime.now().strftime('%H:%M:%S')}]")
        print("  TODO: figure generation not yet implemented")
        if args.linelist:
            print(f"  Linelist: {args.linelist}")
        if args.sitrep:
            print(f"  Sitrep: {args.sitrep}")
        if args.case_data:
            print(f"  Case data: {args.case_data}")

    # ------------------------------------------------------------------
    # Step 5: Generate report (TODO)
    # ------------------------------------------------------------------
    if should_run('generate_report', args.skip_to, args.stop_after):
        print(f"\n── Step 5: Generate report ── [{datetime.now().strftime('%H:%M:%S')}]")
        print("  TODO: report generation not yet implemented")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"  Pipeline part 2 complete")
    print(f"{'='*60}")
    print(f"Output directory: {args.outdir}")
    for f in sorted(args.outdir.iterdir()):
        if f.is_file():
            print(f"  {f.name}")


if __name__ == '__main__':
    main()
