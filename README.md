# BDBV Phylodynamics Pipeline (Part 2) — Nextflow

**BEAST XML generation and postprocessing for Bundibugyo ebolavirus (BDBV) phylodynamics**

A [Nextflow](https://www.nextflow.io/) port of [`bdbv-pipeline-part2`](https://github.com/rmcolq/bdbv-pipeline-part2). It picks up where [`bdbv-pipeline-nf`](https://github.com/rmcolq/bdbv-pipeline-nf) (Part 1) leaves off, taking its final alignment as input. BEAST itself is run outside this pipeline, so there are two stages, used one after the other.

## Pipeline overview

1. **Generate BEAST XMLs** *(`generate_xmls`, with `--alignment`)* — swap the alignment into the templates, update the SkyGrid parameters, adjust operator weights and set the chain length and log interval. One XML per template.
2. *Run BEAST on the XMLs yourself.*
3. **Postprocess logs** *(`postprocess_logs`, with `--beast_dir`)* — combine chain logs, remove burnin and check ESS (**ESS gate**).
4. **Postprocess trees** *(`postprocess_trees`)* — remove burnin, resample and build an MCC tree with TreeAnnotator.

Figures and report generation are not yet implemented in the Python version, so they are not part of this pipeline either.

## Installation

Install [Nextflow](https://www.nextflow.io/) (≥ 23.04) and [Docker](https://www.docker.com/products/docker-desktop), then clone this repository:

```bash
git clone https://github.com/rmcolq/bdbv-pipeline-beast-nf.git
cd bdbv-pipeline-beast-nf
```

The pipeline pulls the `articnetworkorg/bdbv-pipeline-beast` image automatically (built from `Dockerfile` by `.github/workflows/docker.yml`). Other options: `-profile singularity`, `-profile conda` (uses `environment.yml`), or `--wf.container_version TAG` for your own build.

## Running the pipeline

**Generate BEAST XMLs** from the final alignment of Part 1:

```bash
nextflow run . \
    --alignment /path/to/BDBV2026.drc.ADAR_masked.final.aln \
    --outdir output/20261106 \
    --prefix BDBV2026.drc \
    --anchor_date 2025-11-08
```

The XMLs are written to `beast_xml/`. By default one XML is made from each of the SkyGrid and EGC templates. Use `--template_mat templates/GTR_SG_MAT.template.xml` for the MAT SkyGrid template, or `--template_sg false` to skip one.

**Postprocess BEAST output** once the runs have finished:

```bash
nextflow run . \
    --beast_dir /path/to/beast/output \
    --outdir output/20261106 \
    --prefix BDBV2026.drc \
    --burnin 0.1
```

### ESS gate

After combining the logs the pipeline checks ESS for the key parameters (posterior, likelihood, root height, clock rate, SkyGrid parameters, etc.). If any is below `--min_ess` (default 200) it stops before the trees and tells you how to proceed. Review `postprocess_logs/<prefix>.ess_report.txt`, then either:

1. Run more BEAST chains, add their logs and trees to `--beast_dir`, and rerun
2. Continue anyway: rerun with `--force_low_ess -resume`

## Output files

| Folder | Contents |
|--------|----------|
| `beast_xml/` | Generated BEAST XMLs |
| `postprocess_logs/` | Combined log, ESS report and `ess_status.txt` (`pass` or `low`) |
| `postprocess_trees/` | Burnin-removed trees, resampled trees and MCC tree |
| `execution/` | Nextflow timeline, report and trace |

## Parameters

See `nextflow_schema.json` for the full definitions.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--alignment` | Final alignment from Part 1; generates BEAST XMLs | |
| `--beast_dir` | Folder of BEAST `.log` and `.trees` files; runs postprocessing | |
| `--outdir` | Output folder | `output` |
| `--prefix` | Output file prefix | `bdbv_<run_name>` |
| `--run_name` | Label used in the default prefix | today, `YYYYMMDD` |
| `--template_sg`, `--template_mat`, `--template_egc`, `--template_lgc` | Templates; one XML per template given. Use `false` to skip a default | SG and EGC |
| `--anchor_date` | SkyGrid anchor date | `2025-11-08` |
| `--chain_length` | MCMC chain length | `500000000` |
| `--log_every` | Log sampling interval | `50000` |
| `--burnin` | Burnin fraction | `0.1` |
| `--resample` | Target tree count after resampling | `1000` |
| `--min_ess` | Minimum ESS for key parameters | `200` |
| `--force_low_ess` | Continue past a low ESS | `false` |
| `--treeannotator` | TreeAnnotator executable | `treeannotator` |

## Running the Python scripts without Nextflow

The scripts in `bin/` still work on their own (`bin/run.py` is the original single entry point). They only need Python ≥ 3.10, plus [TreeAnnotator](https://beast.community/treeannotator) on `PATH` for the MCC tree; `conda env create -f environment.yml` provides both.

## Data files

- **`data/aliases.csv`** — HZ name canonicalization (copied from bdbv-pipeline)
- **`data/health_area_zone_province_slim.csv`** — HZ → province lookup
- **`templates/`** — BEAST XML templates (GTR_SG, GTR_SG_MAT, GTR_EGC, GTR_LGC)

## Relationship to other repos

- **`bdbv-pipeline-nf`** (Part 1) — Download, clean, align, build ML tree
- **`bdbv-pipeline-part2`** — the original Python-only version of this pipeline
- **`2026_BDBV_phylogeography`** — Phylogeographic analysis (DTA/BSSVS on empirical trees)
