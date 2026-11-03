# tsinfer performance benchmark

Benchmark ancestor generation, ancestor matching and sample matching across
built tsinfer checkouts. The example covers four simulation sample sizes and
1KGP chr20p, with one repetition per dataset and version.

From this directory:

```sh
uv sync --locked
uv run snakemake --cores 196 --resources mem_mb=950000 --printshellcmds
```

The workflow defaults to `config.yaml.example`, the only tracked config. Edit
its paths, or copy it to `config.yaml` and pass `--configfile config.yaml`.
Each dataset supplies a fixed `zarr_path`, independent haploid `sample_size`,
output `prefix`, ancestral state and mask array names, and `num_reps`.
Reps run from 0 to `num_reps - 1` on the same input for every version.
Prefix and sample size together identify a dataset. Output filename `n` counts
haplotypes; simulation input filename `n` counts diploid individuals.

Optional per-dataset `sequence_length` overrides the Zarr contig length. The
simulations retain 28,100,000 bp because their stored contig can describe the
full chromosome. Null or omitted uses the stored length; 1KGP keeps 64,444,167 bp
with masks selecting chr20p. Mask values of true exclude entries. Per-dataset
`precision` preserves the simulation default (null) and 1KGP precision (15).

The dependencies and lockfile come from tsinfer-anc-eval and use Zarr 2, required
by these releases. Build each configured tsinfer checkout for the environment's
Python and NumPy. Workers use the Python running Snakemake and select checkouts
through `PYTHONPATH`. To reuse an existing anc-eval environment, replace
`uv run` above with `uv run --project /path/to/tsinfer-anc-eval --no-sync`.

Run within an exclusive allocation, keeping command resource totals equal to
the config. Each job reserves all configured cores and memory, so stages run
sequentially in fresh processes. Memory is a scheduling reservation; `time_min`
is cluster submission metadata. Sample matching omits post-processing.

Outputs under `data/` (or `output_dir`) use `<prefix>-n<n>-rep<rep>-v<version>-<suffix>`:

- `ancestors/`: `ancestors.zarr` (a Zarr-backed LMDB file) and
  `ancestors.trees`.
- `trees/`: `raw.trees` from sample matching.
- `.progress/`: `<stage>.log`.
- `dataframes/`: `perf.csv`, with one row per inference stage.
- `dataframes/perf-v<version>.csv` and `dataframes/perf_combined.csv`: combined
  tables retaining every repetition.

Version labels are configured as `"0.4"` and `"0.5"`, giving filenames such as
`1kgp-chr20-L0-R2.81e7-n6404-rep0-v0.5-raw.trees`.

All performance tables include `prefix`, `sample_size`, and `rep`.
`num_samples` is the observed haploid count and equals `sample_size`.
Times are seconds and peak memory is bytes. API timers exclude startup, input
loading and output writing; peak memory includes input loading. Configured
`version` and original `provenance_version` are both retained. Choose a new
output directory after changing scientific settings.

Optional focal-site spans, using either built legacy checkout:

```sh
PYTHONPATH=/path/to/built/tsinfer uv run python scripts/extract_spans.py data
```

This writes `data/dataframes/ancestor_spans.csv` with the same run identity
columns. Spans measure the distance between the first and last included
inference sites. Existing grouping analysis tables are retained in
`data/dataframes/analysis/`.
