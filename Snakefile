"""Run each benchmark repetition as three separately measured inference stages."""

import pathlib
import re
import shlex
import sys

configfile: "config.yaml.example"

output_dir = config["output_dir"]
config_path = workflow.configfiles[-1]
threads_per_job = config["resources"]["threads"]
memory_per_job = config["resources"]["mem_mb"]
time_per_job = config["resources"]["time_min"]
worker = pathlib.Path(workflow.basedir) / "anc_perf.py"

# A prefix can be shared by different sample sizes, but each pair is unique.
datasets = {}
for dataset in config["datasets"]:
    key = (dataset["prefix"], dataset["sample_size"])
    if key in datasets:
        raise ValueError(f"Duplicate dataset: {key}")
    if dataset["sample_size"] < 1 or dataset["num_reps"] < 1:
        raise ValueError(f"sample_size and num_reps must be positive: {key}")
    datasets[key] = dataset

run_name = "{prefix}-n{sample_size}-rep{rep}-v{version}"
ancestor_store = f"{output_dir}/ancestors/{run_name}-ancestors.zarr"
ancestor_trees = f"{output_dir}/ancestors/{run_name}-ancestors.trees"
raw_trees = f"{output_dir}/trees/{run_name}-raw.trees"
performance = f"{output_dir}/dataframes/{run_name}-perf.csv"
progress = f"{output_dir}/.progress/{run_name}"

wildcard_constraints:
    version="|".join(re.escape(version) for version in config["versions"]),
    prefix="|".join(re.escape(dataset["prefix"]) for dataset in config["datasets"]),
    sample_size=r"\d+",
    rep=r"\d+",

shell.prefix(
    "set -euo pipefail; "
    f"export OMP_NUM_THREADS={config['helper_threads']} "
    f"OPENBLAS_NUM_THREADS={config['helper_threads']} "
    f"MKL_NUM_THREADS={config['helper_threads']} "
    f"NUMEXPR_NUM_THREADS={config['helper_threads']}; "
)


def zarr_path(wildcards):
    dataset = datasets[(wildcards.prefix, int(wildcards.sample_size))]
    if int(wildcards.rep) >= dataset["num_reps"]:
        raise ValueError(f"Repetition outside num_reps: {wildcards.rep}")
    return dataset["zarr_path"]


def runner(version):
    installation = config["versions"][version]
    return (
        f"PYTHONPATH={shlex.quote(installation)} "
        f"{shlex.quote(sys.executable)} {shlex.quote(str(worker))}"
    )


def performance_inputs(wildcards):
    paths = []
    for dataset in config["datasets"]:
        for rep in range(dataset["num_reps"]):
            path = performance.format(
                version=wildcards.version,
                prefix=dataset["prefix"],
                sample_size=dataset["sample_size"],
                rep=rep,
            )
            paths.append(path)
    return paths


rule all:
    input:
        f"{output_dir}/dataframes/perf_combined.csv",


rule generate_ancestors:
    input:
        zarr=zarr_path,
    output:
        ancestor_store,
    log:
        f"{progress}-generate_ancestors.log",
    threads: threads_per_job
    resources:
        mem_mb=memory_per_job,
        time_min=time_per_job,
    params:
        runner=lambda wildcards: runner(wildcards.version),
        config=config_path,
    shell:
        "{params.runner} generate_ancestors --config {params.config:q} "
        "--output {output:q} --version {wildcards.version:q} "
        "--prefix {wildcards.prefix:q} --sample-size {wildcards.sample_size} "
        "--rep {wildcards.rep} --threads {threads} > {log:q} 2>&1"


rule match_ancestors:
    input:
        zarr=zarr_path,
        ancestors=ancestor_store,
    output:
        ancestor_trees,
    log:
        f"{progress}-match_ancestors.log",
    threads: threads_per_job
    resources:
        mem_mb=memory_per_job,
        time_min=time_per_job,
    params:
        runner=lambda wildcards: runner(wildcards.version),
        config=config_path,
    shell:
        "{params.runner} match_ancestors --config {params.config:q} "
        "--input {input.ancestors:q} --output {output:q} "
        "--version {wildcards.version:q} --prefix {wildcards.prefix:q} "
        "--sample-size {wildcards.sample_size} --rep {wildcards.rep} "
        "--threads {threads} > {log:q} 2>&1"


rule match_samples:
    input:
        zarr=zarr_path,
        ancestors=ancestor_trees,
    output:
        raw_trees,
    log:
        f"{progress}-match_samples.log",
    threads: threads_per_job
    resources:
        mem_mb=memory_per_job,
        time_min=time_per_job,
    params:
        runner=lambda wildcards: runner(wildcards.version),
        config=config_path,
    shell:
        "{params.runner} match_samples --config {params.config:q} "
        "--input {input.ancestors:q} --output {output:q} "
        "--version {wildcards.version:q} --prefix {wildcards.prefix:q} "
        "--sample-size {wildcards.sample_size} --rep {wildcards.rep} "
        "--threads {threads} > {log:q} 2>&1"


rule extract_performance:
    input:
        raw_trees,
    output:
        performance,
    threads: threads_per_job
    resources:
        mem_mb=memory_per_job,
    params:
        runner=lambda wildcards: runner(wildcards.version),
        config=config_path,
    shell:
        "{params.runner} extract --config {params.config:q} --input {input:q} "
        "--output {output:q} --version {wildcards.version:q} "
        "--prefix {wildcards.prefix:q} --sample-size {wildcards.sample_size} "
        "--rep {wildcards.rep} --threads {threads}"


rule aggregate_version:
    input:
        performance_inputs,
    output:
        f"{output_dir}/dataframes/perf-v{{version}}.csv",
    threads: threads_per_job
    resources:
        mem_mb=memory_per_job,
    params:
        runner=lambda wildcards: runner(wildcards.version),
        config=config_path,
    shell:
        "{params.runner} aggregate --config {params.config:q} "
        "--input {input:q} --output {output:q}"


rule aggregate_all:
    input:
        expand(
            f"{output_dir}/dataframes/perf-v{{version}}.csv", version=config["versions"]
        ),
    output:
        f"{output_dir}/dataframes/perf_combined.csv",
    threads: threads_per_job
    resources:
        mem_mb=memory_per_job,
    params:
        runner=runner(next(iter(config["versions"]))),
        config=config_path,
    shell:
        "{params.runner} aggregate --config {params.config:q} "
        "--input {input:q} --output {output:q}"
