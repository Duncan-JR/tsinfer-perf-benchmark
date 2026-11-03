"""Run one benchmark stage, or extract/aggregate its provenance metrics."""

import argparse
import json
import logging
import pathlib

import pandas as pd
import tsinfer
import tskit
import yaml

logger = logging.getLogger(__name__)


def get_dataset(config, prefix, sample_size):
    """Select the dataset identified by the output prefix and haploid count."""
    matches = []
    for dataset in config["datasets"]:
        if dataset["prefix"] == prefix and dataset["sample_size"] == sample_size:
            matches.append(dataset)
    if len(matches) != 1:
        raise ValueError(f"Expected one dataset for {prefix}, sample_size={sample_size}")
    return matches[0]


def extract_performance(path, version, prefix, sample_size, rep, num_threads):
    """Read the three inherited stage metrics and identify their benchmark run."""
    ts = tskit.load(path)
    if ts.num_samples != sample_size:
        raise ValueError(f"Expected {sample_size} haplotypes, found {ts.num_samples}")
    rows = []
    for prov in ts.provenances():
        record = json.loads(prov.record)
        if record["software"]["name"] == "tsinfer":
            command = record["parameters"]["command"]
            if command in ["generate_ancestors", "match_ancestors", "match_samples"]:
                row = dict(record["resources"])
                row.update(
                    command=command,
                    prefix=prefix,
                    sample_size=sample_size,
                    rep=rep,
                    num_samples=ts.num_samples,
                    num_sites=ts.num_sites,
                    num_mutations=ts.num_mutations,
                    num_edges=ts.num_edges,
                    version=version,
                    provenance_version=record["software"]["version"],
                    sequence_length=ts.sequence_length,
                    num_threads=num_threads,
                    hostname=record["environment"]["os"]["node"],
                )
                rows.append(row)
    commands = [row["command"] for row in rows]
    if commands != ["generate_ancestors", "match_ancestors", "match_samples"]:
        raise ValueError(f"Expected all three inference stages in {path}: {commands}")
    return pd.DataFrame(rows)


def aggregate_performance(paths):
    """Combine individual or per-version tables without averaging repetitions."""
    frames = [pd.read_csv(path) for path in paths]
    frame = pd.concat(frames, ignore_index=True)
    return frame.sort_values(["prefix", "sample_size", "version", "rep", "command"])


def run_inference(args, config):
    """Run each stage in a fresh process with its dataset's input settings."""
    dataset = get_dataset(config, args.prefix, args.sample_size)
    installation = pathlib.Path(tsinfer.__file__).resolve().parents[1]
    logger.info("Using tsinfer %s from %s", args.version, installation)
    data = tsinfer.VariantData(
        dataset["zarr_path"],
        ancestral_state=dataset["ancestral_state"],
        sequence_length=dataset.get("sequence_length"),
        site_mask=dataset["site_mask"],
        sample_mask=dataset["sample_mask"],
    )
    if data.num_samples != args.sample_size:
        raise ValueError(
            f"Expected {args.sample_size} haplotypes, found {data.num_samples}"
        )
    logger.info(
        "%s: %s, %s haplotypes, rep %s, %s sites",
        args.command,
        args.prefix,
        data.num_samples,
        args.rep,
        data.num_sites,
    )
    matching = {
        "num_threads": args.threads,
        "recombination_rate": config["recombination_rate"],
        "mismatch_ratio": config["mismatch_ratio"],
        "precision": dataset["precision"],
        "progress_monitor": False,
    }
    if args.command == "generate_ancestors":
        ancestors = tsinfer.generate_ancestors(
            data,
            path=args.output,
            num_threads=args.threads,
            genotype_encoding=config["genotype_encoding"],
            progress_monitor=False,
        )
        ancestors.close()
    elif args.command == "match_ancestors":
        with tsinfer.AncestorData.load(args.input[0]) as ancestors:
            ts = tsinfer.match_ancestors(data, ancestors, **matching)
        ts.dump(args.output)
    elif args.command == "match_samples":
        ancestors_ts = tskit.load(args.input[0])
        ts = tsinfer.match_samples(
            data, ancestors_ts, post_process=config["post_process"], **matching
        )
        ts.dump(args.output)
    logger.info("Completed %s: %s", args.command, args.output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in [
        "generate_ancestors",
        "match_ancestors",
        "match_samples",
        "extract",
        "aggregate",
    ]:
        subparser = commands.add_parser(command)
        subparser.add_argument("--config", required=True, type=pathlib.Path)
        subparser.add_argument("--output", required=True)
        if command != "generate_ancestors":
            input_count = "+" if command == "aggregate" else 1
            subparser.add_argument("--input", required=True, nargs=input_count)
        if command != "aggregate":
            subparser.add_argument("--version", required=True)
            subparser.add_argument("--prefix", required=True)
            subparser.add_argument("--sample-size", required=True, type=int)
            subparser.add_argument("--rep", required=True, type=int)
            subparser.add_argument("--threads", required=True, type=int)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    logging.basicConfig(
        level=config["log_level"],
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    if args.command == "match_samples":
        # Released versions log twice per haplotype at INFO in the worker loop.
        logging.getLogger("tsinfer.inference").setLevel(
            config["sample_matching_log_level"]
        )
    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    if args.command == "extract":
        frame = extract_performance(
            args.input[0],
            args.version,
            args.prefix,
            args.sample_size,
            args.rep,
            args.threads,
        )
    elif args.command == "aggregate":
        frame = aggregate_performance(args.input)
    else:
        run_inference(args, config)
        return
    frame.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
