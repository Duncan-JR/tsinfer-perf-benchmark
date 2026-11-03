"""Extract ancestor spans between their first and last inference sites."""

import argparse
import logging
import pathlib
import re

import numpy as np
import pandas as pd
import tsinfer

logger = logging.getLogger(__name__)


def extract_spans(output_dir):
    """Return focal-site spans for every dataset, version and repetition.

    Read ``ancestors/<prefix>-n<count>-rep<rep>-v<version>-ancestors.zarr``
    using the benchmark's legacy tsinfer runtime (see README.md). These are
    Zarr stores in LMDB files. Only ancestor metadata is loaded.

    Column names follow tsinfer-anc-eval, with prefix, sample_size and rep
    identifying the benchmark run. ``inf_focal_site`` indexes the ancestor
    store's inference sites; ``inferred_node`` is its ancestor ID. For the
    half-open interval [start, end), ``inferred_span`` is the distance in bp
    between positions[end - 1] and positions[start], while
    ``inferred_site_span`` counts included inference sites. A one-site
    ancestor has span zero. Sequence length is not used in this calculation.
    Each focal site contributes one row; ancestors without focal sites,
    including the artificial roots, contribute none.
    """
    output_dir = pathlib.Path(output_dir)
    paths = sorted(output_dir.glob("ancestors/*-ancestors.zarr"))
    if len(paths) == 0:
        raise FileNotFoundError(f"No generated ancestor stores in {output_dir}")
    frames = []
    for path in paths:
        match = re.fullmatch(
            r"(?P<prefix>.+)-n(?P<sample_size>\d+)-rep(?P<rep>\d+)-v(?P<version>.+)-ancestors\.zarr",
            path.name,
        )
        if match is None:
            raise ValueError(f"Invalid ancestor store filename: {path}")
        version = match["version"]
        prefix = match["prefix"]
        sample_size = int(match["sample_size"])
        rep = int(match["rep"])
        logger.info(
            "Reading tsinfer %s, %s, n=%s, rep=%s", version, prefix, sample_size, rep
        )
        with tsinfer.AncestorData.load(str(path)) as ancestors:
            positions = ancestors.sites_position[:]
            start = ancestors.ancestors_start[:]
            end = ancestors.ancestors_end[:]
            times = ancestors.ancestors_time[:]
            focal_sites = ancestors.ancestors_focal_sites[:]
            counts = np.array([len(sites) for sites in focal_sites], dtype=int)
            ancestor_ids = np.arange(ancestors.num_ancestors)
            inferred_node = np.repeat(ancestor_ids, counts)
            inf_focal_site = np.concatenate(focal_sites).astype(int)
            last_site = end - 1
            spans = positions[last_site] - positions[start]
            site_spans = end - start
            frame = pd.DataFrame(
                {
                    "version": version,
                    "prefix": prefix,
                    "sample_size": sample_size,
                    "rep": rep,
                    "num_samples": sample_size,
                    "inferred_node": inferred_node,
                    "inf_focal_site": inf_focal_site,
                    "focal_position": positions[inf_focal_site],
                    "inferred_time": times[inferred_node],
                    "inferred_span": spans[inferred_node],
                    "inferred_site_span": site_spans[inferred_node],
                }
            )
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=pathlib.Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    frame = extract_spans(args.output_dir)
    path = args.output_dir / "dataframes" / "ancestor_spans.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    logger.info("Wrote %s focal-site rows to %s", len(frame), path)


if __name__ == "__main__":
    main()
