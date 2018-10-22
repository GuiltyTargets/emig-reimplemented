#!/usr/bin/env python3

"""Pipeline for emig-reimplementation.

This can be run as a module with ``python -m gene_prioritization.cli``
or through the installed command ``gene_prioritization.
"""

import logging
import os
import time
import rpy2.robjects as ro

from ppi_network_annotation.model.network import Network
from emig_reimplemented.node_scoring import NodeScorer

__all__ = [
    'rank_targets'
]
logger = logging.getLogger(__name__)
HERE = os.path.abspath(os.path.dirname(__file__))


def check_r_installation():
    """Check if R and necessary R packages are installed."""
    try:
        home_dir = os.path.join(HERE, "GENHelper", "R", "RLoader.R")
        ro.r("source('{}')".format(home_dir))
    except Exception:
        raise Exception("R or necessary R packages could not be loaded. "
                        "Please check your installation")


def rank_targets(
        network: Network,
        diff_type: str,
        targets_path: str,
        feature_path: str,
        auc_path: str,
) -> None:
    """Rank proteins based on a given set known target proteins.

    :param network: A PPI network annotated with differential gene expression
    :param str diff_type: Differential expression type chosen by the user; all, down, or up.
    :param targets_path: The path to the file that includes Entrez ids of known targets.
    :param feature_path: The path to the file which will include calculated network features
    :param auc_path: The path to the file to which AUC values of CV will be written.
    """
    start = time.time()

    # Check if necessary R packages are installed
    check_r_installation()

    os.makedirs(os.path.dirname(feature_path), exist_ok=True)
    os.makedirs(os.path.dirname(auc_path), exist_ok=True)

    # Quit if no differentially expressed genes were found
    if len(network.get_upregulated_genes()) + len(network.get_downregulated_genes()) == 0:
        raise Exception("No differentially expressed genes were found. "
                        "Please check your input and the parameters.")

    # Score nodes using different network-based measures
    logger.info('Calculating node scores to write to %s', feature_path)
    feature_extractor = NodeScorer(network)
    feature_extractor.score_nodes(feature_path, diff_type)

    # Prioritize candidates using R
    logger.info("Prioritizing...")
    r_home = os.path.join(HERE, "GENHelper", "R")
    ro.r("targets.path <- '" + targets_path + "'")
    ro.r("features.path <- '" + feature_path + "'")
    ro.r("auc.output.path <- '" + auc_path + "'")
    ro.r("source('{}')".format(os.path.join(r_home, "TargetPrioritization.R")))

    logger.info(f"Target Prioritization - Elapsed time in minutes: {(time.time() - start) / 60}")
