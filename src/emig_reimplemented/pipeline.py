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
import emig_reimplemented.network_feature_extractor as nfe

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


def rank_targets(network: Network,
                 feature_path: str,
                 diff_type: str,
                 targets_path: str,
                 output_auc_path: str
                 ) -> None:
    """Rank genes based on a given set of special genes.

    :param Network network: An instance of Network class.
    :param Params cfp: A Params object containing info on paths, constant fields, and options.
    :param str config_path: The path to the configuration file.
    """
    start = time.time()

    # Check if necessary R packages are installed
    check_r_installation()

    os.makedirs(os.path.dirname(feature_path), exist_ok=True)
    os.makedirs(os.path.dirname(output_auc_path), exist_ok=True)

    # Quit if no differentially expressed genes were found
    if len(network.get_upregulated_genes()) + len(network.get_downregulated_genes()) == 0:
        raise Exception("No differentially expressed genes were found. "
                        "Please check your input and the parameters in .ini file")

    if not feature_path or not os.path.exists(feature_path):
        logger.info('extracting features to %s', feature_path)
        feature_extractor = nfe.NetworkFeatureExtractor(network)
        feature_extractor.extract_features(feature_path, diff_type)

    logger.info("Prioritizing...")

    r_home = os.path.join(HERE, "GENHelper", "R")
    ro.r("targets.path <- '" + targets_path + "'")
    ro.r("features.path <- '" + feature_path + "'")
    ro.r("auc.output.path <- '" + output_auc_path + "'")

    ro.r("source('{}')".format(
        os.path.join(r_home, "TargetPrioritization.R")))

    logger.info("Target Prioritization - Elapsed time in minutes: {}".
                format((time.time() - start) / 60))

# /home/omuslu/Documents/gene-prioritization/data/config.ini

