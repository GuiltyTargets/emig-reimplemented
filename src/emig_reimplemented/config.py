# -*- coding: utf-8 -*-

"""Configuration for gene-prioritization."""

import os
from configparser import ConfigParser

__all__ = [
    'get_config',
]


def get_config(path: str) -> ConfigParser:
    """Read the configuration file and set up necessary parameters.

    :param path: Path to the configuration file.
    :return: A ConfigParser object with necessary fields filled in
    """
    assert os.path.exists(path), f'path does not exist: {path}'

    cfp = ConfigParser()
    cfp.read(path)

    # Read input file paths and join them with the input directory
    input_directory = cfp['paths']['input_directory']
    ppi_path = os.path.join(input_directory, cfp.get('paths', "protein_protein_interaction_graph"))
    data_path = os.path.join(input_directory, cfp.get('paths', 'differential_gene_expression'))
    targets_path = os.path.join(input_directory, cfp.get('paths', "drug_targets"))
    cfp.set('paths', 'ppi_path', ppi_path)
    cfp.set('paths', 'data_path', data_path)
    cfp.set('paths', 'targets_path', targets_path)

    # Read the differential expression type
    diff_type = cfp.get('options', 'diff_type')

    # Read output file paths and join them with the output directory
    output_directory = cfp.get('paths', 'output_directory')
    feature_path = os.path.join(output_directory, f"NetworkAnalysis-{diff_type}.tsv")
    auc_path = os.path.join(output_directory, "auc_emig.tsv")
    cfp.set('paths', 'feature_path', feature_path)
    cfp.set('paths', 'auc_path', auc_path)

    return cfp
