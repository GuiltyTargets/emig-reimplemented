# -*- coding: utf-8 -*-

"""Configuration for gene-prioritization."""

import os
from configparser import ConfigParser

__all__ = [
    'get_config',
]


def get_config(path) -> ConfigParser:
    assert os.path.exists(path), f'path does not exist: {path}'

    cfp = ConfigParser()
    cfp.read(path)

    input_directory = cfp['paths']['input_directory']
    ppi_path = os.path.join(
        input_directory,
        cfp.get('paths', "protein_protein_interaction_graph")
    )
    cfp.set('paths', 'ppi_path', ppi_path)
    # set values for differentiating differentially expressed genes

    data_path = os.path.join(
        input_directory,
        cfp.get('paths', 'differential_gene_expression')
    )
    cfp.set('paths', 'data_path', data_path)

    drug_targets_path = os.path.join(input_directory, cfp.get('paths', "drug_targets"))
    cfp.set('paths', 'drug_targets_path', drug_targets_path)

    output_directory = cfp.get('paths', 'output_directory')

    diff_type = cfp.get('options', 'diff_type')
    feature_path = os.path.join(output_directory, f"NetworkAnalysis-{diff_type}.tsv")
    cfp.set('paths', 'feature_path', feature_path)

    cfp.set('paths', 'auc_emig', os.path.join(output_directory, "auc_emig.tsv"))

    # TODO
    # self.DISEASE_ASSOCIATIONS_PATH = os.path.join(self.INPUT_DIR, paths["disease_associations"])
    # self.CURRENT_DISEASE_IDS_PATH = os.path.join(self.INPUT_DIR, paths["disease_ids"])

    return cfp
