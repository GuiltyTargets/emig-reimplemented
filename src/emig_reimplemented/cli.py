import click
import logging
import os
import time

from emig_reimplemented.pipeline import rank_targets
from emig_reimplemented.config import get_config
from ppi_network_annotation.pipeline import generate_ppi_network

logger = logging.getLogger(__name__)

@click.command()
@click.option('--config-path',
              prompt='Please enter the path to the config file',
              help='Path to config file.')

def main(config_path: str):
    """Run Emig reimplementation

    :param str config_path: The path to the configuration file.
    """
    logging.basicConfig(level=logging.INFO)

    start = time.time()

    click.secho('getting config', color='cyan')
    cfp = get_config(config_path)

    ppi_path = cfp['paths']['ppi_path']
    assert os.path.exists(ppi_path)

    data_path = cfp['paths']['data_path']
    assert os.path.exists(data_path)

    drug_targets_path = cfp['paths']['drug_targets_path']
    assert os.path.exists(drug_targets_path)

    hippie_min_edge_weight = cfp.getfloat('default', 'interaction_confidence_cutoff')
    current_disease_ids_path = cfp.get('paths', 'disease_ids')
    disease_associations_path = cfp.get('paths', 'disease_associations')

    maximum_adjusted_p_value = cfp.getfloat('default', 'maximum_adjusted_p_value')
    maximum_log2_fold_change = cfp.getfloat('default', 'maximum_log2_fold_change')
    minimum_log2_fold_change = cfp.getfloat('default', 'minimum_log2_fold_change')

    entrez_id_header = cfp['dge']['entrez_id']
    log_fold_change_header = cfp['dge']['log2_fold_change']
    adjusted_p_value_header = cfp['dge']['adjusted_p_value']
    split_char = cfp['dge']['split_character']
    base_mean_header = cfp.get('dge', 'base_mean')

    feature_path = cfp.get('paths', 'feature_path')
    auc_emig_output = cfp['paths']['auc_emig']

    diff_type = cfp['options']['diff_type']

    del cfp

    click.secho('generating PPI network', color='cyan')
    network = generate_ppi_network(
        ppi_graph_path=ppi_path,
        gene_expression_file_path=data_path,
        maximum_adjusted_p_value=maximum_adjusted_p_value,
        maximum_log2_fold_change=maximum_log2_fold_change,
        minimum_log2_fold_change=minimum_log2_fold_change,
        entrez_id_header=entrez_id_header,
        log_fold_change_header=log_fold_change_header,
        adjusted_p_value_header=adjusted_p_value_header,
        base_mean_header=base_mean_header,
        split_char=split_char,
        hippie_min_edge_weight=hippie_min_edge_weight,
        current_disease_ids_path=current_disease_ids_path,
        disease_associations_path=disease_associations_path,
    )

    click.secho('ranking targets', color='cyan')
    rank_targets(
        network=network,
        feature_path=feature_path,
        diff_type=diff_type,
        targets_path=drug_targets_path,
        output_auc_path=auc_emig_output,
    )