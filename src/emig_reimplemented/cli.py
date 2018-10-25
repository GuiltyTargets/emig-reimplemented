import click
import logging
import os

from emig_reimplemented.pipeline import rank_targets
from emig_reimplemented.config import get_config
from ppi_network_annotation.pipeline import generate_ppi_network, parse_dge

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

    click.secho('getting config', color='cyan')
    cfp = get_config(config_path)

    # Read input file paths and assert they exist
    ppi_graph_path = cfp['paths']['ppi_path']
    assert os.path.exists(ppi_graph_path)
    dge_path = cfp['paths']['data_path']
    assert os.path.exists(dge_path)
    targets_path = cfp['paths']['targets_path']
    assert os.path.exists(targets_path)

    # Read the thresholds
    min_confidence = cfp.getfloat('default', 'interaction_confidence_cutoff')
    max_padj = cfp.getfloat('default', 'maximum_adjusted_p_value')
    max_l2fc = cfp.getfloat('default', 'maximum_log2_fold_change')
    min_l2fc = cfp.getfloat('default', 'minimum_log2_fold_change')

    # Read the headers of the differential expression file
    entrez_header = cfp['dge']['entrez_id']
    l2fc_header = cfp['dge']['log2_fold_change']
    adjp_header = cfp['dge']['adjusted_p_value']
    base_mean_header = cfp.get('dge', 'base_mean')
    entrez_delimiter = cfp['dge']['split_character']

    # Read the output file paths
    feature_path = cfp.get('paths', 'feature_path')
    auc_path = cfp['paths']['auc_path']

    # Read the differential expression type (up, down, all)
    diff_type = cfp['options']['diff_type']

    click.secho('generating PPI network', color='cyan')
    gene_list = parse_dge(
        dge_path=dge_path,
        entrez_id_header=entrez_header,
        log2_fold_change_header=l2fc_header,
        adj_p_header=adjp_header,
        entrez_delimiter=entrez_delimiter,
        base_mean_header=base_mean_header,
    )
    network = generate_ppi_network(
        ppi_graph_path=ppi_graph_path,
        dge_list=gene_list,
        max_adj_p=max_padj,
        max_log2_fold_change=max_l2fc,
        min_log2_fold_change=min_l2fc,
        ppi_edge_min_confidence=min_confidence,
    )

    click.secho('ranking targets', color='cyan')
    rank_targets(
        network=network,
        diff_type=diff_type,
        targets_path=targets_path,
        feature_path=feature_path,
        auc_path=auc_path,
    )
