#!/usr/bin/env python3

"""Pipeline for emig-reimplementation.

This can be run as a module with ``python -m gene_prioritization.cli``
or through the installed command ``gene_prioritization.
"""

import logging
import os
import time
# import rpy2.robjects as ro
import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn import linear_model

from ppi_network_annotation.model.network import Network
from emig_reimplemented.node_scoring import NodeScorer

__all__ = [
    'rank_targets'
]

logger = logging.getLogger(__name__)
HERE = os.path.abspath(os.path.dirname(__file__))


# def check_r_installation():
#     """Check if R and necessary R packages are installed."""
#     try:
#         home_dir = os.path.join(HERE, "GENHelper", "R", "RLoader.R")
#         ro.r("source('{}')".format(home_dir))
#     except Exception:
#         raise Exception("R or necessary R packages could not be loaded. "
#                         "Please check your installation")

def evaluate_cv_nested(features, labels, n_splits):
    """Do a repeated stratified cross validation.

    :param clf: Classifier object.
    :param features: The feature matrix.
    :param n_splits: Number of folds.
    :return: Dictionary containing numerical results of the classification.
    """
    results = defaultdict(list)
    grid = {
        'C': np.logspace(-4, 4, 20),
    }
    log_reg = linear_model.LogisticRegression(solver='lbfgs')

    # tol, C
    for i in range(10):
        inner_cv = StratifiedKFold(n_splits=n_splits, shuffle=True)
        outer_cv = StratifiedKFold(n_splits=n_splits, shuffle=True)

        for train_idx, test_idx in outer_cv.split(features, labels):
            clf = GridSearchCV(estimator=log_reg, param_grid=grid, cv=inner_cv, iid=False)
            clf.fit(features, labels)

            X_train, X_test, Y_train, Y_test = get_split(features, labels, test_idx, train_idx)
            pred, probs = clf.predict(X_test), clf.predict_proba(X_test)
            results["auc"].append(roc_auc_score(Y_test, probs[:, 1]))

    return results


def evaluate_cv(features, labels, n_splits):
    """Do a repeated stratified cross validation.

    :param clf: Classifier object.
    :param features: The feature matrix.
    :param n_splits: Number of folds.
    :return: Dictionary containing numerical results of the classification.
    """
    results = defaultdict(list)
    for i in range(10):
        outer_cv = StratifiedKFold(n_splits=n_splits, shuffle=True)

        for train_idx, test_idx in outer_cv.split(features, labels):
            X_train, X_test, Y_train, Y_test = get_split(features, labels, test_idx, train_idx)

            clf = linear_model.LogisticRegression(solver='lbfgs')
            clf.fit(X_train, Y_train)

            pred, probs = clf.predict(X_test), clf.predict_proba(X_test)
            results["auc"].append(roc_auc_score(Y_test, probs[:, 1]))

    return results


def get_split(features, labels, test_id, train_id):
    return features[train_id], features[test_id], labels[train_id], labels[test_id]


def rank_targets(
        network: Network,
        diff_type: str,
        targets_path: str,
        # feature_path: str,
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

    # # Check if necessary R packages are installed
    # check_r_installation()

    # os.makedirs(os.path.dirname(feature_path), exist_ok=True)
    os.makedirs(os.path.dirname(auc_path), exist_ok=True)

    # Quit if no differentially expressed genes were found
    if len(network.get_upregulated_genes()) + len(network.get_downregulated_genes()) == 0:
        raise Exception("No differentially expressed genes were found. "
                        "Please check your input and the parameters.")

    # Score nodes using different network-based measures
    # logger.info('Calculating node scores to write to %s', feature_path)
    feature_extractor = NodeScorer(network)
    features = feature_extractor.score_nodes(diff_type)
    features['GeneID'] = features['GeneID'].astype(int)

    targets = pd.read_csv(targets_path, sep='\t', names=['GeneID', 'Label'])
    targets['GeneID'] = targets['GeneID'].astype(int)

    print(features)
    print(targets)

    data_set = features.merge(targets, how='inner', on='GeneID')

    print(len(data_set))

    features = data_set[["Neighborhood", "Interconnectivity", "RandomWalk", "NetworkProp"]]
    labels = data_set["Label"]

    results = evaluate_cv_nested(features, labels, 5)
    with open(auc_path) as f:
        for auc in results:
            f.write("{}\n".format(auc))

    features = data_set[["LocalRadiality"]]
    labels = data_set["Label"]

    results = evaluate_cv(features, labels, 5)
    with open(auc_path.replace('emig', 'isik')) as f:
        for auc in results:
            f.write("{}\n".format(auc))

    logger.info(f"Target Prioritization - Elapsed time in minutes: {(time.time() - start) / 60}")
