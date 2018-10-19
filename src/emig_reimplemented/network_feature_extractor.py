# -*- coding: utf-8 -*-

"""Includes class to calculate the network-based features."""

import logging

import numpy as np
import pandas as pd
from ppi_network_annotation.model.network import Network
from igraph import Vertex
from scipy import sparse
from sklearn.preprocessing import normalize
from tqdm import tqdm

logger = logging.getLogger(__name__)
np.set_printoptions(precision=3)


# TODO: add formula to the documentation

class NetworkFeatureExtractor:
    """Class for calculating features using interaction and differential expression information."""

    def __init__(self, network: Network):
        """Construct the extractor object.

        :param network: The graph using which the features will be calculated.
        :param cfp: An object that includes paths, cutoffs and other necessary information.
        """
        self.ppi_network = network
        self.ppi_network.graph.simplify(combine_edges=min)

    def extract_features(self, feature_path: str, diff_type: str) -> None:
        """Extract all graph-related features and write to a file.

        :param feature_path: Path to write the file.
        :param diff_type: Differential expression type chosen by the user; all, down, or up.
        """
        logger.info("In extract_features()")

        neighborhood_scores = self.score_neighborhood()
        interconnectivity2_scores = self.score_interconnectivity(diff_type, "second-degree")
        random_walk_scores = self.score_by_random_walk(diff_type)
        network_prop_scores = self.score_by_network_propagation(diff_type)

        df = pd.DataFrame({
            "GeneID": self.ppi_network.graph.vs["name"],
            "Neighborhood": neighborhood_scores,
            "Interconnectivity": interconnectivity2_scores,
            "RandomWalk": random_walk_scores,
            "NetworkProp": network_prop_scores,
        })

        logger.info('Writing network to %s', feature_path)
        df.to_csv(feature_path,
                  encoding="utf-8",
                  sep="\t",
                  index=False)

    def score_by_local_radiality(self, diff_type: str = "all"):
        """Score by local radiality algorithm.

        :param str diff_type: Differential expression type chosen by the user; all, down, or up.
        :return list: A list of scores, sorted by node index.
        """
        logger.info("In score_by_local_radiality()")
        diff_expr_genes = self.ppi_network._get_differentially_expressed_genes(diff_type)
        scores = [self._local_radiality(gene, diff_expr_genes) for gene in
                  self.ppi_network.graph.vs]
        return scores

    def _local_radiality(self, gene, diff_expr_genes):
        """Calculate local radiality for one gene.

        :param Vertex gene: Gene whose score will be calculated.
        :param list diff_expr_genes: List of differentially expressed genes
        :return: Local radiality score for gene.
        """
        shortest_paths = self.ppi_network.graph.get_shortest_paths(gene,
                                                                   to=diff_expr_genes)
        path_length_sum = sum([len(path) for path in shortest_paths])
        return path_length_sum / len(diff_expr_genes)

    def score_neighborhood(self) -> list:
        """Score all nodes using neighborhood scoring algorithm.

        :return list: A list of scores, sorted by node index.
        """
        logger.info("In neighborhood_scoring()")
        return list(map(self._neighborhood, self.ppi_network.graph.vs))

    def _neighborhood(self, node: Vertex) -> float:
        """Score a node based on its and its neighbours' log fold change.

        :param Vertex node: Node to be scored.
        :return float: Score of the node.
        """
        node_fc = abs(node["log2_fold_change"])
        sum_fc = 0
        for n in node.neighbors():
            sum_fc += abs(n["log2_fold_change"])
        if len(node.neighbors()) > 0:
            return 0.5 * node_fc + 0.5 * sum_fc / len(node.neighbors())
        else:
            return 0

    def score_interconnectivity(self, diff_type: str = "all",
                                neighbor_type: str = "direct") -> list:
        """Score all nodes based on interconnectivity algorithm.

        :param str diff_type: Differential expression type chosen by the user; all, down, or up.
        :return list: A list of scores, sorted by node index.
        """
        logger.info("In interconnectivity_nodes()")
        icn_mat = self._interconnectivity_edges(diff_type, neighbor_type)
        diff_expr = self.ppi_network._get_differentially_expressed_genes(diff_type)

        icn = np.sum(icn_mat[diff_expr.indices, :], axis=0) / len(diff_expr)
        return list(icn)

    def _interconnectivity_edges(self, diff_type: str = "all",
                                 neighbor_type: str = "direct") -> np.ndarray:
        """Score pairs of nodes based on their shared neighborhood.

        :param str diff_type: Differential expression type chosen by the user; all, down, or up.
        :return np.ndarray: A matrix of scores for pairs.
        """
        key = self._get_diff_expr_key(diff_type)
        nodes = list(self.ppi_network.graph.vs)
        degrees = self.ppi_network.graph.degree(nodes)
        icn_mat = np.zeros([len(nodes), len(nodes)], dtype=float)

        diff_expressed = self._get_diff_expr_vertices(diff_type).indices
        edges = self.ppi_network.graph.es.select(_within=diff_expressed)

        for edge in tqdm(edges, desc="Interconnectivity"):
            icn_score, source, target = self._interconnectivity_edge(degrees, edge, key,
                                                                     neighbor_type)
            icn_mat[source.index, target.index] = icn_score
            icn_mat[target.index, source.index] = icn_score
        return icn_mat

    def _interconnectivity_edge(self, degrees, edge, key, neighbor_type):
        source = self.ppi_network.graph.vs.find(edge.source)
        target = self.ppi_network.graph.vs.find(edge.target)
        icn_score = 0
        if edge != -1 and (source[key] or target[key]):
            inter = self.ppi_network.get_neighborhood_overlap(source, target, neighbor_type)
            mult_degrees = degrees[source.index] * degrees[target.index]
            if mult_degrees > 0:
                icn_score = (2 + len(inter)) / np.sqrt(mult_degrees)

        return icn_score, source, target

    def score_by_random_walk(self, diff_type: str = "all") -> list:
        """Score nodes using random walk algorithm (Koehler et al).

        :param str diff_type: Differential expression type chosen by the user; all, down, or up.
        :return list: List of scores, sorted by node index.
        """
        logger.info("In random_walk()")
        self._random_walk_init(diff_type)

        adj = sparse.coo_matrix(np.array(self.ppi_network.graph.get_adjacency().data,
                                         dtype="float64"))
        adj = normalize(adj, norm="l1", axis=0)  # column normalized
        return self._walk_randomly(adj, "random_walk_score", 0.5)

    def _random_walk_init(self, diff_type: str = "all") -> None:
        """Initialize the graph for random walk algorithm.

         By setting attribute "random_walk_score" to 1/no_of_diff_expressed
         for differentially expressed genes.
        :param str diff_type: Differential expression type chosen by the user; all, down, or up.
        """
        self.ppi_network.graph.vs["random_walk_score"] = 0

        if diff_type == "up":
            prob = 1 / len(self.ppi_network.graph.vs.select(up_regulated_eq=True))
            self.ppi_network.graph.vs.select(up_regulated_eq=True)["random_walk_score"] = prob
        elif diff_type == "down":
            prob = 1 / len(self.ppi_network.graph.vs.select(down_regulated_eq=True))
            self.ppi_network.graph.vs.select(down_regulated_eq=True)["random_walk_score"] = prob
        else:
            prob = 1 / len(self.ppi_network.graph.vs.select(is_diff_expressed_eq=True))
            self.ppi_network.graph.vs.select(is_diff_expressed_eq=True)["random_walk_score"] = prob

    def score_by_diffusion(self):
        """Score nodes using heat diffusion algorithm.

        :return list: A list of scores, sorted by node index.
        """
        logger.info("In score_by_diffusion()")
        self.ppi_network.graph.vs["diffusion_score"] = np.absolute(
            self.ppi_network.graph.vs["log2_fold_change"])

        adj = sparse.coo_matrix(np.array(self.ppi_network.graph.get_adjacency().data,
                                         dtype="float64"))
        adj = normalize(adj, norm="l1", axis=0)  # column normalized
        return self._walk_randomly(adj, "diffusion_score", 0.5)

    def score_by_network_propagation(self, diff_type: str = "all") -> list:
        """Score nodes using network propagation algorithm.

        :param str diff_type: Differential expression type chosen by the user; all, down, or up.
        :return list: A list of scores, sorted by node index.
        """
        logger.info("In propagate_network()")
        self._propagate_network_init(diff_type)

        adj = sparse.dok_matrix(np.array(self.ppi_network.graph.get_adjacency().data,
                                         dtype="float64"))
        # normalized by the degrees of source and target nodes
        adj = self._normalize_by_degrees(adj)
        return self._walk_randomly(adj, "network_prop_score", 0.5)

    def _propagate_network_init(self, diff_type: str = "all") -> None:
        """Initialize the graph for network propagation algorithm.

        By setting attribute "network_prop_score" to 1 for differentially
        expressed genes.
        :param str diff_type: Differential expression type chosen by the user; all, down, or up.
        """
        self.ppi_network.graph.vs["network_prop_score"] = 0
        vertices = self.ppi_network.graph.vs

        if diff_type == "up":
            vertices.select(up_regulated_eq=True)["network_prop_score"] = 1
        elif diff_type == "down":
            vertices.select(down_regulated_eq=True)["network_prop_score"] = 1
        else:
            vertices.select(is_diff_expressed_eq=True)["network_prop_score"] = 1

    def _normalize_by_degrees(self,
                              adj: sparse.dok_matrix) -> sparse.dok_matrix:
        """Normalize an adjacency matrix based on the node degrees(Vanunu et al).

        :param adj: Adjacency matrix to be normalized.
        :return: Normalized adjacency matrix.
        """
        row_sums = np.sum(adj, axis=0)
        dia = row_sums + 1
        norm_adj = sparse.dok_matrix(np.zeros(adj.shape))
        for key in adj.keys():
            norm_adj[key] = adj[key] / np.sqrt(dia[0, key[0]] * dia[0, key[1]])
        return norm_adj

    def _walk_randomly(self, adj, score_type, alpha):
        """ Randomly walk on the network while updating the visitation probabilities.

        :param adj: Normalized adjacency matrix.
        :param score_type: One of random_walk_score, diffusion_score, or network_prop_score.
        :param alpha: Probability of restarting the walk.
        :return: Vector of updated visitation probabilities.
        """
        # initialize for first round
        p0 = self.ppi_network.graph.vs[score_type]
        pt1 = p0
        pt2 = self._update_visitation_probabilities(p0, pt1, adj, alpha)
        while self._l1_norm(pt1, pt2) > 10 ** -6:
            pt1 = pt2
            pt2 = self._update_visitation_probabilities(p0, pt1, adj, alpha)
        return pt2

    def _update_visitation_probabilities(self, p0, p1, adj, alpha):
        """Update the visitation probabilities.

        :param p0: scores at time point 0.
        :param p1: scores at time point t.
        :param alpha: Weighting factor.
        :return: p2: scores at time point t+1.
        """
        p1 = np.array(p1, dtype="float64")
        p0 = np.array(p0, dtype="float64")
        p2 = (1 - alpha) * adj.dot(p1) + alpha * p0
        return p2

    def _l1_norm(self, v1, v2) -> float:
        """Calculate the L1 norm of two vectors.

        :param v1: Vector 1.
        :param v2: Vector 2.
        :return: L1 norm of v1 and v2.
        """
        return sum(
            abs(a - b)
            for a, b in zip(v1, v2)
        )


    def _get_diff_expr_vertices(self, diff_type):
        if diff_type == "up":
            return self.ppi_network.graph.vs.select(up_regulated_eq=True)

        if diff_type == "down":
            return self.ppi_network.graph.vs.select(down_regulated_eq=True)

        return self.ppi_network.graph.vs.select(is_diff_expressed_eq=True)

    def _get_diff_expr_key(self, diff_type):
        if diff_type == "up":
            return "up_regulated"

        if diff_type == "down":
            return "down_regulated"

        return "is_diff_expressed"
