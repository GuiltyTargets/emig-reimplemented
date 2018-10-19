"""Module to test network_feature_extractor module."""
import unittest

from ppi_network_annotation.model.network import Network
from emig_reimplemented import network_feature_extractor as feat
from igraph import Graph
import numpy as np
from sklearn.preprocessing import normalize


class NetworkFeatureExtractorTest(unittest.TestCase):
    """Class that tests the network_feature_extractor module."""

    def __init__(self, *args, **kwargs):
        """Initialize the class by setting up a dummy graph.

        :param args:
        :param kwargs:
        """
        super(NetworkFeatureExtractorTest, self).__init__(*args, **kwargs)
        self.tg = Graph()
        self.tg.add_vertices(7)
        self.tg.add_edges([(0, 1), (0, 3), (1, 2), (1, 3), (3, 4), (4, 5)])
        self.tg.vs["log2_fold_change"] = [2, 1, -1, -2, 1, 2, 1]
        self.tg.vs["is_diff_expressed"] = [True, False, False, True, False,
                                           True, False]
        self.tg = Network(
            self.tg,
            maximum_adjusted_p_value=0.05,
            maximum_log2_fold_change=-1,
            minimum_log2_fold_change=1
        )
        # print(self.tg)

    def test_neighborhood(self):
        """Test neighborhood scoring method."""
        print("+test_neighborhood")
        extractor = feat.NetworkFeatureExtractor(self.tg)
        self.assertAlmostEqual(extractor.score_neighborhood(),
                               [1.75, 1.3333333333333335, 1, 1.6666666666666665,
                                1.5, 1.5, 0])
        self.assertAlmostEqual(extractor._neighborhood(self.tg.graph.vs[0]), 1.75)
        self.assertAlmostEqual(extractor._neighborhood(self.tg.graph.vs[1]),
                               1.3333333333333335)
        self.assertAlmostEqual(extractor._neighborhood(self.tg.graph.vs[2]), 1)
        self.assertAlmostEqual(extractor._neighborhood(self.tg.graph.vs[3]),
                               1.6666666666666665)
        self.assertAlmostEqual(extractor._neighborhood(self.tg.graph.vs[4]), 1.5)
        self.assertAlmostEqual(extractor._neighborhood(self.tg.graph.vs[4]), 1.5)

    def test_interconnectivity(self):
        """Test interconnectivity method."""
        print("+test_interconnectivity")
        mat = np.zeros([len(self.tg.graph.vs), len(self.tg.graph.vs)], dtype=float)
        mat[0, 1] = 1.224744871
        mat[1, 0] = 1.224744871
        mat[0, 3] = 1.224744871
        mat[3, 0] = 1.224744871
        mat[1, 3] = 1.0
        mat[3, 1] = 1.0
        mat[3, 4] = 0.816496581
        mat[4, 3] = 0.816496581
        mat[4, 5] = 1.414213562
        mat[5, 4] = 1.414213562

        extractor = feat.NetworkFeatureExtractor(self.tg)
        icn_mat = extractor._interconnectivity_edges()

        np.testing.assert_array_almost_equal(icn_mat, mat)

        self.assertEqual(extractor.score_interconnectivity(),
                         [0.40824829046386307, 0.7415816237971965, 0.0,
                          0.40824829046386307, 0.74357004776694036, 0.0,
                          0.0])

    def test_random_walk(self):
        """Test random walk method."""
        print("+test_random_walk")
        extractor = feat.NetworkFeatureExtractor(self.tg)
        extractor._random_walk_init()
        self.assertEqual(self.tg.graph.vs["random_walk_score"],
                         [0.3333333333333333, 0, 0, 0.3333333333333333, 0,
                          0.3333333333333333, 0])

        self.assertEqual(extractor._l1_norm([0, 1, 2], [2, 3, 4]), 6)

        p0 = self.tg.graph.vs["random_walk_score"]
        adj = self.tg.graph.get_adjacency()
        adj = np.array(adj.data, dtype="float64")
        adj = normalize(adj, norm="l1", axis=0)
        self.assertEqual(extractor._update_visitation_probabilities(p0, p0, adj,
                                                                    0.5).tolist(),
                         [0.2222222222222222, 0.1388888888888889, 0., 0.25,
                          0.2222222222222222, 0.16666666666666666,
                          0.0])

        scores = extractor.score_by_random_walk()
        self.assertEqual(scores.tolist(),
                         [0.23260918606716072, 0.11454908138395957,
                          0.019091493191259246, 0.28110591916884364,
                          0.14878217937527216, 0.20386214081350457,
                          0.0])

    def test_network_propagation(self):
        """Test network propagation method."""
        print("+test_network_propagation")
        extractor = feat.NetworkFeatureExtractor(self.tg)
        extractor._propagate_network_init()
        self.assertEqual(self.tg.graph.vs["network_prop_score"],
                         [1, 0, 0, 1, 0, 1, 0])
        self.assertEqual(
            extractor.score_by_network_propagation().tolist(),
            [0.6176211255734259, 0.17463663716790856,
             0.030871704149187506, 0.6402662646984134,
             0.20293206698067687, 0.5414233502248208, 0.0])
