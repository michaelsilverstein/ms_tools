"Tools for working with phylogenetic trees"

from skbio.tree import TreeNode
from skbio.diversity.beta import weighted_unifrac, unweighted_unifrac
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from argparse import ArgumentParser, RawTextHelpFormatter

def setup_unifrac(tree_file, features, weighted=True):
    """
    Setup unifrac function to pass to `pairwise_distances`
    :param tree: skbio-readable tree
    :param features: Features to keep in the tree
    :param weighted: True for weighted_unifrac (default), False for unweighted_unifrac
    :return: unifrac funtion
    """
    # Load tree
    print('Loading tree...')
    tree = TreeNode.read(tree_file)

    # Shear to OTUs of interest
    print('Shearing tree...')
    sheared = tree.shear(features).root_at_midpoint()

    # Get metric function
    fnc = weighted_unifrac if weighted else unweighted_unifrac
    def unifrac(u, v):
        return fnc(u, v, features, sheared)
    return unifrac

if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', help='Path to OTU count table', required=True, type=str)
    parser.add_argument('-d', help='Delimiter for OTU table [Default: "\t"]', type=str, default='\t')
    