"Tools for working with phylogenetic trees"

from skbio.tree import TreeNode
from skbio.diversity.beta import weighted_unifrac, unweighted_unifrac

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