"Consumer resource model utilities"

import numpy as np
import pandas as pd

def TrophicResourceMatrix(resource_class_sizes, transition_matrix, sparsity=.2):
    """
    Generate the resource matrix, D, with provided transition probabilites between resource classes
    Inputs:
    | resource_class_sizes: List of sizes for each resource class
    | transition_matrix: Transition probabilities from each resource to the others
    | sparsity: Controls sparsity of matrix {0, 1}
    """
    n_classes = len(resource_class_sizes)
    n_resources = sum(resource_class_sizes)
    transition_matrix = np.array(transition_matrix)
    assert transition_matrix.shape == (n_classes, n_classes), 'Transitions must be provided to and from all resource classes.'
    
    # Initialize D matrix (originally in X out)
    resource_names = ['R%d'%r for r in range(n_resources)]
    type_names = [x for i, n in enumerate(resource_class_sizes) for x in ['T%d' % i] * n]
    resource_index = pd.MultiIndex.from_arrays((type_names, resource_names), names=['resource_type', 'resource'])
    D = pd.DataFrame(np.zeros((n_resources, n_resources)), resource_index, resource_index)
    resource_types = D.index.get_level_values('resource_type').unique()
    
    # Sample D matrix for each resource type
    for resource_type, transitions, resource_class_size in zip(resource_types, transition_matrix, resource_class_sizes):
        # Generate transition probabilities based on transitions prob and class size
        alpha = np.array([x for t, s in zip(transitions, resource_class_sizes) for x in [t / s] * s])
        # Sample from dirichlet after imposing sparsity
        Db = np.random.dirichlet(alpha / sparsity, resource_class_size)

        # Add to D
        D.loc[resource_type] = Db
    
    # Transpose to out X in to match community simulator input
    return D.T

def TrophicConsumerMatrix(consumer_class_sizes, consumer_class_preferences, resource_class_sizes, sparsity=None, n=None):
    """
    Generate consumer matrix, c, given resource type preferences for each consumer family
    
    Inputs:
    | consumer_class_sizes: List of size of each consumer family
    | consumer_class_preferences: Resource preferences of each family
    | sparsity: Controls sparsity of matrix {0, 1}
    | n: Number of resource preferences for each consumer
    
    ** ONLY `sparsity` or `n` can be provided
    """
    # Check that only sparsity or n has been provided
    sampling_mode = 'sparsity' if sparsity is not None else 'n'
    if (sparsity is not None) & (n is not None):
        raise ValueError('Only `sparsity` or `n` can be provided, not both.')
    
    # Consumer features
    n_consumers = sum(consumer_class_sizes)
    n_families = len(consumer_class_sizes)
    
    # Resource features
    consumer_names = ['S%d' % i for i in range(n_consumers)]
    family_names = [x for i, n in enumerate(consumer_class_sizes) for x in ['F%d' % i] * n]
    families = ['F%d'%i for i in range(n_families)]
    consumer_index = pd.MultiIndex.from_arrays([family_names, consumer_names])

    "Resource parameters"
    n_resource_classes = len(resource_class_sizes)
    n_resources = sum(resource_class_sizes)
    resource_names = ['R%d'%r for r in range(n_resources)]
    type_names = [x for i, n in enumerate(resource_class_sizes) for x in ['T%d' % i] * n]
    resource_index = pd.MultiIndex.from_arrays((type_names, resource_names))

    assert np.array(consumer_class_preferences).shape == (n_families, n_resource_classes), 'Consumer preferences must be provided for each resource class.'

    # Initialize c matrix
    c = pd.DataFrame(np.zeros((n_consumers, n_resources)), consumer_index, resource_index)

    for species_class_size, consumer_preferences, family in zip(consumer_class_sizes, consumer_class_preferences, families):
        # Mask which points to choose
        consumer_class_mask = [x for p, s in zip(consumer_preferences, resource_class_sizes) for x in [p] * s]

        if sampling_mode == 'sparsity':
            # Sample preferences from this families preference set and impose sparsity
            c_class = ((np.random.rand(species_class_size, n_resources) * consumer_class_mask) > sparsity).astype(int)
        
        else:
            # Get indices of class
            class_idx = np.where(consumer_class_mask)[0]
            # Sample idxs
            c_class_idx = c_class_idx = np.array([np.random.choice(class_idx, n, replace=False) for _ in range(species_class_size)])
            c_class = np.zeros((species_class_size, n_resources))
            for s, idx in zip(c_class, c_class_idx):
                s[idx] = 1
            
        # Add to C matrix
        c.loc[family] = c_class
    
    return c.rename_axis(index=['family', 'species'], columns=['resource_type', 'resource'])