"ML helper functions"

from sklearn.model_selection import cross_val_predict, LeaveOneOut
import scipy.stats as sps

def loo_score(estimator, X, y, return_predictions=False, **kwargs):
    """
    Compute Leave-One-Out regression score by predicting each value with LOO and
    then calculating R^2 with observations
    Inputs:
    | estimator: scikit-learn estimator
    | X: feature data
    | y: target data
    | return_predictions: Whether to return just score (Default: False) or to include
        predicted y values as well
    | All other keyword arguments are passed to `cross_val_predict()`
    """
    loo_y = cross_val_predict(estimator, X, y, cv=LeaveOneOut(), **kwargs)
    r2 = sps.pearsonr(y, loo_y)[0] ** 2
    if return_predictions:
        return_this = r2, loo_y
    else:
        return_this = r2
    return return_this

def training_score(estimator, X, y, **kwargs):
    """
    Compute training error
    Inputs:
    | estimator: scikit-learn estimator
    | X: feature data
    | y: target data
    | All other keyword arguments are passed to `cross_val_predict()`
    """
    return estimator.fit(X, y, **kwargs).score(X, y)