"Utility functions"

def check_len_n(obj, attr, n):
    "Check that `attr` from `obj` has length of n"
    val = getattr(obj, attr)
    if len(val) != n:
        raise ValueError(f'{attr} must be of length {n}')