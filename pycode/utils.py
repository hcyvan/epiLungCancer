import pandas as pd


def match(a, b):
    """
    get the index of a in b
    :param a:
    :param b:
    :return:
    """
    return pd.Index(pd.Series(b)).get_indexer(a)
