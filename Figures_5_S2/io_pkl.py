import collections
import dill as pickle
import fnmatch
import os


def load_pkl(file):
    """Loads data from pickle

    Parameters
    ----------
    file : str

    Returns
    -------
    data
    """
    f = open(file, 'rb')
    data = pickle.load(f)
    f.close()
    return data


def save_pkl(data, file):
    """Saves data to a pickle

    Parameters
    ----------
    data : object
    file : str
    """
    f = open(file, 'wb')
    pickle.dump(data, f)
    f.close()
