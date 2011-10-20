""" Utility functions """

import inspect
import matplotlib.pyplot as plt
import numpy as np


__suffix = ''


def set_suffix(suffix):
    global __suffix
    __suffix = suffix

def whosparent():
    """Return parent function name of caller"""

    return inspect.stack()[2][3]

def saveplot(suffix=''):
    """Save a plot using caller's name"""

    if suffix:
        suffix = '-' + suffix
    plt.savefig('plots/%s%s%s.pdf' % (whosparent(), suffix, __suffix))

def title(text):
    plt.title(text + '\n(%s)' % __suffix)

mylog = np.vectorize(lambda x: np.log10(x) if x > 0 else 0)
