import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *


def run_plugin_local(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    plugin_local can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('plugin_local')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.set_local_option('PLUGIN_LOCAL', 'PRINT', 1)
    psi4.plugin('plugin_local.so')

# Integration with driver routines
procedures['energy']['plugin_local'] = run_plugin_local


def localize():
    psi4.plugin('plugin_local.so')
    pass
