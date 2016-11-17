"""
A module for computing fast convolutions with Python types by precomputing
Karatsuba's algorithm for a given case.
"""

__version__ = '0.1'

import sys

def _convolution(l1, l2):
    """
    Classic formula implemented for testing purposes only.
    """
    N = len(l1)-1
    l = []
    for k in range(2*N+2):
        l.append(sum( l1[i]*l2[k-i] for i in range(max(0, k-N), min(k+1, N+1))))
    return l

if sys.version_info[0] >= 3:
    from .py3 import make_plan as __mkp
else:
    from .py2 import make_plan as __mkp

def make_plan(l1, l2, plan=None, raw=False, stats = None):
    """
    Return a plan for performing a convolution of two sequences of objects
    (numbers, symbols, etc.) of a given length.

    The plan is a function of two arguments (lists); the lists may have a
    different size (most likely larger) than l1 and l2 as long as the mapping
    from the plan remains consistent.

    Parameters
    ----------
    l1, l2 : list like (generally using range(n) will be enough)
             Lists of integers; the length of both lists must be the same;
             furthermore, the length must be a power of 2.
             Zero-padding can be done by replacing unneeded coefficients with
             None.
             The integer numbers are the index of the corresponding coefficient
             in the lists to be convolved (later).

    plan :   list of booleans; the length of the list must be twice the length
             of l1 or l2. Each boolean is True if the corresponding coefficient
             has to be computed and returned or False if it is not needed (and
             not returned).
    
    raw :    mostly for study purpose; instead of returning a function, the code
             of the computed function is returned (as a string); it may also be
             used for converting it into another language.
    
    stats :  a dictionary in which some statistics will be inserted (number of
             multiplications, additions, etc.

    Returns
    -------
    out : function
          A function (of two arguments) for performing the convolution.
    """
    stats = {} if stats == None else stats
    return __mkp(l1, l2, plan=plan, raw=raw, stats=stats)

__all__ = ['make_plan']
