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

def make_plan(l1, l2, plan=None, raw=False, stats = {}):
    """
    DOC
    """
    return __mkp(l1, l2, plan=plan, raw=raw, stats=stats)





l = [1,2,3,4,5,6,7,8]
p = [11,12,13,14,15,16,17,18]
k1 = make_plan(range(8), range(8))
k2 = make_plan(range(8), range(8), plan=[True]*15+[False])
print(_convolution(l,p))
print(k1(l,p))
print(k2(l,p))

print("\nExemple pour make_plan(range(4), range(4), raw=True) :")
print(make_plan(range(4), range(4), raw=True))
