"""
A module for computing fast convolutions with Python types by precomputing
Karatsuba's algorithm for a given case.
"""

__version__ = '0.1.1'

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
    from .py3 import compile_plan
else:
    from .py2 import compile_plan

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
    # Internal functions
    class _atom:
        def __init__(self, op):
            self.is_zero = False
            self.children = []
            self.parents = []
            self.ref = None
            self.content = None
            self.op = op
            self.top = None
    def _add(l1, l2):
        r = []
        for a, b in zip(l1, l2):
            if b.is_zero: # both a=0 and b=0 or only b=0
                r.append(a)
            elif a.is_zero:
                r.append(b)
            else:
                c = _atom("add")
                c.children = [a, b]
                a.parents.append(c)
                b.parents.append(c)
                r.append(c)
        return r
    def _sub(l1, l2):
        r = []
        for a, b in zip(l1, l2):
            if b.is_zero: # both a=0 and b=0 or only b=0
                r.append(a)
            elif a.is_zero:
                c = _atom("neg")
                c.children = [b]
                b.parents.append(c)
                r.append(c)
            else:
                c = _atom("sub")
                c.children = [a, b]
                a.parents.append(c)
                b.parents.append(c)
                r.append(c)
        return r
    def _karatsuba(l1, l2):
        if len(l1) == 1:
            if l1[0].is_zero or l2[0].is_zero:
                x = _atom("0")
                x.is_zero = True
                x.content = "0"
            else:
                x = _atom("mul")
                x.children = [l1[0], l2[0]]
                l1[0].parents.append(x)
                l2[0].parents.append(x)
            y = _atom("0")
            y.is_zero = True
            y.content = "0"
            return [x, y]
        n = len(l1)
        m = n // 2
        b, a = l1[:m], l1[m:]
        d, c = l2[:m], l2[m:]
        X = _karatsuba(a, c)
        Y = _karatsuba(b, d)
        Z = _karatsuba(_sub(a, b), _sub(c, d))
        T = _sub(_add(X, Y), Z)
        return Y[:m] + _add( Y[m:] + X[:m], T) + X[m:]
    def _parse_tree(t, parse):
        if (t.op != "origin" or t.top != None) and not t.is_zero:
            if id(t) not in map(id, parse):
                parse.append(t)
        for c in t.children: _parse_tree(c, parse)
    # Main function
    try:
        l1 = [None if x==None else int(x) for x in l1]
        l2 = [None if x==None else int(x) for x in l2]
    except:
        raise TypeError
    n1, n2 = len(l1), len(l2)
    if n1 != n2:
        raise ValueError("Input lists have different lengths")
    if plan == None:
        plan = [True]*(2*n1)
    elif len(plan) != (2*n1):
        raise ValueError("Optional 'plan' parameter should be twice"
                        + " as long as input lists")
    if n1 ^ (1<<(n1).bit_length()>>1) != 0:
        raise ValueError("Length of input lists should be a power of 2")
    x, y = [], []
    for k in range(len(l1)):
        if l1[k] != None:
            a = _atom("origin")
            a.is_zero = False
            a.content = "u["+str(l1[k])+"]"
            a.ref = a.content
            x.append(a)
        else:
            a = _atom("0")
            a.is_zero = True
            a.content = "0"
            x.append(a)
        if l2[k] != None:
            a = _atom("origin")
            a.is_zero = False
            a.content = "v["+str(l2[k])+"]"
            a.ref = a.content
            y.append(a)
        else:
            a = _atom("0")
            a.is_zero = True
            a.content = "0"
            y.append(a)
    n, parse = 0, []
    for i, k in enumerate(_karatsuba(x, y)):
        if plan[i]:
            k.top = n
            _parse_tree(k, parse)
            n += 1
    n, out = 0, []
    stats = {} if stats == None else stats
    stats['mul'], stats['add'], stats['sub'], stats['neg'] = 0, 0, 0, 0
    while len(parse) > 0:
        parse.sort(key=lambda x: len([c for c in x.children if c.ref==None]))
        i = 0
        while len([c for c in parse[i].children if c.ref==None]) == 0:
            if parse[i].op == "mul":
                parse[i].content = "*".join(
                        [c.ref for c in parse[i].children])
                stats['mul'] += 1
            elif parse[i].op == "add":
                parse[i].content = "+".join(
                        [c.ref for c in parse[i].children])
                stats['add'] += 1
            elif parse[i].op == "sub":
                parse[i].content = (parse[i].children[0].ref + "-"
                                    + parse[i].children[1].ref)
                stats['sub'] += 1
            elif parse[i].op == "neg":
                parse[i].content = ("-" + parse[i].children[0].ref)
                stats['neg'] += 1
            elif parse[i].op == "origin":
                parse[i].content = parse[i].content
            if parse[i].top == None:
                l = "t.append(" + parse[i].content + ")"
                parse[i].ref = "t["+str(n)+"]"
                l += " #  -->  " + parse[i].ref
                out.append(l)
                n += 1
            else:
                l = "r["+str(parse[i].top)+"]="+parse[i].content
                parse[i].ref = "r["+str(parse[i].top)+"]"
                out.append(l)
            i += 1
            if i == len(parse): break
        parse = parse[i:]
    d =  "def %s(u, v):\n"
    d += "    t, r = [], [0]*"+str(sum(plan))+"\n"
    d += "    " + "\n    ".join(out) + "\n"
    d += "    return r    # Stats: "+str(stats)
    if raw: return (d%"convolution")
    return compile_plan("convolution", d)

__all__ = ['make_plan']
