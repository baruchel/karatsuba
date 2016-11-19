# karatsuba
Fast convolution algorithms with Python types

A module for performing repeated convolutions involving high-level Python objects (which includes large integers, rationals, SymPy terms, Sage objects, etc.). By relying on Karatsuba's algorithm, the function is faster than available ones for such purpose. When working with sequences of 32 terms (or more), the current module gets even significantly faster than the well-known `numpy.convolve` function (as long as high-level Python objects are involved rather than CPU types).

## Background

Fast convolutions are achievable with low-level types by using FFT (Fast Fourier Transform) based algorithms, and very good modules (usually wrappers for C or Fortran code) are available for this task. But as soon as coefficients to be convolved are arbitrary-size integers, rationals, or any high-level types, these low-level functions can't be used; usually the best thing to do is to write two or three lines of code for performing the convolution by using the very basic recursive formula; a Python programmer can also use the `convolve` function from the Numpy module with arrays having `object` as their `dtype`, but in this very specific case it will not be much faster.

It is well-known that Karatsuba's algorithm may have some benefit here, since it requires less multiplications (which are likely to be the slowest operations with such objects), but from a practical point of view, naive implementations of the algorithm are often slower than the direct formula for technical reasons: many recursive calls, inefficient memory allocations, etc.

The strategy followed by the current module is then to compute plans for a given kind of convolution: a first call to the module will pre-compute the convolution with symbols (rather than values) by using Karatsuba's algorithm, build a tree of all arithmetical operations, clean the tree if needed, then return an _ad-hoc_ function for later use. _Of course, the initial call is slow (from a computational point of view), but repeated calls (in a loop for instance) may be twice as fast as any other available function._

Since the module builds plans, some additional features are provided in order to let the user call the plans without having to reshape the input or output lists: the plan will be clever enough to discard unwanted values (and avoid computing them) or to pick the input coefficients at various locations in a list (for instance, if the coefficients to be convolved are not at the beginning of the list, if they are in reversed order, if odd or even coefficients are to be taken, etc.).

## Installation

The module is compatible with all versions of Python (Python2, Python3, Pypy2, Pypy3).

Just type:

    pip install karatsuba

or (for a system-wide installation):

    sudo pip install karatsuba

## Usage

For importing the module, the following line is enough:

    from karatsuba import *

The simplest use is (for performing the convolution of two sequences of 16 terms):

    k = make_plan(range(16), range(16))

    from fractions import Fraction
    l = [Fraction(1,x) for x in range(1,17)]
    print(k(l, l))

The previous plan will return 32 terms, but it is well known that the last term is 0, which is why some functions (like `numpy.convolve`) only return 31 terms; the last term can indeed be discarded:

    k = make_plan(range(16), range(16), plan=[True]*31+[False])

which will avoid (later) the need of removing the last term if this is required.

The rule for customizing the output is very easy: each boolean value in the `plan` optional argument tells if the coefficient of the corresponding degree is wanted or not (of course, disabling a coefficient will both shorten the output list and discard some useless computation).

Sometimes, an algorithm may require extracting the input coefficients from another list; but this may have some computational cost. If this has to be performed repeatedly, it is better to tune the plan: the two initial arguments map any index to the given degree. Thus `[3,2,1,0]` (instead of `range(4)`) tells that the four coefficients should be read in reverse order. Another meaningful example could be: `[4,5,6,7]` in order to use four terms _starting from the index 4 in the input list_. Since two lists are involved, two different mappings may be provided.

If the sequences to be convolved have a size which is not a power of two, the expected result may be achieved by replacing unexisting indices with `None`; for instance for two sequences of 7 terms:

    k = make_plan([0,1,2,3,4,5,6,None], [0,1,2,3,4,5,6,None], [True]*13+[False]*3) 

If the user wants to study (or tune) the returned function instead of using it, the `raw` option can be set to `True`; in this case, the plan will not be returned as a function but as a string (containing Python materials).
    

## Performance

Some tests are provided below. For several kinds of lists, a simple plan (corresponding to the same behaviour than `numpy.convolve`) is built; then the plan is applied repeatedly (by convolving the sequence with itself) a great number of time. The same convolution is performed with `numpy.convolve` and with the following piece of code.

    def convolution(l1, l2):
        N = len(l1)-1
        l = []
        for k in range(2*N+1):
            l.append(sum( l1[i]*l2[k-i] for i in range(max(0, k-N), min(k+1, N+1))))
        return l

The whole computation was performed with Python 3; results are:
