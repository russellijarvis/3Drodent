#    This file is part of pyEntropy
#
#    pyEntropy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    pyEntropy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with pyEntropy. If not, see <http://www.gnu.org/licenses/>.
#
#    Copyright 2009, 2010 Robin Ince
import numpy as np
from nose.tools import assert_raises
from numpy.testing import *
#from pyentropy.utils import *

def setup():
    global x, x2, y, a1, b1, a2, b2

    x = np.arange(3**3)
    x2 = np.atleast_2d(x).T
    y = np.array([[0, 0, 0, 0],
       [0, 0, 0, 1],
       [0, 0, 0, 2],
       [0, 0, 1, 0],
       [0, 0, 1, 1],
       [0, 0, 1, 2],
       [0, 0, 2, 0],
       [0, 0, 2, 1],
       [0, 0, 2, 2],
       [0, 1, 0, 0],
       [0, 1, 0, 1],
       [0, 1, 0, 2],
       [0, 1, 1, 0],
       [0, 1, 1, 1],
       [0, 1, 1, 2],
       [0, 1, 2, 0],
       [0, 1, 2, 1],
       [0, 1, 2, 2],
       [0, 2, 0, 0],
       [0, 2, 0, 1],
       [0, 2, 0, 2],
       [0, 2, 1, 0],
       [0, 2, 1, 1],
       [0, 2, 1, 2],
       [0, 2, 2, 0],
       [0, 2, 2, 1],
       [0, 2, 2, 2]])

    a1 = np.array([8, 9, 7, 9, 3, 3, 9, 7, 9, 2])
                         #8 occurs 1/10, 9 occurs 4/10 3 occurs 2/10 7 occurs 2/10
    b1 = np.array([0, 0, 1, 2, 0, 0, 0, 2, 1, 4]) / 10.0
    a2 = np.array([0, 1, 7, 1, 3, 3, 1, 7, 1, 2])
    b2 = np.array([1, 4, 1, 2, 0, 0, 0, 2, 0, 0]) / 10.0


def test_decimalise():
    assert_equal(decimalise(y.T,4,3),x)
def test_prob_naive():
    assert_equal(prob(a1,10), b1)
"""
Count number of occurrences of each value in array of non-negative ints.

The number of bins (of size 1) is one larger than the largest value in x. If minlength is specified, there will be at least this number of bins in the output array (though it will be longer if necessary, depending on the contents of x). Each bin gives the number of occurrences of its index value in x. If weights is specified the input array is weighted by it, i.e. if a value n is found at position i, out[n] += weight[i] instead of out[n] += 1.
Parameters :	

For example
In [35]: a=np.bincount(a1)

In [36]: a
Out[36]: array([0, 0, 1, 2, 0, 0, 0, 2, 1, 4])
                0 occurs 0 times, 1 occurs 0 times, 2 occurs 1  3 occurs 2, 4,5,6occcur zero times.
In [37]: a1
Out[37]: array([8, 9, 7, 9, 3, 3, 9, 7, 9, 2])


x : array_like, 1 dimension, nonnegative ints

    Input array.

weights : array_like, optional

    Weights, array of the same shape as x.

minlength : int, optional

    New in version 1.6.0.

    A minimum number of bins for the output array.

Returns :	

out : ndarray of ints

    The result of binning the input array. The length of out is equal to np.amax(x)+1.

Raises :	

ValueError

    If the input is not 1-dimensional, or contains elements with negative values, or if minlength is non-positive.

TypeError

    If the type of the input is float or complex.


"""
setup()

pr_vs=np.divide(np.bincount(a1.astype(float)),len(a1.astype(float)))

entrop=sum(pr_vs[:]*np.log2(pr_vs[:]))
print entrop
def test_prob_naive_missed_responses():
    assert_equal(prob(a2,10), b2)    
    
    
def teardown():
    global x, x2, y, a1, b1, a2, b2
    del x, x2, y, a1, b1, a2, b2
    
def test_dec2base_1d():
    assert_equal(dec2base(x,3,4),y)

def test_dec2base_2d():
    assert_equal(dec2base(x2,3,4),y)

def test_dec2base_noncol():
    assert_raises(ValueError, dec2base, x2.T, 3, 4)
    
def test_base2dec():
    assert_equal(base2dec(y,3),x)

def test_decimalise_error():
    assert_raises(ValueError, decimalise, y, 3, 4)


    
def test_pt_bayescount():
    # values match original bayescount.m file
    for n,r in [(100000, 5.0), (50, 5.0), (30, 6.0),
                (12, 7.0), (10, 8.0), (8, 9.0), (7, 10.0)]:
        yield check_pt_bayes, n, r
        
def check_pt_bayes(n, r):
    assert_equal(pt_bayescount(b1,n),r)
    
if __name__ == '__main__':
    run_module_suite()
