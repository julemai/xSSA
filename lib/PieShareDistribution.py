#!/usr/bin/env python
from __future__ import print_function

import numpy as np

__all__ = ['PieShareDistribution']

def PieShareDistribution(M, N, remainder=True, randomnumbers=None):
    """
        This function produces M sets of N random numbers that sum up to C.
        The N random numbers are distributed identically with the PDF f(x) = N * (1-x)^(N-1).

        Definition
        ----------
        def PieShareDistribution(M, N, C, remainder=True)


        Input           Format        Description
        -----           -----         -----------
        M               integer       Number of random number sets. Each set is on length N.
        N               integer       Number of random numbers that sum up to C.                                         
        remainder       bool          If False, only (N-1) random numbers are returned.
                                      The one that is not returned is the remainder:
                                             r_N = C-r_1-r_2-...-r_(N-1)
                                      If True, remainder is calculated and returned.
                                      Default: True
        randomnumbers   1D/2D array   Sampled uniform distributed random numbers that will be 
                                      used to derive the weights
                                      M x (N-1) random numbers need to be provided to retrieve 
                                      M x N weights
                                      Default: None

        Output          Format        Description
        -----           -----         -----------
        rand            1D/2D array   Matrix (M>1) or Vector (M=1) of random number sets
                                      Dim 1: M = number of random number sets
                                      Dim 2: N = number of random numbers per set
                                                

        Description
        -----------
        We assume that there are N alternatives/options that need to be weighted such that the weights sum up to C.
        The steps to obtain the N weights w_1, w_2, ..., w_N are:
        1) Sample (N-1) independent random numbers r_i with i=1,(N-1) that are uniformly distributed between [0,1].
        2) Derive the first N-1 weights w_i one after the other:
           For i=1,(N-1):
                w_i = ( 1 - sum_{j=1}^{i-1} w_j ) * ( 1 - (1 - r_i)**(1/(N-i)) )
        3) Derive the remainder:
                w_N = ( 1 - sum_{j=1}^{N-1} w_j )
        The probability density function of the N weights are all the same. The PDF is f(x) = N * (1-x)^(N-1).


        Restrictions
        ------------
        None

        Examples
        --------

        Without remainder

        >>> np.random.seed(seed=12345)
        >>> M = 3    # number of random number sets
        >>> N = 5    # number of random numbers that need to sum up to 1.0
        >>> rand = PieShareDistribution(M, N, remainder=True)
        >>> print('random set #1: '+str(rand[0,:]))
        random set #1: [ 0.48492752  0.06133198  0.04384398  0.08384854  0.32604797]
        >>> print('random set #2: '+str(rand[1,:]))
        random set #2: [ 0.18915093  0.21120041  0.4866893   0.07378246  0.03917689]
        >>> print('random set #3: '+str(rand[2,:]))
        random set #3: [ 0.29212136  0.21071729  0.24744715  0.24005194  0.00966226]
        >>> print('sum: '+str(np.sum(rand[:,:],axis=1)))
        sum: [ 1.  1.  1.]

        No remainder

        >>> np.random.seed(seed=12345)
        >>> M = 3    # number of random number sets
        >>> N = 5    # number of random numbers that need to sum up to 1.0
        >>> rand = PieShareDistribution(M, N, remainder=False)
        >>> print('random set #1: '+str(rand[0,:]))
        random set #1: [ 0.48492752  0.06133198  0.04384398  0.08384854]
        >>> print('random set #2: '+str(rand[1,:]))
        random set #2: [ 0.18915093  0.21120041  0.4866893   0.07378246]
        >>> print('random set #3: '+str(rand[2,:]))
        random set #3: [ 0.29212136  0.21071729  0.24744715  0.24005194]

        Provide random numbers

        >>> np.random.seed(seed=12345)
        >>> M = 3    # number of random number sets
        >>> N = 5    # number of random numbers that need to sum up to 1.0
        >>> rr = np.random.rand(M,N-1)
        >>> rand = PieShareDistribution(M, N, remainder=True, randomnumbers=rr)
        
        License
        -------
        This file is part of the the PieShareDistribution Python package.

        The PieShareDistribution Python package is free software: you can redistribute it and/or modify
        it under the terms of the GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        The PieShareDistribution Python package is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with the PieShareDistribution project (cf. gpl.txt and lgpl.txt).
        If not, see <http://www.gnu.org/licenses/>.

        Copyright 2019 Juliane Mai - juliane.mai@uwaterloo.ca


        History
        -------
        Written,  Juliane Mai, June 2019
    """

    if randomnumbers is None:
        # sample uniform random numbers ~ U[0,1]
        rr = np.random.rand(M,N-1)

    else:
        if len(np.shape(randomnumbers)) == 1:
            if np.shape(randomnumbers)[0] != N-1:
                raise ValueError("PieShareDistribution: N-1 random numbers need to be provided" )
            else:
                randomnumbers = randomnumbers[np.newaxis,:]
        elif len(np.shape(randomnumbers)) == 2:
            if ( np.shape(randomnumbers)[0] != M or np.shape(randomnumbers)[1] != N-1):
                raise ValueError("PieShareDistribution: random numbers provided need to have M rows and (N-1) columns" )
        else:
            raise ValueError("PieShareDistribution: random numbers provided need to be either 1D or 2D arrays" )

        if (np.any(randomnumbers > 1.0) or np.any(randomnumbers < 0.0)):
            raise ValueError("PieShareDistribution: random numbers provided need to be uniform distributed between 0 and 1")

        rr = randomnumbers

    # initialize matrix of weights
    ww = np.zeros((M,N))

    # derive the first (N-1) weights
    for ii in range(0,N-1):
        ss = np.sum(ww[:,0:ii],axis=1)
        ww[:,ii] = ( 1. - ss ) * ( 1. - (1. - rr[:,ii])**(1./(N-ii-1)) )

    if remainder:
        # derive remainder
        ss = np.sum(ww[:,0:N-1],axis=1)
        ww[:,N-1] = ( 1. - ss )
    else:
        # remove last column because remainder is not returned
        ww = ww[:,0:N-1]

    # just return a 1D array
    if M == 1:
        return ww[0,:]

    return ww
    


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # M = 3000
    # N = 10

    # weights = PieShareDistribution(M,N)

