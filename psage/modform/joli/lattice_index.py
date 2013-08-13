#*****************************************************************************
#       Copyright (C) 2013 Nils-Peter Skoruppa <nils.skoruppa@gmail.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
CREATES LATTICE INDEX FOR JACOBI FORMS
Provides classes and various functions to do computations
in the context of Jacobi forms of lattice index over $Q$.
"""
_VERSION = '$Id$'

from lattice import Lattice_class
from sage.rings.integer import Integer
from sage.matrix.constructor import matrix


def LatticeIndex( x):
    """
    Return an nstance of Lattice_class initialized by $x$.

    INPUT
        x can be:

        o a list of integers representing the elements
          of a Gram matrix above the diagonal read from left
          to right and top to bottom
        o a Gram matrix
        o a string representing an integral lattice
    """
    
    # x is a list of the entries of a Gram matrix above the diagonal
    try:
        return Lattice_class( x)
    except:
        pass

    # x is a Gram matrix
    try:
        n = x.nrows()
        xl = [x[i,j] for i in range(n) for j in range(i,n)]
        return Lattice_class( xl)
    except:
        pass

    # x is a string describing a lattice

    try:
        # A_n: o-o-...-o
        assert 'A_' == x[:2]
        n = Integer( x[2:])
        assert n >= 1
        G = matrix( ZZ, n, n, 2)
        for i in range(n-1):
            G[i,i+1] = G[i+1,i] = -1
        return LatticeIndex( G)
    except:
        pass

    try:
        # D_n: o-o-...-o=8
        assert 'D_' == x[:2]
        n = Integer( x[2:])
        assert n >= 3
        G = matrix( ZZ, n, n, 2)
        for i in range(n-2):
            G[i,i+1] = G[i+1,i] = -1
        G[n-1,n-3] = G[n-3,n-1] = -1
        return LatticeIndex( G)
    except:
        pass

    try:
        # E_n:     o
        #          |
        #      o-o-o-o-...-o
        assert 'E_' == x[:2]
        n = Integer( x[2:])
        assert n in [6,7,8]
        G = matrix( ZZ, n, n, 2)
        for i in range(n-2):
            G[i,i+1] = G[i+1,i] = -1
        G[n-1,2] = G[2,n-1] = -1
        return LatticeIndex( G)
    except:
        pass

    try:
        # Z^n
        assert 'Z^' == x[:2]
        n = Integer( x[2:]) 
        return LatticeIndex( matrix( ZZ, n, n, 1))
    except:
        pass

    raise NotImplementedError( '%s: not recognized as representing a lattice'%x)
