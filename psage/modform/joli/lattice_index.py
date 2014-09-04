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
from sage.rings.integer_ring import IntegerRing
from sage.matrix.constructor import matrix


def LatticeIndex( x):
    """
    Return an instance of Lattice_class initialized by $x$.

    INPUT
        x can be:

        o a list of integers representing the elements
          of a Gram matrix above the diagonal read from left
          to right and top to bottom
        o a Gram matrix
        o a string representing an integral lattice:
          A_n, B_n, C_n, D_n, E_6, E_7, E_8, G_2, F_4, Z^n
    """

    # x is a list of the entries of a Gram matrix above the diagonal
    try:
        return Lattice_class( x)
    except:
        pass

    # x is a Gram matrix
    try:
        n = x.nrows()
        assert x == x.transpose()
        xl = [x[i,j] for j in range(n) for i in range(j,n)]
        return Lattice_class( xl)
    except:
        pass

    # x is a string describing a lattice

    try:
        # A_n: o-o-...-o
        assert 'A_' == x[:2]
        n = Integer( x[2:])
        assert n >= 1
        G = matrix( IntegerRing(), n, n, 2)
        for i in range(n-1):
            G[i,i+1] = G[i+1,i] = -1
        return LatticeIndex( G)
    except:
        pass

    try:
        # B_n: o-o-...-o=>o
        assert 'B_' == x[:2]
        n = Integer( x[2:])
        assert n >= 2
        G = matrix( IntegerRing(), n, n, 2)
        G[n-1,n-1] = 1
        for i in range(n-1):
            G[i,i+1] = G[i+1,i] = -1
        return LatticeIndex( G)
    except:
        pass

    try:
        # C_n: o-o-...-o<=o
        assert 'C_' == x[:2]
        n = Integer( x[2:])
        assert n >= 2
        G = matrix( IntegerRing(), n, n, 2)
        G[n-1,n-1] = 4
        for i in range(n-2):
            G[i,i+1] = G[i+1,i] = -1
        G[n-1,n-2] = G[n-2,n-1] = -2
        return LatticeIndex( G)
    except:
        pass

    try:
        # D_n: o-o-...-o=8
        assert 'D_' == x[:2]
        n = Integer( x[2:])
        assert n >= 3
        G = matrix( IntegerRing(), n, n, 2)
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
        G = matrix( IntegerRing(), n, n, 2)
        for i in range(n-2):
            G[i,i+1] = G[i+1,i] = -1
        G[n-1,2] = G[2,n-1] = -1
        return LatticeIndex( G)
    except:
        pass

    try:
        # F_4: o-o=>o-o
        assert 'F_4' == x
        n = 4
        G = matrix( IntegerRing(), n, n, 4)
        G[2,2] = G[3,3] = 2
        G[0,1] = G[1,0] = -2
        G[1,2] = G[2,1] = -2
        G[2,3] = G[3,2] = -1
        return LatticeIndex( G)
    except:
        pass

    try:
        # G_2: o<3o
        assert 'G_2' == x
        G = matrix( IntegerRing(), 2, 2, [2,-3,-3,6])
        return LatticeIndex( G)
    except:
        pass

    try:
        # Z^n
        assert 'Z^' == x[:2]
        n = Integer( x[2:]) 
        return LatticeIndex( matrix( IntegerRing(), n, n, 1))
    except:
        pass

    # try:
    #     # R^#(h): dual of a lattice spanned by the root system R
    #     # rescaled by the Coxeter number h of R
    #     y = x[:-5]
    #     assert '^#(h)' == y
    #     L = LatticeIndex( y)
    #     G = L.gram_matrix()^-1
    #     h = ???
    #     return LatticeIndex( h*G)
    # except:
    #     pass


    raise NotImplementedError( '%s: not recognized as representing a lattice'%x)
