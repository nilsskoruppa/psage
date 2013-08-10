#*****************************************************************************
#       Copyright (C) 2013 Nils-Peter Skoruppa <nils.skoruppa@gmail.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
Provides classes and various functions to do computations
in the context of Jacobi forms of lattice index over $Q$.

Copyright (C) Nils-Peter Skoruppa 2013
"""
_VERSION = '$Id$'


load '/home/nils/Sandbox/Devel/fqm-devel/cn_group/finite_quadratic_module.sage'
load 'lattice_data.sage'
load 'lattice.sage'


def dim_Jac( k, L, h):
    """
    Return the dimension of the space $J_{k,L}(\varepsilon^h)$.

    INPUT
       k -- a half integer
       L -- an instance of Lattice_class
       h -- an integer
    """
    g = k - h/2
    try:
        g = Integer(g)
    except:
        return Integer(0)

    # preparations
    o_inv = L.o_invariant()
    V = L.shadow_vectors()
    V2 = L.shadow_vectors_of_order_2()
    n = L.rank()

    if k < n/2:
        return Integer(0)
    if n/2 <= k and k < 2+n/2:
        raise NotImplementedError( 'weight %d: not yet implemented'%k)

    # scalar term
    scal = (L.det() + len(V2)*(-1)^(g+o_inv))*(k-n/2-1)/24

    # elliptic S-term
    ell_S = (L.chi(2)*QQbar.zeta(4)^g).real()/4 + kronecker(12,2*g+2*n+1)/6
    
    # elliptic R-term
    ell_R = (L.chi(-3)*QQbar.zeta(24)^(n+2)*QQbar.zeta(6)^g).real()*(-1)^g/(3*QQbar(sqrt(3)))
    
    # parabolic term
    B = lambda t: t-floor(t)-1/2
    par = -sum(B(h/24-L.beta(x)) for x in V)/2
    par -= (-1)^(g+o_inv) * sum(B(h/24-L.beta(x)) for x in V2)/2
    
    return Integer(scal + ell_S + ell_R + par)
