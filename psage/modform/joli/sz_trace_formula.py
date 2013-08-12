#*****************************************************************************
#       Copyright (C) 2007 Nils-Peter Skoruppa <nils.skoruppa@gmail.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
TRACES OF HECKE OPERATORS ON THE SPACE OF NEWFORMS ON $\Gamma_0(m)$

DESCRIPTION
    This is an implementation of the trace formulas $s_{k,m}(l,n)$
    as given in the article [S-Z].

REFERENCES
    [S-Z] N-P. Skoruppa and D. Zagier,
          Jacobi forms and a certain space of modular forms,
          Inv. Math. 94 (1988), 113--146
AUTHORS
    -- Nils-Peter Skoruppa (2007-08-15)
"""
_VERSION = '$Id$'


from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.interfaces.gp import gp
from sage.rings.arith import sigma
from sage.misc.functional import is_odd


def _gegenbauer_pol( h, s):
    r"""
    OUTPUT
        The coefficient of $x^{2*h}$ in the Taylor expansion
        of  $1/(1 - \sqrt{s} x + x^2)$  around  $0$.

    INPUT
        h -- a positive integer
        s -- element of a ring
    """
    h = Integer( h)
    if h < 0:
        raise ValueError, "first argument %s must be a nonnegative integer" % h
    if 0 == h:
        return s**0
    if 1 == h:
        return s - s**0
    i = 2; a = s - s**0; b = s**0
    while i <= h:
        c = b; b = a; a = (s-2)*b - c; i += 1
    return a


def _ab( m):
    r"""
    OUTPUT
        The unique pair of integers $a,b$ such that $m=a^2b$
        with squarefree $b$.

    INPUT
        m -- a nonzero integer
    """
    m = Integer(m)
    b = m.squarefree_part()
    a = (m//b).isqrt()
    return a, b


def hurwitz_kronecker_class_no_x( n, D):
    r"""
    OUTPUT
        The generalized Hurwitz-Kronecker class number $H_n(D)$
        as defined in [S-Z].

    INPUT
        n -- an integer $\ge 1$
        D -- an integer $\le 0$
    """
    n = Integer(n)
    D = Integer(D)
    if n <= 0:
        raise ValueError, "%s: must be an integer >=1" % n
    if D > 0:
        raise ValueError, "%s: must be an integer <= 0" % D
    g = D.gcd( n)
    a,b = _ab(g)
    h = g*b
    if 0 != D%h:
        return 0
    Dp = D//h
    H1 = Rational( gp.qfbhclassno( -Dp))
    if 0 == H1:
        return H1
    return g * Dp.kronecker( n//g) * H1


def _alpha( n):
    r"""
    OUTPUT
        Value at $n$ of the multiplicative arithmetic function
        which takes values $-1$ for $p$, $p^2$ and $+1$ for $p^3$,
        and $0$ otherwise, where $p$ is a prime, cf. [S-Z].

    INPUT
        n -- an integer $\ge 1$
    """
    return sum( moebius(n//d//d)*moebius(d) for d in _ab(n)[0].divisors())


def _sz_s( k, m, l, n):
    r"""
    OUTPUT
        The function  $s_{k,m}(l,n)$, i.e.~the trace
        of $T(l) \circ W_n$ on the "certain" space
        $\mathcal{M}_{2k-2}^{\text{cusp}}(m)$ of modular forms
        as defined in [S-Z].

    INPUT
        l -- index of the Hecke operator, rel. prime to the level $m$
        n -- index of the Atkin-Lehner involution (exact divisor of $m$)
        k -- integer ($2k-2$ is the weight)
        m -- level
    """
    k = Integer(k)
    m = Integer(m)
    l = Integer(l)
    n = Integer(n)
    if k < 2:
        return 0
    if 1 != m.gcd(l):
        raise ValueError, \
              "gcd(%s,%s) != 1: not yet implemented" % (m,l)

    ellT = 0
    for np in n.divisors():
        ellT -= sum( l**(k-2) \
                     * _gegenbauer_pol( k-2, s*s*np/l) \
                     * hurwitz_kronecker_class_no_x( m/n, s*s*np*np - 4*l*np) \
                     for s in range( -(Integer(4*l*np).isqrt()//np), \
                                     1+(Integer(4*l*np).isqrt()//np)) \
                     if (n//np).gcd(s*s).is_squarefree())/2

    parT = -sum( min( lp, l//lp)**(2*k-3) \
                 * _ab(n)[0].gcd( lp+l//lp) \
                 * _ab(m/n)[0].gcd( lp-l//lp) \
                 for lp in l.divisors())/2
    
    if 2 == k and (m//n).is_square():
        corT = sigma(n,0)*sigma(l,1)
    else:
        corT = 0
        
    return Integer( ellT + parT + corT)


def __check( k, m, l, n):
    k = Integer(k)
    m = Integer(m)
    l = Integer(l)
    n = Integer(n)

    if m <= 0:
        raise ValueError, "%s: m must be a positive integer" % m

    if l < 1 or 1 != l.gcd( m):
        raise ValueError, "%s: l must be a positive integer prime to %d" \
              % ( l, m)

    if n <=0 or 0 != m%n or 1 != n.gcd( m//n):
        raise ValueError, "%s: n must be an exact divisor of %d" \
              ( n, m)

    return k, m, l, n

    
####################################################################
# Exported Functions
####################################################################


def trace_special_cusp_forms( k, m, l, n):
    r"""
    OUTPUT
        The function  $s_{k,m}(l,n)$, i.e.~the trace
        of $T(l) \circ W_n$ on the "certain" space
        $\mathcal{M}_{2k-2}^{\text{cusp}}(m)$ of modular forms
        as defined in [S-Z].

    INPUT
        l -- index of the Hecke operator, rel. prime to the level $m$
        n -- index of the Atkin-Lehner involution (exact divisor of $m$)
        k -- integer ($2k-2$ is the weight)
        m -- level
    """
    k, m, l, n = __check( k, m, l, n)
    if is_odd(k) or k <= 0:
        return 0
    return _sz_s( k//2+1, m, l, n)


def trace_cusp_forms( k, m, l, n):
    r"""
    OUTPUT
        The trace of  $T(l) \circ W_n$  on $S_k(Gamma_0(m)$).

    INPUT
        l -- Hecke operator index (rel. prime to the level $m$)
        n -- Atkin-Lehner involution index (exact divisor of the level $m$)
        k -- weight
        m -- a level ($\ge 1$)
    """
    k, m, l, n = __check( k, m, l, n)
    if is_odd(k) or k <= 0:
        return 0
    return sum( _sz_s( k//2+1, mp, l, n.gcd( mp)) \
                    for mp in m.divisors() if (n//mp).is_squarefree())


def trace_new_cusp_forms( k, m, l, n):
    r"""
    OUTPUT
        Rhe trace of  $T(l) \circ W_n$  on $S_k^{\text{new}}(Gamma_0(m)$).

    INPUT
        l -- Hecke operator index (rel. prime to the level $m$)
        n -- Atkin-Lehner involution index (exact divisor of the level $m$)
        k -- weight
        m -- a level ($\ge 1$)
    """
    k, m, l, n = __check( k, m, l, n)
    if is_odd(k) or k <= 0:
        return 0
    return sum( _alpha( m//mp) * _sz_s( k//2+1, mp, l, n.gcd(mp)) \
                for mp in m.divisors())


# def trace_eisenstein_series( l, n, k, m):
#     r"""
#     OUTPUT
#         The trace of  $T(l) \circ W_n$  on $M_k^{\text{Eis.}}(Gamma_0(m)$).

#     INPUT
#         l -- Hecke operator index (rel. prime to the level $m$)
#         n -- Atkin-Lehner involution index (exact divisor of the level $m$)
#         k -- weight
#         m -- a level ($\ge 1$)
#     """
#     k, m, l, n = __check( k, m, l, n)
#     if is_odd(k) or k < 0:
#         return 0
#     if 0 == k and 1 != n:
#         return 0
#     if 0 == k and 1 == n:
#         return sigma(l,0)
#     raise NotImplementedError, "requested traces not yet implemented"


# def trace_new_eisenstein_series( l, n, k, m):
#     r"""
#     OUTPUT
#         The trace of  $T(l) \circ W_n$  on $M_k^{\text{Eis., new}}(Gamma_0(m)$).

#     INPUT
#         l -- Hecke operator index (rel. prime to the level $m$)
#         n -- Atkin-Lehner involution index (exact divisor of the level $m$)
#         k -- weight
#         m -- a level ($\ge 1$)
#     """
#     k, m, l, n = __check( k, m, l, n)
#     if is_odd(k) or k < 0:
#         return 0
#     if 0 == k:
#         if 1 == n == m:
#             return sigma(l,0)
#         return 0
#     if False == is_square(m) or (2 == k and 1 == m):
#         return 0
#     raise NotImplementedError, "requested traces not yet implemented"    


# def trace_modular_forms( l, n, k, m):
#     r"""
#     OUTPUT
#         The trace of  $T(l) \circ W_n$  on $M_k(Gamma_0(m)$).

#     INPUT
#         l -- Hecke operator index (rel. prime to the level $m$)
#         n -- Atkin-Lehner involution index (exact divisor of the level $m$)
#         k -- weight
#         m -- a level ($\ge 1$)
#     """
#     k, m, l, n = __check( k, m, l, n)
#     return trace_eisenstein_series( l, n, k, m) \
#            + trace_cusp_forms( l, n, k, m)


# def trace_new_modular_forms( l, n, k, m):
#     r"""
#     OUTPUT
#         The trace of  $T(l) \circ W_n$  on $M_k^{\text{new}}(Gamma_0(m)$).

#     INPUT
#         l -- Hecke operator index (rel. prime to the level $m$)
#         n -- Atkin-Lehner involution index (exact divisor of the level $m$)
#         k -- weight
#         m -- a level ($\ge 1$)
#     """
#     k, m, l, n = __check( k, m, l, n)
#     return trace_new_eisenstein_series( l, n, k, m) \
#            + trace_new_cusp_forms( l, n, k, m)
