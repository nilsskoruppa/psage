#*****************************************************************************
#       Copyright (C) 2007 Nils-Peter Skoruppa <nils.skoruppa@gmail.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
TRACES OF HECKE OPERATORS ON SPACES OF JACOBI FORMS OF INTEGRAL SCALAR INDEX

DESCRIPTION
    This is an implementation of the trace formulas
    as given in the article [S-Z] for holomorphic Jacobi
    forms and as generalized to skew-holomorphic ones in [S-Zh].

REFERENCES
    [S-Z] N-P. Skoruppa and D. Zagier,
          Jacobi forms and a certain space of modular forms,
          Inv. Math. 94 (1988), 113--146
   [S-Zh] N-P. Skorupp and Hai-Gang Zhou,
          Jacobi forms and a certain space of modular forms II,
          to appear
          
AUTHORS
    -- Nils-Peter Skoruppa (2007-08-15)

TODOS
    o traces on Eisenstein series
    o traces on the full space of Jacobi forms
"""
_VERSION = '$Id$'


from sage.rings.integer import Integer
from sz_trace_formula import (_sz_s, _ab, __check)


def __check_jac_params( k, m, kind):
    k = Integer(k)
    m = Integer(m)

    if m <= 0:
        raise ValueError, "index %s: must be a positive integer" % m

    if kind != 'holomorphic' and kind != 'skew_holomorphic':
        raise ValueError, "kind %s: must be %s or %s" \
              % ( kind, 'holomorphic', 'skew_holomorphic')

    return k, m, Integer(-1) if 'holomorphic' == kind else Integer(+1)
 


def trace_jacobi_cusp_forms( k, m, l, n, kind = 'holomorphic'):
    r"""
    OUTPUT
        Trace of  $T(l) \circ W_n$  on the space
        of holomorphic or skew-holomorphic Jacobi cusp forms
        of weight $k$, index $m$ on $\mathoperator{SL}(2,Z)$.

    INPUT
        k   -- weight
        m   -- index
        l   -- a Hecke operator index (rel. prime to the index $m$)
        n   -- an Atkin-Lehner involution index (i.e. exact divisor of $m$)
        kind -- the string 'holomorphic' or 'skew-holomorphic'

    NOTE

    EXAMPLES
    """
    k, m, eps = __check_jac_params( k, m, kind)
    k, m, l, n = __check( k, m, l, n)
    if 1 == k == eps:
        raise NotImplementedError, \
              "dimensions for skew-holomorphic cusp forms wt 1: not implemented"
    if k < 1:
        return 0
    return (_sz_s(k,m,l,n) - eps * (-1)**(k%2) * _sz_s(k,m,l,m//n))/2


 
def trace_jacobi_eienstein_series( k, m, l, n, kind = 'holomorphic'):
    r"""
    OUTPUT
        Trace of  $T(l) \circ W_n$  on the space
        of holomorphic or skew-holomorphic Jacobi Eisenstein series
        of weight $k$, index $m$ on $\mathoperator{SL}(2,Z)$.

    INPUT
        k   -- weight
        m   -- index
        l   -- a Hecke operator index (rel. prime to the index $m$)
        n   -- an Atkin-Lehner involution index (i.e. exact divisor of $m$)
        kind -- the string 'holomorphic' or 'skew-holomorphic'

    NOTE

    EXAMPLES
    """
    k, m, eps = __check_jac_params( k, m, kind)
    k, m, l, n = __check( k, m, l, n)
    raise NotImplementedError()


def trace_jacobi_forms( k, m, l, n, kind = 'holomorphic'):
    r"""
    OUTPUT
        Trace of  $T(l) \circ W_n$  on the space
        of holomorphic or skew-holomorphic Jacobi cusp forms
        of weight $k$, index $m$ on $\mathoperator{SL}(2,Z)$.

    INPUT
        k   -- weight
        m   -- index
        l   -- a Hecke operator index (rel. prime to the index $m$)
        n   -- an Atkin-Lehner involution index (i.e. exact divisor of $m$)
        kind -- the string 'holomorphic' or 'skew-holomorphic'

    NOTE

    EXAMPLES
    """
    return trace_jacobi_eienstein_series( k, m, l, n, kind = 'holomorphic') \
        + trace_jacobi_cusp_forms( k, m, l, n, kind = 'holomorphic')




# def dimension_jac_eis_series( k, m, kind = 'holomorphic'):
#     r"""
#     OUTPUT
#         Dimension of the space of Jacobi Eisenstein series
#         of weight $k$, index $m$ on $\mathsymbol{SL}(2,Z)$.

#     INPUT
#         k   -- weight
#         m   -- index
#         kind -- the string 'holomorphic' or 'skew-holomorphic'

#     NOTE

#     EXAMPLES
#     """
#     k, m, eps = __check_jac_params( k, m, kind)
#     if -1 != eps:
#         raise NotImplementedError, \
#               "dimensions for skew-holomorphic Eis. series: not implemented"
#     if k <= 1:
#         return 0
#     a = _ab( m)[0];
#     if 2 == k:
#         return a//2 + 1 - sigma(a,0)
#     if 0 == k%2:
#         return a//2 + 1
#     if 1 == k%2:
#         return (a-1)//2


# def dimension_jac_cusp_forms( k, m, eps = -1):
#     r"""
#     OUTPUT
#         Dimension of the space of Jacobi cusp form
#         of weight $k$, index $m$ on $\mathoperator{SL}(2,Z)$.

#     INPUT
#         k   -- weight
#         m   -- index
#         eps -- $-1$ for holomorphic and $+1$ for skew-holomorphic Jacobi forms

#     EXAMPLES
#     """
#     k, m, eps = __check_jac_params( k, m, eps)
#     return trace_jacobi_cusp_forms( 1, 1, k, m, eps)


# def dimension_jac_forms( k, m, eps = -1):
#     r"""
#     OUTPUT
#         Dimension of the space of all Jacobi forms
#         of weight $k$, index $m$ on $\mathoperator{SL}(2,Z)$.

#     INPUT
#         k   -- weight
#         m   -- index
#         eps -- $-1$ for holomorphic and $+1$ for skew-holomorphic Jacobi forms

#     EXAMPLES
#     """
#     k, m, eps = __check_jac_params( k, m, eps)
#     return dimension_jac_eis_series( k, m, eps) \
#            + dimension_jac_cusp_forms( k, m, eps );


# def poincare_poly_jac_forms( m, var):
#     r"""
#     OUTPUT
#         The polynomial $p(x)$ such that the Poincare series of
#         $\bigosum_{k \in Z} J_{k,m}^\eps$ equals $p(x)/(1-x^4)(1-x^6)$.

#     INPUT
#         m   -- index

#     NOTE
#        Multiplication with elliptic modular forms if $eps=-1$ or
#        turns $J := \bigosum_{k \in Z} J_{k,m}^\eps$ into a
#        graded free module over
#        the ring $M$ of elliptic modular forms on the full modular group.
#        The (Hilbert-)Poincare series of $M$ is accordingly of the form
#        $p(x)/(1-x^4)(1-x^6)$ with a polynomial $p(x)$. The degree of $p(x)$
#        is always $\le 12$ since [TODO]
#        The $k$-th coefficient of $p(x)$ equals the number of generators
#        of weight $k$ in a minimal system of homogeneous generators
#        of $J$ over $M$.

#     EXAMPLES
#     """
#     k, m, eps = __check_jac_params( 0, m, -1)
#     x = ZZ[var].gen();
#     a = dict([ (k,dimension_jac_forms(k,m)) for k in range(-10,13)])
#     return sum( (a[k]-a[k-4]-a[k-6]+a[k-10])*x**k for k in range(13))


# def poincare_poly_jac_cusp_forms( m, var):
#     r"""
#     OUTPUT
#         The polynomial $p(x)$ such that the Poincare series of
#         $\bigosum_{k \in Z} J_{k,m}^{\text{cusp},-1}$ equals
#         $p(x)/(1-x^4)(1-x^6)$.

#     INPUT
#         m   -- index

#     NOTE
#        Multiplication with elliptic modular forms turns
#        $S := \bigosum_{k \in Z} J_{k,m}^{\text{cusp},-1}$ into a
#        graded free module of rank $2m$ over the ring $M$
#        of elliptic modular forms
#        on the full modular group.  The (Hilbert-)Poincare series of
#        $M$ is accordingly of the form $p(x)/(1-x^4)(1-x^6)$ with a
#        polynomial $p(x)$. The degree of $p(x)$ is always $\le 12$
#        since [TODO] The $k$-th coefficient of $p(x)$ equals the number
#        of generators of weight $k$ in a minimal system of homogeneous
#        generators of $S$ over $M$.

#     EXAMPLES
#     """
#     k, m, eps = __check_jac_params( 0, m, -1)
#     x = ZZ[var].gen();
#     a = dict([ (k,dimension_jac_cusp_forms(k,m)) for k in range(-10,13)])
#     return sum( (a[k]-a[k-4]-a[k-6]+a[k-10])*x**k for k in range(13))
