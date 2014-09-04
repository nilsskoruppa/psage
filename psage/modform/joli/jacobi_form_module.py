#*****************************************************************************
#       Copyright (C) 2013 Nils-Peter Skoruppa <nils.skoruppa@gmail.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
CLASS OF JACOBI FORMS OF LATTICE INDEX OVER Q AS MODULE OVER THE RING OF MODULAR FORMS
"""
_VERSION = '$Id$'

import dimensions
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sage_object import SageObject
from sage.rings.big_oh import O
from sage.misc.cachefunc import cached_method
from sage.rings.rational import Rational

class JoliModule_class (SageObject):
    """
    Model the module
    $J_{parity,L}(\varepsilon^h) = \bigoplus_{k-h/2 parity} J_{k,L}(\varepsilon^h)$
    over the ring of modular forms on $\SL(2\ZZ)$,
    """

    def __init__( self, L, h, parity = 'odd', uterm = None):
        self.__index = L
        self.__character = h
        self.__parity = parity
        self.__par = 1 if 'odd' == parity else 0
        ## The unknown dimensions occur in weight n/2, n/2 + 1/2, n/2 + 1, n/2 +3/2
        ## where n = rank L.
        ## Here we have the unknown dimension for
        ## k = n/2 + t, k - h/2 parity, i.e. for 2t = h-n + parity mod 4
        self.__uterm = uterm
        n = L.rank()
        self.__special_weight = n/2 + ((h - n + 2*self.__par)%4)/2

    
    def _latex_( self):
        return 'The module $J_{%s,L}(\\varepsilon^{%d})$,\nwhere L is the lattice with Gram matrix\n%s'\
            %(self.parity(), self.character(), latex( self.index().gram_matrix()))


    def _repr_( self):
        return 'The module J_{%s,L}(varepsilon^{%d}),\nwhere L is the lattice with Gram matrix\n%s'\
            %(self.parity(), self.character(), self.index().gram_matrix())


    def index( self):
        return self.__index


    def character( self):
        return self.__character


    def parity( self):
        return self.__parity


    def special_weight( self):
        return self.__special_weight


    def rank( self):
        L = self.index()
        o_inv = L.o_invariant()
        V2 = L.shadow_vectors_of_order_2()
        return Integer( (L.det() + len(V2)*(-1)**(self.__par + o_inv))/2)


    def set_uterm( self, d):
        self.__uterm = d


    @cached_method
    def dimension( self, k):
        try:
            t = Integer((k - self.special_weight())/2)
        except:
            raise ValueError( '%d: not in the range of grades'%k)
        if k < self.index().rank()/2:
            return Integer(0)
        elif self.special_weight() == k:
            if self.__uterm: return self.__uterm
            else: raise NotImplementedError( 'weight %s: not implemented'%k)
        else:
            return dimensions.dim_Joli( k, self.index(), self.character())


    def Poincare_pol( self, var = 'x', uterm = 0):
        """
        Return $s$ and the polynomial $p(x)$ such that
        \[
            x^s p(x)/(1-x^4)(1-x^6)
            =
            uterm
            +
            x^{s}
            \sum_{k=1}^\infty
            \dim J_{s + 2k,L}(\varepsilon^h) x^{2k}
        \]
        where $s$ is the (unique) half integer
        such that $s \equiv h/2 + parity \bmod 2\ZZ$
        and $n/2 \le s < n/2 +2$.

        OUTPUT
           We return the pair s, p(x)

        REMARK
           Note that we only have dimension formulas for weights $k \ge n/2 + 2$.
           Accordingly, to have a correct result, one should input uterm
           equal to the dimension of $J_{s,L}(\varepsilon^h)$.

        IMPLEMENTATION
           We use the Supplement to Theorem 4.1 in
           Boylan, H. and Skoruppa, N-P.: Jacobi forms of lattice index. I: Basic Theory.
           Namely, the Poincare polynomial has always degree $< n/2 + 12$, where
           $n$ is the rank of $self.index()$.
        """
        n  = self.index().rank()
        s = self.special_weight()
        v = PolynomialRing( IntegerRing(), var).gens()[0]
        top = Rational(n/2+12-s).ceil()
        a = dict([ (s+l,self.dimension(s+l)) for l in range(-10,top,2) if l != 0])
        a[s] = uterm
        Poincare_pol = s,\
            sum( (a[s+l]-a[s+l-4]-a[s+l-6]+a[s+l-10])* v**l for l in range(0,top,2))

        # L = self.index()
        # h = self.character()
        # s = self.special_weight()
        # N = 50 #TODO: does 13 suffice?
        # P = uterm + sum( self.dimension( s + 2*k) * x**(2*k) for k in range( 1, N)) + O(x**N)
        # P *= (1-x**4)*(1-x**6)
        # Poincare_pol = s, P.truncate()

        return Poincare_pol



def JoliModule( L, h, parity = 'odd', uterm = None):
    """
    Return an instance of class JoliModule_class.
    """
    return JoliModule_class( L, h, parity = parity, uterm = uterm)
