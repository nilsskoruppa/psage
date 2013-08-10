#*****************************************************************************
#       Copyright (C) 2013 Nils-Peter Skoruppa <nils.skoruppa@gmail.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
_"""
Provides classes and various functions to do computations
in the context of Jacobi forms of lattice index over $Q$.

Copyright (C) Nils-Peter Skoruppa 2013
"""
VERSION = '$Id$'


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



class JF_Module_class (SageObject):
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
        self.__Poincare_pol = None

    
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
        o_inv = L.o_invariant()
        V2 = L.shadow_vectors_of_order_2()
        return Integer( (L.det() + len(V2)*(-1)^(self.__par + o_inv))/2)


    def set_uterm( self, d):
        self.__uterm = d


    def dimension( self, k):
        try:
            t = Integer((k - self.special_weight())/2)
        except:
            raise ValueError( '%d: not in the range of grades'%k)
        if k < self.index().rank()/2:
            return Integer(0)
        elif self.special_weight() == k:
            if self.__uterm: return self.__uterm
            else: raise NotImplementedError( 'weight %d: not implemented'%k)
        else:
            return dim_Jac( k, self.index(), self.character())


    def Poincare_pol( self, x, uterm = 0):
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
        """
        if self.__Poincare_pol:
            return self.__Poincare_pol

        L = self.index()
        h = self.character()
        s = self.special_weight()
        N = 50
        P = uterm + sum( dim_Jac( s + 2*k, L, h)*x^(2*k) for k in range( 1, N)) + O(x^N)
        P *= (1-x^4)*(1-x^6)

        self.__Poincare_pol = s, P.truncate()
        return self.__Poincare_pol



class JacobiForm_class (SageObject):
    """
    A template class for a Jacobi form.
    """

    def __init__( self, k, L, h):
        self.__weight = k
        self.__index = L
        self.__character = h

    def _repr_(self):
        return 'A Jacobi form in J_{%d,L}( epsilon^%d), where L has Gram matrix\n%s'\
            %(self.weight(), self.character(), self.index().gram_matrix()) 

    def _latex_(self):
        return 'A Jacobi form in $J_{%d,L}( \epsilon^%d)$, where $L$ has Gram matrix\n%s'\
            %(self.weight(), self.character(), latex(self.index().gram_matrix())) 

    def weight( self):
        return self.__weight

    def index( self):
        return self.__index

    def character( self):
        return self.__character

    def _is_valid_C_index( self, D, r):
        """
        Return True if
        $D \ge 0$,
        $D + \beta(r) - h/24 is integral$,
        and otherwise False.
        """
        if D < 0 or False == self.index().is_in_shadow( r):
            return False
        try:
            a = Integer
            a = Integer( D + self.index().beta(r) - self.character()/24)
            return True
        except:
            return False

    def _coefficient( self, D, r):
        raise NotImplementedError( 'not implemented')

    def coefficient( self, D, r):
        if 0 != D:
            return 0
        if self._is_valid_C_index( D, r):
            return self._coefficient( D, r)
        raise ValueError( '(%d,%s): not a valid index pair' %( D,r))
    


class JacobiForm_singular_class (JacobiForm_class):
    """
    A class for Jacobi forms initialized by invariants.

    EXAMPLES
        sage: phi = JacobiForm_singular_class( L, 9, inv); phi
        A Jacobi form in J_{1,L}( epsilon^9), where L has Gram matrix
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: phi.coefficient( 0, L.a_shadow_vector())        
        1
        sage: L = Lattice_class( [1,0,0,1,0,1])
        sage: inv = lambda x: kronecker_symbol( -4, Integer(8*x[0]*x[1]*x[2]))
        sage: phi.coefficient( 0, L.a_shadow_vector()); phi
        1
        A Jacobi form in J_{1,L}( epsilon^9), where L has Gram matrix
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: phi.coefficient( 0, L.a_shadow_vector())
        1

    """
    def __init__( self, L, h, inv):
        super( JacobiForm_singular_class, self).__init__( L.rank()/2, L, h) 
        self.__lambda = inv

    def _coefficient( self, D, r):
        return self.__lambda(r)


class JacobiForm_pullback_class (JacobiForm_class):
    """
    An class for Jacobi forms initialized by pullback of a given Jacobi form.
    """
    def __init__( self, a, L, M, phi):
        super( JacobiForm_singular_class, self).__init__( phi.weight(), L, phi.character())

    def __coefficient( self, D, r):
        """

        NOTE
            The coefficient is the sum of the C_phi(D-M.beta(s-a(r)),s),
            where  s  is in  M^bullet  such that
                o D >= M.beta(s-a(r)),
                o (s-a(r))a = 0.
        """
        ar = a*r
        b = 30 #TODO determine reasonable b
        l = L.shadow_vectors_in_shell( b)
        l += [-s for s in l] 
        l.append( vector( [0]*L.rank())
        l1 = map( lambda x: x-ar, l)
        l1 = filter( lambda x: M.beta(x) <= D and 0 == x*a, l)
        return sum( phi.coefficient( D-M.beta(s-ar),s) for s in l)

