#*****************************************************************************
#       Copyright (C) 2013 Nils-Peter Skoruppa <nils.skoruppa@gmail.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
CLASSES OF JACOBI FORMS OF LATTICE INDEX OVER Q
"""
_VERSION = '$Id$'


from sage.rings.integer import Integer
from sage.matrix.constructor import vector
from sage.misc.functional import is_odd



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
        return 'A Jacobi form in `J_{%d,L}( \epsilon^%d)`, where `L` has Gram matrix\n%s'\
            %(self.weight(), self.character(), latex(self.index().gram_matrix())) 

    def weight( self):
        return self.__weight

    def index( self):
        return self.__index

    def character( self):
        return self.__character

    def is_odd( self):
        return is_odd( Integer(self.weight() - self.character()/2))

    def _is_valid_C_index( self, D, r):
        """
        Return True if
        `D \ge 0`,
        `D + \beta(r) - h/24 is integral`,
        and otherwise False.
        """
        if D < 0 or False == self.index().is_in_shadow( r):
            return False
        try:
            a = Integer( D + self.index().beta(r) - self.character()/24)
            return True
        except:
            return False

    def _coefficient( self, D, r):
        raise NotImplementedError( 'not implemented')

    def coefficient( self, D, r):
        if self._is_valid_C_index( D, r):
            return self._coefficient( D, r)
        raise ValueError( '(%d,%s): not a valid index pair' %( D,r))

    def derivative( self):
        """
        Return the derivative `q\frac {d}{dq} \phi - \frac {k}{12} E_2 \phi`
        of this Jacobi form `\phi` (=\code{self}).
        """
        pass



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
        """
        Create the singular Jacobi form on the lattice `L`, whose `0`th coefficient
        equals the fiunction \code{inv}.
        """
        super( JacobiForm_singular_class, self).__init__( L.rank()/2, L, h) 
        self.__lambda = inv


    def _coefficient( self, D, r):
        return self.__lambda(r) if 0 == D else 0



class JacobiForm_pullback_class (JacobiForm_class):
    """
    A class for Jacobi forms initialized by pullback of a given Jacobi form.
    """
    def __init__( self, alpha, phi):
        self.__alpha = alpha
        sel.__phi = phi
        L = alpha.l_lattice()
        super( JacobiForm_pullback_class, self).__init__( phi.weight(), L, phi.character())


    def _coefficient( self, D, r):
        """

        NOTE
            The coefficient is the sum of the `C_{\phi}(D + L.\beta(r) - M.\beta(s),s)`,
            where  `s`  is in  `M^bullet`  such that
                o `M.\beta(s) <= D + L.beta(r)`,
                o `M.\beta(s-\alpha(r),\alpha(L)) = 0.
            Here `L` is self and `M` the codomain of `\alpha`.
        """
        alpha = self.__alpha
        L = alpha.l_module
        M = alpha.r_module
        D_r = D + L.beta(r)
        sv = M.shadow_vectors_in_shell( D_r)
        ar = alpha(r)
        # pick
        def is_ortogonal( x):
            for y in L.basis():
                if M.beta(x,y) != 0:
                    return False
            return True
        x_sv = filter( lambda s: is_ortogonal( s-ar), sv)
        return sum( self.__phi.coefficient( D_r - M.beta(s), s) for s in x_sv)



def JoliForm( **kwargs):
    """
    Return an instance of a class of Jacobi forms
    according to the initializing sequence of keyword arguments.

    For the moment, `kwargs` can be:
    o invariant = function --- an invariant of a Weil/shadow representations `W(L)`
    o embedding = an_embedding, form = phi --- an embedding alpha of a a lattice L and a Jacobi form on L
    """
    try:
        return JacobiForm_singular_class( kwargs['invariant'])
    except:
        pass

    try:
        return JacobiForm_pullback_class( kwargs['embedding'], kwargs['form']) 
    except:
        pass

    raise NotImplementedError()
