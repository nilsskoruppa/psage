#*****************************************************************************
#       Copyright (C) 2013 Nils-Peter Skoruppa <nils.skoruppa@gmail.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
LATTICE CLASS FOR JACOBI FORMS OF LATTICE INDEX

Provides classes and various functions to do computations
in the context of Jacobi forms of lattice index over $Q$.
"""
_VERSION = '$Id$'


from psage.modules.finite_quadratic_module import *

# load '/home/nils/Sandbox/Articles/Present/Joli-I/Computations/lattice_data.sage'

from sage.all import *

sage.misc.preparser.load(sage.misc.preparser.base64.b64decode("Jy9ob21lL25pbHMvU2FuZGJveC9BcnRpY2xlcy9QcmVzZW50L0pvbGktSS9Db21wdXRhdGlvbnMvbGF0dGljZV9kYXRhLnNhZ2Un"),globals(),False)


class Lattice_class (SageObject):
    """
    A class representing a lattice $L$.

    EXAMPLES
        sage: L = Lattice_class( [2,1,2]); L
        The ZZ-lattice (ZZ^2, x^tGy), where G = 
        [2 1]
        [1 2]
        sage: L.is_even()
        True
        sage: L.gram_matrix()
        [2 1]
        [1 2]
        sage: L.dual_vectors()
        {(2/3, -1/3), (0, 0), (4/3, -2/3)}
        sage: L.representatives(2)
        {}
        sage: L.representatives(1)
        {(2/3, -1/3), (4/3, -2/3)}
        sage: L.values()          
        {0: {(0, 0)}, 1: {(2/3, -1/3), (4/3, -2/3)}}
        sage: L.det()                 
        3

        sage: g = lambda n: [1] if 1 == n else [1] + (n-1)*[0] + g(n-1) 
        sage: Z8 = Lattice_class( g(8))
        sage: a,A8 = Z8.ev()
        sage: M8 = A8.fqm()
        sage: M8.jordan_decomposition().genus_symbol()
        '2^2'
    """
    
    def __init__( self, q):
        """
        We initialize by a list of integers $[a_1,...,a_N]$. The
        lattice self is then $L = (L,\beta) = (G^{-1}\ZZ^n/\ZZ^n, G[x])$,
        where $G$ is the symmetric matrix
        $G=[a_1,...,a_n;*,a_{n+1},....;..;* ... * a_N]$,
        i.e.~the $a_j$ denote the elements above the diagonal of $G$.
        """
        self.__form = q
        # We compute the rank
        N = len(q)
        assert is_square( 1+8*N)
        n = Integer( (-1+isqrt(1+8*N))/2)
        self.__rank = n
        # We set up the Gram matrix
        self.__G = matrix( ZZ, n, n)
        i = j = 0
        for a in q:
            self.__G[i,j] = self.__G[j,i] = a
            if j < n-1:
                j += 1
            else:
                i += 1
                j = i
        # We compute the level
        Gi = self.__G**-1
        self.__Gi = Gi
        a = lcm( map( lambda x: x.denominator(), Gi.list()))
        I = Gi.diagonal()
        b =  lcm( map( lambda x: (x/2).denominator(), I))
        self.__level = lcm( a, b)
        # We compute a shadow vector
        self.__shadow_vector = vector((a%2)/2 for a in self.__G.diagonal())*Gi
        # We define a basis
        M = Matrix( ZZ, n, n, 1)
        self.__basis = M.rows()
        # We prepare a cache
        self.__dual_vectors = None
        self.__values = None
        self.__chi = {}


    def _latex_( self):
        return 'The ZZ-lattice $(ZZ^{%d}, x\'Gy)$, where $G = %s$' % (self.__rank, latex(self.__G))

    
    def _repr_( self):
        return 'The ZZ-lattice (ZZ^%d, x^tGy), where G = \n%r' % (self.__rank, self.__G)


    def rank( self):
        return self.__rank


    def level( self):
        """
        Return the smallest integer $l$ such that $lG[x]/2 \in \Z$
        for all $x$ in $L^*$.
        """
        return self.__level


    def _basis( self):
        return self.__basis


    def is_positive( self):
        pass


    def is_even( self):
        I = self.gram_matrix().diagonal()
        for a in I:
            if is_odd(a):
                return False
        return True


    def is_maximal( self):
        pass


    def shadow_level( self):
        """
        Return the level of $L_ev$, where $L_ev$ is the kernel
        of the map $x\mapsto e(G[x]/2)$.

        REMARK
            Lemma: If $L$ is odd, if $s$ is a shadow vector, and if $N$ denotes the
            denominator of $G[s]/2$, then the level of $L_ev$
            equals the lcm of the order of $s$, $N$ and the level of $L$.

            (For the proof:  the requested level is the smallest integer $l$
            such that $lG[x]/2$ is integral for all $x$ in $L^*$ and $x$ in $s+L^*$
            since $L_ev^* = L^* \cup s+L^*$. This $l$ is the smallest integer
            such that $lG[s]/2$, $lG[x]/2$ and $ls^tGx$ are integral for all $x$ in $L^*$.)
        """
        if self.is_even():
            return self.level()
        s = self.a_shadow_vector()
        N = denominator((s*self.gram_matrix()*s)/2)
        h = lcm([denominator(x) for x in s])
        return lcm( [h, N, self.level()])


    def gram_matrix( self):
        return self.__G


    def beta(self, x, y = None):
        if None == y:
            return (x*self.gram_matrix()*x)/2
        else:
            return x*self.gram_matrix()*y


    def det( self):
        return self.gram_matrix().det()


    def dual_vectors( self):
        """
        Return a set of representatives
        for $L^#/L$.
        """
        D,U,V = self.gram_matrix().smith_form()
        # hence D = U * self.gram_matrix() * V
        if None == self.__dual_vectors:
            W = V*D**-1
            self.__dual_vectors = [ W*vector(v) for v in mrange( D.diagonal())]
        return self.__dual_vectors


    def is_in_dual( self, y):
        """
        Checks whether $y$ is in the dual of self.
        """
        return self._is_int_vec( self.gram_matrix()*y)


    def a_shadow_vector( self):
        """
        Return a vector in the shadow $L^\bullet$ of $L$.

        REMARK
            As shadow vector one can take the diagnal of the Gram
            matrix reduced modulo 2 and multiplied by $1/2G^{-1}$.

            Note that this vector $s$ is arbitrary (in particular,
            it does not satisfy in general $2s \in L$).
        """
        return self.__shadow_vector


    def shadow_vectors( self):
        """
        Return a list of representatives for the vectors in $L^\bullet/L$.
        """
        if self.is_even():
            return self.dual_vectors() 
        R = self.dual_vectors()
        s = self.a_shadow_vector()
        return [s+r for r in R]


    def shadow_vectors_of_order_2( self):
        """
        Return the list of representatives for those
        vectors $x$ in $L^\bullet/L$ such that
        $2x=0$.
        """
        R =  self.shadow_vectors()
        return [r for r in R if self._is_int_vec( 2*r)]
        

    def is_in_shadow( self, r):
        """
        Checks whether $r$ is a shadow vector.
        """
        c = self.a_shadow_vector()
        return self.is_in_dual( r - c)


    def o_invariant( self):
        """
        Return $0$ if $L$ is even, otherwise
        the parity of $G[2*s]$, where $s$ is a shadow vector such that
        $2s$ is in $L$ (so that $(-1)^parity = e(G[2*s]/2)$.
        
        REMARK
            The o_invariant equals the parity of the Gram matrix. 
        """
        V = self.shadow_vectors_of_order_2()
        s = V[0]
        return Integer(4*(s*self.gram_matrix()*s))%2
 

    def values( self):
        """
        Return a dictionary
        N \mapsto Set of representatives r in $L^\bullet/L$ which satisfy
        $\beta(r) \equiv n/l \bmod \ZZ$, where $l$ is the shadow_level
        (i.e. the level of $L_ev$).

        REMARK
            $\beta(r)+ \Z$ does only depend on $r+L$ if $r$ is a shadow vector,
            otherwise this is not necessarily true. Hence we return only
            shadow vectors here.
        """
        if None == self.__values:
            self.__values = {}
            G = self.gram_matrix()
            l = self.shadow_level()
            R =  self.dual_vectors()
            if not self.is_even():
                s = self.a_shadow_vector()
                R = [s+r for r in R]
            for r in R:
                N = Integer((r*G*r)*l/2)
                N = N%l
                if not self.__values.has_key(N):
                    self.__values[N] = [r]
                else:
                    self.__values[N] += [r]
        return self.__values


    def representatives( self, N):
        """
        Return the subset of representatives $r$
        in $L^\bullet/L$ which satisfy $\beta(r) \equiv N/level \bmod \ZZ$.
        """
        v = self.values()
        return v.get(N%self.shadow_level(), [])
        

    def chi( self, t):
        """
        Return the value of the Gauss sum
        $\sum_{x\in L^\bullet/L} e(tG[x]/2)/\sqrt D$,
        where $D = |L^\bullet/L|$.
        """
        t = Integer(t)%self.shadow_level()
        if self.__chi.has_key(t):
            return self.__chi[t]
        l = self.shadow_level()
        z = QQbar.zeta(l)
        w = QQbar(sqrt(self.det()))
        V = self.values()
        self.__chi[t] = sum(len(V[a])*z**(t*a) for a in V)/w
        return self.__chi[t]


    def ev( self):
        """
        Return the pair $a, L_{ev}$, where $L_{ev}$
        is the kernel of $L\rightarrow \{\pm 1\}$, 
        $x\mapsto e(\beta(x))$,
        and where $a$ are the vectors in $L$ which
        correspond to the canonical basis for $L_{ev}$.

        REMARK
            We have to find the kernel of the map
            $x\mapsto G[x] \equiv \sum_{j\in S}x_2 \bmod 2$.
            Here $S$ is the set of indices $i$ such that
            the $i$-diagonal element of $G$ is odd.
            A basis is given by
            $e_i$ ($i not\in S$) and $e_i + e_j$ ($i \in S$, $i\not=j$)
            and $2e_j$, where $j$ is any fixed index in $S$.
        """
        if self.is_even():
            return matrix( ZZ, self.rank(), self.rank(), 1), self
        D = self.gram_matrix().diagonal()
        S = [ i for i in range( len(D)) if is_odd(D[i])]
        j = min(S)
        e = self._basis()
        n = len(e)
        form = lambda i: e[i] if i not in S else 2*e[j] if i == j else e[i]+e[j] 
        a = [ form(i) for i in range( n)]
        A = matrix( a).transpose()
        return A, Lattice_class( [ self.beta( a[i],a[j]) for i in range(n) for j in range(i,n)])
 

    def fqm( self):
        """
        Return a pair $(f,M)$, where $M$ is the discriminant module
        of $L_ev$ up to isomorphism and $f:L\rightarrow M$ the
        canonical map.

        TODO
            Currently returns only M.
        """
        if self.is_even():
            return FiniteQuadraticModule( self.gram_matrix())
        else:
            a,Lev = self.ev()
            return FiniteQuadraticModule( Lev.gram_matrix())


    def ZN_embeddings( self):
        """
        Return a list of all embeddings of $L$ into $\ZZ^N$ modulo
        the action of $O(\ZZ^N)$.
        """
        def cs_range( f, subset = None):
            """
            For a symmetric semi-positive integral matrix $F$,
            return a list of all integral $n$-vectors $v$ such that
            $x^tfx - (v*x)^2 >= 0$ for all $x$.
            """
            n = f.dimensions()[0]
            b = vector( isqrt(s) for s in f.diagonal())
            zv = vector([0]*n)
            box = [b - vector(t) for t in mrange( (2*b + vector([1]*n)).list()) if b - vector(t) > zv]
            if subset:
                box = [v for v in box if v in subset]
            rge = [v for v in box if min( (f - matrix(v).transpose()*matrix(v)).eigenvalues()) >=0]
            return rge

        def embs( f, subset = None):
            """
            For a semipositive matrix $f$ return a list of all integral
            matrices $M$ such that $f=M^t * M$ modulo the left action
            of $\lim O(\ZZ^N)$ on the set of these matrices.
            """
            res = []
            csr = cs_range( f, subset)
            for v in csr:
                fv = f - matrix(v).transpose()*matrix(v)
                if 0 == fv:
                    res.append( matrix(v))
                else:
                    tmp = embs( fv, csr[ csr.index( v):])
                    for e in tmp:
                        res.append( e.insert_row(0, v))
            return res 

        return embs( self.gram_matrix())


    def vectors_in_shell ( self, bound, max = 10000):
        """
        Return the list of all $x\not=0$ in self
        modulo multiplication by $-1$ such that
        $\beta(x) <= bound$.
        """
        return self._vectors_in_shell ( self.gram_matrix(), 2*bound, max)


    def dual_vectors_in_shell ( self, bound, max = 10000):
        """
        Return the list of all $x\not=0$ in the dual of self
        modulo multiplication by $-1$ such that
        $\beta(x) <= bound$.
        """
        l = self.level()
        F = l*self.__Gi
        return map( lambda x: self.__Gi*x, self._vectors_in_shell ( F, 2*l*bound, max))


    def shadow_vectors_in_shell (self, bound, max = 10000):
        """
        Return the list of all $x\not=0$ in the bullet of self
        modulo multiplication by $-1$ such that
        $\beta(x) <= bound$.
        """
        if self.is_even():
            return self.dual_vectors_in_shell ( bound, max)
        a, Lev = L.ev()
        l = Lev.dual_vectors_in_shell ( bound, max)
        return filter( lambda x: self.is_in_shadow(x) ,[a*y for y in l])


    @staticmethod
    def _vectors_in_shell( F, bound, max = 1000):
        """
        For an (integral) positive symmetric matrix $F$
        return a list of vectors such that $x^tFx <= bound$.
        The method returns only nonzero vestors and for each pair $v$ and $-v$
        only one.

        NOTE
            A wrapper for the pari method qfminim().
        """
        # assert bound >=1, 'Pari bug: bound should be >= 2'
        if bound < 1: return []
        p = F._pari_()
        po = p.qfminim( bound, max, flag = 0).sage()
        if max == po[0]:
            raise ArithmeticError( 'Increase max')
        return po[2].columns()

    
    @staticmethod
    def _is_int_vec( v):
        for c in v:
            if c.denominator() > 1:
                return False
        return True
