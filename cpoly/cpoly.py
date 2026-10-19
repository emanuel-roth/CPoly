from sage.all import RootSystem as BaseRootSystem
from sage.all import CartanType
from sage.all import CartanMatrix
from sage.all import matrix, vector, diagonal_matrix, QQ
from sage.combinat.subset import Subsets
from itertools import combinations

def greater_than(A, B):
    r"""Compares two lists, and distinguishes when elements of :math:`B` are in :math:`A`

    INPUT: Tuples :math:`A` and :math:`B`

    OUTPUT: ``True`` if elements of :math:`B` are in :math:`A`, and ``False`` otherwise

    EXAMPLE::

        sage: from cpoly import RootSystem, greater_than

        sage: A2 = RootSystem(['A', 2])
        sage: A = A2.facets()[0]
        sage: B = A2.facets()[1]
        sage: greater_than(A,B)
        True

    """
    return all(x in B for x in A)

class RootSystem(BaseRootSystem):

    def weyl_chambers(self):
        r"""Returns the Weyl chambers of the root system. The Weyl chambers are represented by the images of the tuple of fundamental weights under the action of the Weyl group.  More precisely, if :math:`\Lambda_1, \ldots, \Lambda_n` are the fundamental weights and ``W`` is the Weyl group, this method returns :math:`\{(w(\Lambda_1), \ldots, w(\Lambda_n)) : w \in W\}`. Each chamber is represented as a tuple whose entries are the fundamental weights that span it as a cone.

        OUTPUT: :math:`\{(w(\Lambda_1), \ldots, w(\Lambda_n)) : w \in W\}`

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: A2.weyl_chambers()
            [(Lambda[1], Lambda[2]),
            (-Lambda[1] + Lambda[2], -Lambda[1]),
            (-Lambda[2], Lambda[1] - Lambda[2]),
            (-Lambda[1] + Lambda[2], Lambda[2]),
            (Lambda[1], Lambda[1] - Lambda[2]),
            (-Lambda[2], -Lambda[1])]

        """
        WL = self.weight_lattice()
        W = WL.weyl_group()
        Lambda = WL.fundamental_weights()
        return [tuple(w.action(Lambda[i]) for i in Lambda.keys()) for w in W]

    def facets(self):
        r"""Returns the facets of the root system. The facets are the faces of each Weyl chamber, including the trivial facet :math:`\{0\}`. Each facet is represented as a tuple of its verticies, i.e., the fundamental weights that span it as a cone.

        OUTPUT: Facets as a list of tuples of their vertices

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: A2.facets()
            [(),
            (Lambda[1],),
            (-Lambda[1] + Lambda[2],),
            (-Lambda[2],),
            (Lambda[2],),
            (-Lambda[1],),
            (Lambda[1] - Lambda[2],),
            (Lambda[1], Lambda[2]),
            (-Lambda[1] + Lambda[2], -Lambda[1]),
            (-Lambda[2], Lambda[1] - Lambda[2]),
            (-Lambda[1] + Lambda[2], Lambda[2]),
            (Lambda[1], Lambda[1] - Lambda[2]),
            (-Lambda[2], -Lambda[1])]

        """
        WL = self.weight_lattice()
        W = WL.weyl_group()
        Lambda = list(WL.fundamental_weights().values())
        fundamental_facets = [tuple(S) for r in range(1, len(Lambda) + 1) for S in Subsets(Lambda, r)]
        facets = [()]
        seen = set()
        for facet in fundamental_facets:
            for w in W:
                image = tuple(w.action(weight) for weight in facet)
                key = tuple(sorted(image, key=str))
                if key not in seen:
                    seen.add(key)
                    facets.append(image)
        return facets

    def facets_one_dimensional(self):
        r"""Returns the one-dimensional facets of the root system.
                
        The facets are the faces of each Weyl chamber, of which one-dimensional facets are spanned by one vertex.
        
        OUTPUT: One-dimensional facets as a list of tuples of their vertices

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: A2.facets_one_dimensional()
            [(Lambda[1],),
             (-Lambda[1] + Lambda[2],),
             (-Lambda[2],),
             (Lambda[2],),
             (-Lambda[1],),
             (Lambda[1] - Lambda[2],)]
        
        """
        WL = self.weight_lattice()
        W = WL.weyl_group()
        Lambda = list(WL.fundamental_weights().values())

        facets = []
        seen = set()

        for weight in Lambda:
            for w in W:
                image = (w.action(weight),)
                key = tuple(sorted(image, key=str))

                if key not in seen:
                    seen.add(key)
                    facets.append(image)

        return facets
    
    def weyl_chambers_smaller_than_facet(self, facet):
        r"""Returns the Weyl chambers smaller than given facet
                                        
        OUTPUT: The Weyl chambers smaller than given facet
        
        EXAMPLE::
        
            sage: from cpoly import RootSystem, greater_than
        
            sage: A2 = RootSystem(['A', 2])
            sage: Facet = A2.facets()[1]
            sage: A2.weyl_chambers_smaller_than_facet(Facet)
            [(Lambda[1], Lambda[2]), (Lambda[1], Lambda[1] - Lambda[2])]
        
        """
        WL = self.weight_lattice()
        W = WL.weyl_group()
        Lambda = list(WL.fundamental_weights().values())
        return [tuple(w.action(lam) for lam in Lambda) for w in W if greater_than(facet, tuple(w.action(lam) for lam in Lambda))]

    def euclidean_product(self, v, w):
        r"""Compute the Euclidean inner product of two elements of the root or weight lattice. The bilinear form is obtained from the Cartan matrix ``A`` and its symmetrizer ``D``.

        INPUT:
            - ``v`` -- an element of the root or weight lattice
            - ``w`` -- an element of the root or weight lattice
            
        OUTPUT: An element of `\QQ` representing the Euclidean inner product of ``v`` and ``w``

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: V = A2.root_space()
            sage: alpha = V.simple_roots()
            sage: A2.euclidean_product(alpha[1], alpha[2])
            -1
        
        """
        R = self.root_space()
        v = R(v)
        w = R(w)

        A = matrix(QQ, self.cartan_matrix())
        D = diagonal_matrix(QQ, self.cartan_type().symmetrizer())
        G = D * A
        return v.to_vector() * G * w.to_vector()

    def weight_to_root_space(self, weight):
        r"""Express a weight as an element of the root space. Given a weight ``\lambda``, this method computes its coordinates with respect to the basis of simple roots and returns the corresponding element of the root space.

        INPUT: ``weight`` -- an element of the weight lattice
        
        OUTPUT: An element of the root space representing ``weight`` in the basis of simple roots

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: Lambda = W.fundamental_weights()
            sage: A2.weight_to_root_space(Lambda[1])
            2/3*alpha[1] + 1/3*alpha[2]
        
        """
        R = self.root_space()
        A = matrix(QQ, self.cartan_matrix())
        c = weight.to_vector()
        r = A.transpose().inverse() * c
        simple_roots = list(R.simple_roots().values())
        return sum(
            r[i] * simple_roots[i]
            for i in range(len(simple_roots))
        )

    def is_domain_range_comp_poly(self, comp_poly):
        r"""Tests whether a function is evaluated on the Weyl chambers, into the root space, similar to a complementary polyhedron. A valid ``comp_poly`` is a dictionary whose keys are Weyl chambers
        of ``self`` and whose values are elements of the root space of ``self``. We identify the root space with its dual through the canonical isomorphism from the inner product, so this is isomorphic to Behrend's definition.

        INPUT: ``comp_poly`` -- a dictionary mapping Weyl chambers to elements of the root space

        OUTPUT: ``True`` if ``comp_poly`` has the required form, and ``False`` otherwise
        
        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: A2.is_domain_range_comp_poly(comp_poly)
            True
        
        """
        if not isinstance(comp_poly, dict):
            return False

        chambers = self.weyl_chambers()
        R = self.root_space()

        return (set(comp_poly.keys()) == set(chambers) and all(R(v) in R for v in comp_poly.values()))

    def has_A1(self, comp_poly):
        r"""Tests whether a complementary polyhedron ``comp_poly`` has property ``A1``, as laid out in Behrend's paper.

        INPUT: ``comp_poly`` -- a dictionary mapping Weyl chambers to elements of the root space
        
        OUTPUT: ``True`` if ``comp_poly`` fulfills ``A1``, and ``False`` otherwise

        EXAMPLE::
                        
            sage: from cpoly import RootSystem, greater_than
                        
            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: A2.has_A1(comp_poly)
            True

        """
        for facet in self.facets_one_dimensional():
            facet_vector = self.weight_to_root_space(facet[0])
            chambers = self.weyl_chambers_smaller_than_facet(facet)
            for chamber1, chamber2 in combinations(chambers, 2):
                value1 = self.euclidean_product(comp_poly[chamber1], facet_vector)
                value2 = self.euclidean_product(comp_poly[chamber2], facet_vector)
                if value1 != value2:
                    return False
        return True

    def has_A2(self, comp_poly):
        r"""Tests whether a complementary polyhedron ``comp_poly`` has property ``A2``, as laid out in Behrend's paper.
        
        INPUT: ``comp_poly`` -- a dictionary mapping Weyl chambers to elements of the root space
                
        OUTPUT: ``True`` if ``comp_poly`` fulfills ``A2``, and ``False`` otherwise

        EXAMPLE::
                                
            sage: from cpoly import RootSystem, greater_than
                                
            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: A2.has_A2(comp_poly)
            True
        
        """
        roots = self.root_space().roots()
        for facet in self.facets_one_dimensional():
            chambers = self.weyl_chambers_smaller_than_facet(facet)
            for chamber1, chamber2 in combinations(chambers, 2):
                for root in roots:
                    positive = all(self.euclidean_product(root, self.weight_to_root_space(weight)) > 0 for weight in chamber1)
                    negative = all(self.euclidean_product(root, self.weight_to_root_space(weight)) < 0 for weight in chamber2)
                    if positive and negative:
                        if self.euclidean_product(comp_poly[chamber1], root) > self.euclidean_product(comp_poly[chamber2], root):
                            return False
        return True

    def is_comp_poly(self, comp_poly):
        r"""Tests whether ``comp_poly`` satisfies the conditions for being a complementary polyhedron. The following three conditions are checked:

            1. ``comp_poly`` is a dictionary whose keys are Weyl chambers and whose values are elements of the root space.
            2. ``comp_poly`` satisfies property :math:`A1`.
            3. ``comp_poly`` satisfies property :math:`A2`.

        INPUT: ``comp_poly`` -- a dictionary mapping Weyl chambers to elements of the root space

        OUTPUT: ``True`` if ``comp_poly`` satisfies all three conditions, and ``False`` otherwise
        
        EXAMPLE::
                                        
            sage: from cpoly import RootSystem, greater_than
                                        
            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: A2.is_comp_poly(comp_poly)
            True
                
        """
        return (self.is_domain_range_comp_poly(comp_poly) and self.has_A1(comp_poly) and self.has_A2(comp_poly))

    def R_facet(self, facet):
        r"""Returns the roots that are nonnegative with respect to every vertex of ``facet``.
    
        INPUT: ``facet`` -- a facet represented as a tuple of vertices in the weight lattice
    
        OUTPUT: A list of roots ``alpha`` that are nonnegative with respect to ``facet``
    
        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: Facet = A2.facets()[1]
            sage: A2.R_facet(Facet)
            [alpha[1], alpha[2], alpha[1] + alpha[2], -alpha[2]]
    
        """
        roots = self.root_space().roots()
        return [root for root in roots if all(self.euclidean_product(root, self.weight_to_root_space(vertex)) >= 0 for vertex in facet)]

    def U_facet(self, facet):
        r"""Returns the roots that are positive with respect to every vertex of ``facet``.

        INPUT: ``facet`` -- a facet represented as a tuple of vertices in the weight lattice

        OUTPUT: A list of roots ``alpha`` that are positive with respect to ``facet``

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: Facet = A2.facets()[1]
            sage: A2.U_facet(Facet)
            [alpha[1], alpha[1] + alpha[2]]

        """
        roots = self.root_space().roots()
        return [root for root in roots if all(self.euclidean_product(root, self.weight_to_root_space(vertex)) > 0 for vertex in facet)]

    def degree(self, comp_poly, facet):
        r"""Returns the degree associated to ``facet``. The degree is computed by choosing a Weyl chamber smaller than ``facet`` and summing the Euclidean products of its image under ``comp_poly`` with all roots in ``U_facet(facet)``.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron
            - ``facet`` -- a facet

        OUTPUT: A rational number

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[3]
            sage: A2.degree(comp_poly, Facet)
            -3
        
        """
        chambers = self.weyl_chambers_smaller_than_facet(facet)
        chamber = chambers[0]
        roots = self.U_facet(facet)
        return sum(self.euclidean_product(comp_poly[chamber], alpha) for alpha in roots)

    def is_semistable(self, comp_poly, check=True):
        r"""Tests whether ``comp_poly`` is semistable. A complementary polyhedron is semistable if the degree of every one-dimensional facet is less than or equal to zero.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron
            - ``check`` -- checks if ``comp_poly`` is a complementary polyhedron or not


        OUTPUT: ``True`` if ``comp_poly`` is a semistable complementary polyhedron, and ``False`` otherwise

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: A2.is_semistable(comp_poly)
            False
        
        """
        if check and not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")
        return all(self.degree(comp_poly, facet) <= 0 for facet in self.facets_one_dimensional())

    def is_stable(self, comp_poly, check=True):
        r"""Tests whether ``comp_poly`` is stable. A complementary polyhedron is stable if the degree of every one-dimensional facet is less than zero.
    
        INPUT:
            -``comp_poly`` -- a complementary polyhedron
            - ``check`` -- checks if ``comp_poly`` is a complementary polyhedron or not
    
        OUTPUT: ``True`` if ``comp_poly`` is a semistable complementary polyhedron, and ``False`` otherwise

        EXAMPLE::
        
            sage: from cpoly import RootSystem, greater_than
        
            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: A2.is_stable(comp_poly)
            False
        
        """
        if check and not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")
        return all(self.degree(comp_poly, facet) < 0 for facet in self.facets_one_dimensional())

    def psi(self, facet, vertex):
        r"""Returns the roots ``alpha`` such that :math:`(\alpha, \text{vertex}) = 1` and :math:`(\alpha, \mu) = 0` for every other vertex ``mu`` in ``facet``.

        INPUT:
            - ``facet`` -- a facet represented as a tuple of vertices
            - ``vertex`` -- a vertex of ``facet``

        OUTPUT: A list of roots satisfying the above conditions.

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: Lambda = A2.weight_lattice().fundamental_weights()
            sage: Facet = A2.facets()[3]
            sage: A2.psi(Facet, -Lambda[2])
            [-alpha[2], -alpha[1] - alpha[2]]
        
        """
        if vertex not in facet:
            raise ValueError("vertex must be an element of facet")
        roots = self.root_space().roots()
        return [alpha for alpha in roots if self.euclidean_product(alpha, self.weight_to_root_space(vertex)) == 1 and all(self.euclidean_product(alpha, self.weight_to_root_space(mu)) == 0 for mu in facet if mu != vertex)]

    def numerical_invariant(self, comp_poly, facet, vertex):
        r"""Returns the degree associated to ``facet``. The degree is computed by choosing a Weyl chamber smaller than ``facet`` and summing the Euclidean products of its image under ``comp_poly`` with all roots in ``U_facet(facet)``.
        
        INPUT:
            - ``comp_poly`` -- a complementary polyhedron
            - ``facet`` -- a facet
            - ``vertex`` -- a vertex of ``facet``
        
        OUTPUT: A rational number

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[3]
            sage: A2.numerical_invariant(comp_poly, Facet, -Lambda[2])
            -3
        
        """
        if vertex not in facet:
            raise ValueError("vertex must be an element of facet")
        chambers = self.weyl_chambers_smaller_than_facet(facet)
        chamber = chambers[0]
        roots = self.psi(facet, vertex)
        return sum(self.euclidean_product(comp_poly[chamber], alpha) for alpha in roots)

    def F_facet(self, comp_poly, facet):
        r"""Returns the restriction of ``comp_poly`` to the Weyl chambers smaller than ``facet``.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron
            - ``facet`` -- a facet

        OUTPUT: The restriction of ``comp_poly`` to the Weyl chambers smaller than ``facet``

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[3]
            sage: A2.F_facet(comp_poly, Facet)
            {(-Lambda[2], -Lambda[1]): alpha[1] + alpha[2],
             (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2]}
        
        """
        if not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")
        chambers = self.weyl_chambers_smaller_than_facet(facet)
        return {chamber: comp_poly[chamber] for chamber in chambers}

    def y_facet(self, comp_poly, facet):
        r"""Computes the y value of a facet with respect to a complementary polyhedron.
        
        INPUT:
            - ``comp_poly`` -- a complementary polyhedron
            - ``facet`` -- a facet
        
        OUTPUT: The restriction of ``comp_poly`` to the Weyl chambers smaller than ``facet``

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[3]
            sage: A2.y_facet(comp_poly, Facet)
            3*alpha[1] + 6*alpha[2]
        
        """
        if not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")

        R = self.root_space()

        return sum(
            (
                2 * self.numerical_invariant(comp_poly, facet, vertex)
                / self.euclidean_product(
                    self.weight_to_root_space(vertex),
                    self.weight_to_root_space(vertex)
                )
            ) * self.weight_to_root_space(vertex)
            for vertex in facet
        )

    def roots_facet_perp(self, facet):
        r"""Returns the roots orthogonal to every vertex of ``facet``, forming a root system in the perpendicular space to ``facet``.

        INPUT: ``facet`` -- a facet

        OUTPUT: A list of roots perpendicular to ``facet``

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: Facet = A2.facets()[3]
            sage: A2.roots_facet_perp(Facet)
            [alpha[1], -alpha[1]]

        """
        roots = self.root_space().roots()
        return [
            alpha for alpha in roots
            if all(
                self.euclidean_product(
                    alpha,
                    self.weight_to_root_space(vertex)
                ) == 0
                for vertex in facet
            )
        ]

    def ambient_vector_to_root_space(self, v):
        r"""Converts an ambient vector to an element of the root space. This method provides a common interface for vectors that may be represented in the root space, root lattice, or weight lattice.
        
        INPUT: ``v`` -- an element of the root space, root lattice, or weight lattice

        OUTPUT: An element of the root space representing ``v``

        EXAMPLES::
        
            sage: from cpoly import RootSystem, greater_than
        
            sage: A2 = RootSystem(['A', 2])
            sage: V = A2.root_space()
            sage: W = A2.weight_lattice()
            sage: alpha = V.simple_roots()
            sage: Lambda = W.fundamental_weights()
            sage: A2.ambient_vector_to_root_space(alpha[1])
            alpha[1]

            sage: A2.ambient_vector_to_root_space(Lambda[1])
            2/3*alpha[1] + 1/3*alpha[2]
        
        """
        R = self.root_space()

        if v.parent() is R:
            return v

        try:
            return R(v)
        except (TypeError, ValueError):
            return self.weight_to_root_space(v)


    def project_to_facet_perp(self, v, facet):
        r"""Orthogonally projects ``v`` onto the orthogonal complement of the linear span of ``facet``.

        INPUT: ``v`` -- an element of the root space, root lattice, or weight lattice

        OUTPUT: ``v`` -- an element of the root space, root lattice, or weight lattice

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: Lambda = W.fundamental_weights()
            sage: Facet = A2.facets()[3]
            sage: A2.project_to_facet_perp(Lambda[1], Facet)
            1/2*alpha[1]
        
        """
        R = self.root_space()
        v = self.ambient_vector_to_root_space(v)

        if not facet:
            return v

        facet_vectors = [
            self.ambient_vector_to_root_space(x)
            for x in facet
        ]

        #Gram matrix with respect to the root-system inner product
        G = matrix(QQ, [
            [
                self.euclidean_product(x, y)
                for y in facet_vectors
            ]
            for x in facet_vectors
        ])

        b = vector(QQ, [
            self.euclidean_product(x, v)
            for x in facet_vectors
        ])

        # G c = b
        coefficients = G.solve_right(b)

        #Projection onto span(facet)
        projection = R.zero()
        for i in range(len(facet_vectors)):
            projection += coefficients[i] * facet_vectors[i]

        #Orthogonal component.
        return R.from_vector(
            v.to_vector() - projection.to_vector()
        )

    def project_chamber_to_facet_perp(self, chamber, facet):
        r"""Orthogonally projects a Weyl chamber smaller than ``facet`` onto the perpendicular space to ``facet``.

        INPUT:
            - ``chamber`` -- a Weyl chamber
            - ``facet`` -- a facet in the Weyl chamber
        
        OUTPUT: The projected Weyl chamber in the perpendicular space of ``facet``
        
        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: Chamber = A2.weyl_chambers()[1]
            sage: Facet = A2.facets()[2]
            sage: A2.project_chamber_to_facet_perp(Chamber, Facet)
            (-1/2*alpha[1] - 1/2*alpha[2],)

        """
        if not greater_than(facet, chamber):
            raise ValueError("chamber must be smaller than facet")

        return tuple(
            self.project_to_facet_perp(vertex, facet)
            for vertex in chamber
            if vertex not in facet
        )

    def simple_roots_facet_perp(self, facet):
        r"""Returns the simple roots of the root system in :math:`P^\perp`.
        
        INPUT: ``facet`` -- a facet in the root system

        OUTPUT: Simple roots in the root system of the perpendicular space

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: Facet = A2.facets()[2]
            sage: A2.simple_roots_facet_perp(Facet)
            [-alpha[1] - alpha[2]]

        """
        roots = self.roots_facet_perp(facet)

        chambers = self.weyl_chambers_smaller_than_facet(facet)
        if not chambers:
            raise ValueError("facet has no smaller Weyl chambers")

        chamber = chambers[0]

        #A vector in the interior of the projected chamber
        h = sum(
            (
                self.project_to_facet_perp(
                    self.weight_to_root_space(vertex), facet
                )   
                for vertex in chamber
                if vertex not in facet
            ),
        self.root_space().zero()
        )

        positive_roots = [
            alpha for alpha in roots
            if self.euclidean_product(alpha, h) > 0
        ]

        return [
            alpha for alpha in positive_roots
            if not any(
                alpha - beta in positive_roots
                for beta in positive_roots
            )
        ]

    def cartan_matrix_facet_perp(self, facet):
        r"""Returns the Cartan matrix of the root system in :math:`P^\perp`.

        INPUT: ``facet`` -- a facet in the root system
        
        OUTPUT: The Cartan matrix in the root system of the perpendicular space
        
        EXAMPLE::
        
            sage: from cpoly import RootSystem, greater_than
        
            sage: A2 = RootSystem(['A', 2])
            sage: Facet = A2.facets()[2]
            sage: A2.cartan_matrix_facet_perp(Facet)
            [2]
        
        """
        simple_roots = self.simple_roots_facet_perp(facet)

        return matrix(QQ, [
            [
                2 * self.euclidean_product(alpha, beta)
                / self.euclidean_product(alpha, alpha)
                for beta in simple_roots
            ]
            for alpha in simple_roots
        ])

    def root_system_facet_perp(self, facet):
        r"""Returns the projected root system in :math:`P^\perp`.
        
        INPUT: ``facet`` -- a facet in the root system
                
        OUTPUT: The root system of the perpendicular space
                
        EXAMPLE::
                
            sage: from cpoly import RootSystem, greater_than
                
            sage: A2 = RootSystem(['A', 2])
            sage: Facet = A2.facets()[2]
            sage: A2.root_system_facet_perp(Facet)
            Root system of type ['A', 1] relabelled by {1: 0}
                
        """
        C = self.cartan_matrix_facet_perp(facet)
        CM = CartanMatrix(C)
        return RootSystem(CM.cartan_type())

    def ambient_to_facet_root_space(self, v, facet):
        r"""Expresses an ambient vector in the root space of the induced root system on :math:`P^\perp`, where :math:`P` is the span of ``facet``.

        INPUT:
            - ``v`` -- a vector in the ambient root space, assumed to lie in :math:`P^\perp``
            - ``facet`` -- a facet of the ambient Weyl chamber, represented as a tuple of vertices in the weight lattice

        OUTPUT: A vector in the root space of the root system induced on :math:`P^\perp`, having the same coordinates with respect to the induced simple roots as ``v`` has with respect to the corresponding ambient vectors. If :math:`P^\perp` is zero-dimensional, the induced root system has no simple roots and the zero vector of its root space is returned.

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: Lambda = W.fundamental_weights()
            sage: Facet = A2.facets()[2]
            sage: v = A2.project_to_facet_perp(Lambda[1], Facet)
            sage: A2.ambient_to_facet_root_space(v, Facet)
            -1/2*alpha[0]
        
        """
        self_P = self.root_system_facet_perp(facet)
        simple_roots = self.simple_roots_facet_perp(facet)

        if not simple_roots:
            return self_P.root_space().zero()

        G = matrix(QQ, [
            [
                self.euclidean_product(alpha, beta)
                for beta in simple_roots
            ]
            for alpha in simple_roots
        ])

        b = vector(QQ, [
            self.euclidean_product(alpha, v)
            for alpha in simple_roots
        ])

        coefficients = vector(QQ, G.solve_right(b))

        return self_P.root_space().from_vector(coefficients)

    def facet_chamber_correspondence(self, facet):
        r"""Returns the correspondence between Weyl chambers smaller than ``facet`` and Weyl chambers of the induced root system. The keys are chambers of the induced root system and the values are the corresponding ambient chambers.

        INPUT: ``facet`` -- a facet of the ambient root system, represented as a tuple of vertices in the weight lattice
    
        OUTPUT: A dictionary whose keys are Weyl chambers of the induced root system in :math:`P^\perp` and whose values are the corresponding Weyl chambers of the ambient root system containing ``facet`` as a face
    
        EXAMPLE::
    
            sage: from cpoly import RootSystem, greater_than
    
            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: Lambda = W.fundamental_weights()
            sage: facet = A2.facets()[2] 
            sage: A2.facet_chamber_correspondence(facet)
            {(-Lambda[0],): (-Lambda[1] + Lambda[2], Lambda[2]),
             (Lambda[0],): (-Lambda[1] + Lambda[2], -Lambda[1])}
        
        """
        self_P = self.root_system_facet_perp(facet)

        ambient_chambers = self.weyl_chambers_smaller_than_facet(facet)
        induced_chambers = self_P.weyl_chambers()

        simple_roots = self.simple_roots_facet_perp(facet)

        #Determine the signs of the ambient projected chambers
        ambient_signs = {}

        for chamber in ambient_chambers:
            projected = self.project_chamber_to_facet_perp(
                chamber, facet
            )

            h = sum(
                projected,
                self.root_space().zero()
            )

            ambient_signs[chamber] = tuple(
                self.euclidean_product(alpha, h) > 0
                for alpha in simple_roots
            )

        #Determine the signs of the induced chambers
        #We use the fundamental weights spanning each induced chamber
        #Their sum is an interior vector
        induced_signs = {}

        P = self_P.root_space()
        induced_simple_roots = list(
            P.simple_roots().values()
        )

        #We cannot use self_P.euclidean_product(), since self_P
        #is the ordinary Sage RootSystem. Compute the induced
        #Euclidean form directly
        A = matrix(QQ, self_P.cartan_matrix())
        D = diagonal_matrix(
            QQ,
            self_P.cartan_type().symmetrizer()
            )
        G = D * A

        for chamber in induced_chambers:
            h = sum(
                (
                    self_P.weight_to_root_space(weight)
                    for weight in chamber
                ),
                P.zero()
            )

            h_vec = h.to_vector()

            induced_signs[chamber] = tuple(
                (alpha.to_vector() * G * h_vec) > 0
                for alpha in induced_simple_roots
            )

        return {
            chamber_P: next(
                chamber
                for chamber in ambient_chambers
                if ambient_signs[chamber] == induced_signs[chamber_P]
            )
            for chamber_P in induced_chambers
        }

    def comp_poly_facet(self, comp_poly, facet):
        r"""Constructs the complementary polyhedron ``comp_poly`` projected down to :math:`P^\perp``, where :math:`P` is the ``facet``.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron of the ambient root system, represented as a dictionary whose keys are Weyl chambers and whose values are points in the ambient root space
            - ``facet`` -- a facet of the ambient Weyl chamber decomposition, represented as a tuple of vertices in the weight lattice

        OUTPUT:

            A dictionary representing the complementary polyhedron induced on ``P^\perp``

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[1]
            sage: A2.comp_poly_facet(comp_poly,Facet)
            {(-Lambda[0],): alpha[0], (Lambda[0],): 1/2*alpha[0]}

        """
        if not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")

        correspondence = self.facet_chamber_correspondence(facet)

        return {
            chamber_P: self.ambient_to_facet_root_space(
                self.project_to_facet_perp(
                    comp_poly[chamber],
                    facet
                ),
                facet
            )
            for chamber_P, chamber in correspondence.items()
        }

    def has_J1(self, comp_poly, facet, check=True):
        r"""Tests whether ``comp_poly`` has property ``J1`` with respect to ``facet``, which means all the numerical invariants are zero.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron, represented as a dictionary mapping Weyl chambers to elements of the root space
            - ``facet`` -- a facet of the ambient root system, represented as a tuple of vertices in the weight lattice
            - ``check`` -- checks if ``comp_poly`` is a complementary polyhedron or not

        OUTPUT: ``True`` if all numerical invariants associated to ``facet`` are zero, and ``False`` otherwise.

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[1]
            sage: A2.has_J1(comp_poly, Facet)
            True

        """
        if check and not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")

        return all(
            self.numerical_invariant(comp_poly, facet, vertex) == 0
            for vertex in facet
        )

    def has_C1(self, comp_poly, facet, check=True):
        r"""Tests whether ``comp_poly`` has property ``C1`` with respect to ``facet``, which means all the numerical invariants are greater than zero.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron, represented as a dictionary mapping Weyl chambers to elements of the root space
            - ``facet`` -- a facet of the ambient root system, represented as a tuple of vertices in the weight lattice
            - ``check`` -- checks if ``comp_poly`` is a complementary polyhedron or not

        OUTPUT: ``True`` if all numerical invariants associated to ``facet`` are zero, and ``False`` otherwise.

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[1]
            sage: A2.has_C1(comp_poly, Facet)
            False

        """
        if check and not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")

        return all(
            self.numerical_invariant(comp_poly, facet, vertex) > 0
            for vertex in facet
        )

    def has_J2(self, comp_poly, facet, check=True):
        r"""Tests whether ``comp_poly`` has property ``J2`` with respect to ``facet``, which meansns the complementary polyhedron induced on the orthogonal complement of ``facet`` is stable.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron of the ambient root system
            - ``facet`` -- a facet of the ambient root system, represented as a tuple of vertices in the weight lattice
            - ``check`` -- checks if ``comp_poly`` is a complementary polyhedron or not

        OUTPUT: ``True`` if the complementary polyhedron induced on :math:`P^\perp` is stable, and ``False`` otherwise

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[1]
            sage: A2.has_J2(comp_poly, Facet)
            False

        """
        if check and not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")

        comp_poly_P = self.comp_poly_facet(comp_poly, facet)
        self_P = self.root_system_facet_perp(facet)

        return self_P.is_stable(comp_poly_P, check=False)

    def has_C2(self, comp_poly, facet, check=True):
        r"""Tests whether ``comp_poly`` has property ``C2`` with respect to ``facet``, which means the complementary polyhedron induced on the orthogonal complement of ``facet`` is semistable.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron of the ambient root system
            - ``facet`` -- a facet of the ambient root system, represented as a tuple of vertices in the weight lattice
            - ``check`` -- checks if ``comp_poly`` is a complementary polyhedron or not

        OUTPUT: ``True`` if the complementary polyhedron induced on :math:`P^\perp` is stable, and ``False`` otherwise

        EXAMPLE::

            sage: from cpoly import RootSystem, greater_than

            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: Facet = A2.facets()[1]
            sage: A2.has_C2(comp_poly, Facet)
            False

        """
        if check and not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")

        comp_poly_P = self.comp_poly_facet(comp_poly, facet)
        self_P = self.root_system_facet_perp(facet)

        return self_P.is_semistable(comp_poly_P, check=False)

    def canonical_facets(self, comp_poly):
        r"""Returns the canonical facets satisfying both properties ``C1`` and ``C2``. A facet satisfies ``C1`` if all of its numerical invariants are greater than 0. It satisfies ``C2`` if the complementary polyhedron induced on the orthogonal complement of the facet is semistable.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron of the ambient root system
            - ``facet`` -- a facet of the ambient root system, represented as a tuple of vertices in the weight lattice

        OUTPUT: A list of facets satisfying both ``J1`` and ``J2``, if ``check`` is ``True``, we validate that ``comp_poly`` is a complementary polyhedron.

        EXAMPLE::
            sage: from cpoly import RootSystem, greater_than
            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): 1/2*alpha[2], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1]+alpha[2], (-Lambda[2], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[1] + Lambda[2], Lambda[2]): 1/2*(alpha[1]+alpha[2]), (Lambda[1], Lambda[1] - Lambda[2]): alpha[2], (-Lambda[2], -Lambda[1]): alpha[1]+alpha[2]}
            sage: A2.canonical_facets(comp_poly)
            (Lambda[2],)
        """
        if not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")
    
        for facet in self.facets():
            if self.has_C1(comp_poly, facet, check=False) and self.has_C2(comp_poly, facet, check=False):
                return facet
        return None

    def JH_facets(self, comp_poly):
        r"""Returns the Jordan-Hölder facets satisfying both properties ``J1`` and ``J2`` for a semistable complementary polyhedron. A facet satisfies ``J1`` if all of its numerical invariants vanish. It satisfies ``J2`` if the complementary polyhedron induced on the orthogonal complement of the facet is stable.

        INPUT:
            - ``comp_poly`` -- a complementary polyhedron of the ambient root system
            - ``facet`` -- a facet of the ambient root system, represented as a tuple of vertices in the weight lattice

        OUTPUT: A list of facets satisfying both ``J1`` and ``J2``.

        EXAMPLE::
            sage: from cpoly import RootSystem, greater_than
            sage: A2 = RootSystem(['A', 2])
            sage: W = A2.weight_lattice()
            sage: V = A2.root_space()
            sage: Lambda = W.fundamental_weights()
            sage: alpha = V.simple_roots()
            sage: comp_poly = {(Lambda[1], Lambda[2]): -alpha[1], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1], (-Lambda[2], Lambda[1] - Lambda[2]): -alpha[1], (-Lambda[1] + Lambda[2], Lambda[2]): alpha[1], (Lambda[1], Lambda[1] - Lambda[2]): -alpha[1], (-Lambda[2], -Lambda[1]): alpha[1]}
            sage: A2.JH_facets(comp_poly)
            [(-Lambda[2],), (Lambda[2],)]
        
        """
        if not self.is_comp_poly(comp_poly):
            raise ValueError("comp_poly must be a complementary polyhedron")

        if not self.is_semistable(comp_poly):
            raise ValueError("comp_poly must be semistable")

        facets = self.facets()

        # First find one J1/J2 facet and determine its dimension.
        target_size = None
        first_facet = None

        for facet in facets:
            if self.has_J1(comp_poly, facet, check=False) and self.has_J2(comp_poly, facet, check=False):
                first_facet = facet
                target_size = len(facet)
                break

        # No J1/J2 facets exist (though this never happens).
        if target_size is None:
            return []

        # Only search facets of the same dimension.
        return [
            facet
            for facet in facets
            if len(facet) == target_size
            and self.has_J1(comp_poly, facet, check=False)
            and self.has_J2(comp_poly, facet, check=False)
        ]

    def is_polystable(self, comp_poly, check=True):
            r"""Tests whether ``comp_poly`` is polystable. A complementary polyhedron is polystable if the degree of every one-dimensional facet is less than or equal to zero, and degree 0 facets have degree 0 oppositve facets
    
            INPUT:
                - ``comp_poly`` -- a complementary polyhedron
                - ``check`` -- checks if ``comp_poly`` is a complementary polyhedron or not
    
    
            OUTPUT: ``True`` if ``comp_poly`` is a semistable complementary polyhedron, and ``False`` otherwise
    
            EXAMPLE::
    
                sage: from cpoly import RootSystem, greater_than
    
                sage: A2 = RootSystem(['A', 2])
                sage: W = A2.weight_lattice()
                sage: V = A2.root_space()
                sage: Lambda = W.fundamental_weights()
                sage: alpha = V.simple_roots()
                sage: comp_poly = {(Lambda[1], Lambda[2]): -alpha[1], (-Lambda[1] + Lambda[2], -Lambda[1]): alpha[1], (-Lambda[2], Lambda[1] - Lambda[2]): -alpha[1], (-Lambda[1] + Lambda[2], Lambda[2]): alpha[1], (Lambda[1], Lambda[1] - Lambda[2]): -alpha[1], (-Lambda[2], -Lambda[1]): alpha[1]}
                sage: A2.is_polystable(comp_poly)
                True
            
            """
            if check and not self.is_comp_poly(comp_poly):
                raise ValueError("comp_poly must be a complementary polyhedron")

            facets = self.facets_one_dimensional()

            # Compute every degree only once.
            degrees = {
                facet: self.degree(comp_poly, facet)
                for facet in facets
            }

            for facet, degree in degrees.items():
                # Polystability implies semistability.
                if degree > 0:
                    return False

                # A degree-zero facet must have a degree-zero opposite facet.
                if degree == 0:
                    opposite = (-facet[0],)

                    if opposite not in degrees:
                        return False

                    if degrees[opposite] != 0:
                        return False
            return True