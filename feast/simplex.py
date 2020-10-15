"""Finite elements on simplices."""

import sympy
from itertools import product
from .symbolic import one, zero
from .finite_element import FiniteElement, make_integral_moment_dofs
from .polynomials import polynomial_set, Hcurl_polynomials, Hdiv_polynomials
from .functionals import (
    PointEvaluation,
    DotPointEvaluation,
    TangentIntegralMoment,
    NormalIntegralMoment,
    IntegralMoment,
)


class Lagrange(FiniteElement):
    """Lagrange finite element."""

    def __init__(self, reference, order):
        if order == 0:
            dofs = [
                PointEvaluation(
                    tuple(
                        sympy.Rational(1, reference.tdim + 1)
                        for i in range(reference.tdim)
                    ),
                )
            ]
        else:
            dofs = []
            for v in reference.reference_vertices:
                dofs.append(PointEvaluation(v))
            for edim in range(1, 4):
                for vs in reference.sub_entities(edim):
                    entity = reference.sub_entity_types[edim](
                        vertices=tuple(reference.reference_vertices[i] for i in vs)
                    )
                    for i in product(range(1, order), repeat=edim):
                        if sum(i) < order:
                            dofs.append(
                                PointEvaluation(
                                    tuple(
                                        o
                                        + sum(
                                            a[j] * b / order
                                            for a, b in zip(entity.axes, i)
                                        )
                                        for j, o in enumerate(entity.origin)
                                    )
                                )
                            )

        super().__init__(
            reference, polynomial_set(reference.tdim, 1, order), dofs, reference.tdim, 1
        )

    names = ["Lagrange", "P"]


class VectorLagrange(FiniteElement):
    """Vector Lagrange finite element."""

    def __init__(self, reference, order):
        scalar_space = Lagrange(reference, order)
        dofs = []
        if reference.tdim == 1:
            directions = [1]
        else:
            directions = [
                tuple(one if i == j else zero for j in range(reference.tdim))
                for i in range(reference.tdim)
            ]
        for p in scalar_space.dofs:
            for d in directions:
                dofs.append(DotPointEvaluation(p.point, d))

        super().__init__(
            reference,
            polynomial_set(reference.tdim, reference.tdim, order),
            dofs,
            reference.tdim,
            reference.tdim,
        )

    names = ["vector Lagrange", "vP"]


class NedelecFirstKind(FiniteElement):
    """Nedelec first kind Hcurl finite element."""

    def __init__(self, reference, order):
        poly = polynomial_set(reference.tdim, reference.tdim, order - 1)
        poly += Hcurl_polynomials(reference.tdim, reference.tdim, order)

        dofs = make_integral_moment_dofs(
            reference,
            edges=(TangentIntegralMoment, Lagrange, order - 1, 0),
            faces=(IntegralMoment, VectorLagrange, order - 2, 0),
            volumes=(IntegralMoment, VectorLagrange, order - 3, 0),
        )

        super().__init__(reference, poly, dofs, reference.tdim, reference.tdim)

    names = ["Nedelec", "Nedelec1", "N1curl"]


class RaviartThomas(FiniteElement):
    """Raviart-Thomas Hdiv finite element."""

    def __init__(self, reference, order):
        poly = polynomial_set(reference.tdim, reference.tdim, order - 1)
        poly += Hdiv_polynomials(reference.tdim, reference.tdim, order)

        dofs = make_integral_moment_dofs(
            reference,
            facets=(NormalIntegralMoment, Lagrange, order - 1, 0),
            cells=(IntegralMoment, VectorLagrange, order - 2, 0),
        )

        super().__init__(reference, poly, dofs, reference.tdim, reference.tdim)

    names = ["Raviart-Thomas", "RT", "N1div"]


class NedelecSecondKind(FiniteElement):
    """Nedelec second kind Hcurl finite element."""

    def __init__(self, reference, order):
        poly = polynomial_set(reference.tdim, reference.tdim, order)

        dofs = make_integral_moment_dofs(
            reference,
            edges=(TangentIntegralMoment, Lagrange, order, 0),
            faces=(IntegralMoment, RaviartThomas, order - 1, 1),
            volumes=(IntegralMoment, RaviartThomas, order - 2, 1),
        )

        super().__init__(reference, poly, dofs, reference.tdim, reference.tdim)

    names = ["Nedelec2", "N2curl"]


class BDM(FiniteElement):
    """Brezzi-Douglas-Marini Hdiv finite element."""

    def __init__(self, reference, order):
        poly = polynomial_set(reference.tdim, reference.tdim, order)

        dofs = make_integral_moment_dofs(
            reference,
            facets=(NormalIntegralMoment, Lagrange, order, 1),
            cells=(IntegralMoment, NedelecFirstKind, order - 1, 1),
        )

        super().__init__(reference, poly, dofs, reference.tdim, reference.tdim)

    names = ["Brezzi-Douglas-Marini", "BDM", "N2div"]
