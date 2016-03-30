
import logging
from itertools import product

from .lpsolver import lp
from .formula import Atom, Radical, Formula

logger = logging.getLogger(__name__)


def _weight(element):
    if element in (Atom('N'), Atom('O'), Atom('P')):
        return 0.4
    elif isinstance(element, Radical):
        return 40.0
    return 1.0


def _weighted_formula(form):
    for e, mf in form.items():
        if e == Atom('H'):
            continue

        yield e, mf, _weight(e)


class UnbalancedReactionError(ValueError):
    """Raised when an unbalanced reaction is provided."""


def predict_rpair(reaction, compound_formula, solver):
    elements = set()
    for compound, value in reaction.compounds:
        if not compound.name in compound_formula:
            logger.warning('Reaction pairs are not available when formula'
                           ' is missing: {}, {}'.format(reaction, compound))
            return {}, {}

        f = compound_formula[compound.name]
        elements.update(e for e, _, _ in _weighted_formula(f))

    p = solver.create_problem()

    objective = 0
    for c1, c2 in product((c for c, _ in reaction.left),
                          (c for c, _ in reaction.right)):
        p.define(('m', c1, c2), lower=0)
        p.define(('q', c1, c2), lower=0)
        p.define(('omega', c1, c2), types=lp.VariableType.Binary)

        objective += (100 * p.var(('m', c1, c2)) +
                      90 * p.var(('q', c1, c2)) +
                      -0.1 * p.var(('omega', c1, c2)))

        f1 = compound_formula[c1.name]
        f2 = compound_formula[c2.name]

        p.define(*(('gamma', c1, c2, e) for e in elements),
                 types=lp.VariableType.Binary)
        p.define(*(('delta', c1, c2, e) for e in elements), lower=0)

        # Eq 12
        p.define(('x', c1, c2), types=lp.VariableType.Binary)
        p.define(('y', c1, c2), types=lp.VariableType.Binary)
        p.add_linear_constraints(
            p.var(('y', c1, c2)) <= 1 - p.var(('x', c1, c2)))

        # Eq 6, 9
        delta_wsum = p.expr({
            ('delta', c1, c2, e): _weight(e) for e in elements})
        p.add_linear_constraints(
            p.var(('m', c1, c2)) <= delta_wsum,
            p.var(('q', c1, c2)) <= delta_wsum)

        objective -= p.expr({('gamma', c1, c2, e): 1 for e in elements})

    p.set_objective(objective)

    # Eq 3
    zs = {}
    for c1, v in reaction.left:
        f = dict(compound_formula[c1.name].items())
        for e in elements:
            mf = f.get(e, 0)
            delta_sum = 0
            for c2, _ in reaction.right:
                delta_sum += p.var(('delta', c1, c2, e))
            try:
                p.add_linear_constraints(delta_sum == v * mf)
            except ValueError:
                raise UnbalancedReactionError('Unable to add constraint')

        # Eq 8, 11
        x_sum = p.expr({('x', c1, c2): 1 for c2, _ in reaction.right})
        y_sum = p.expr({('y', c1, c2): 1 for c2, _ in reaction.right})
        p.add_linear_constraints(x_sum <= 1, y_sum <= 1)

        # Eq 13
        zs[c1] = 0
        for e, mf, w in _weighted_formula(compound_formula[c1.name]):
            zs[c1] += w * mf

    # Eq 2
    for c2, v in reaction.right:
        f = dict(compound_formula[c2.name].items())
        for e in elements:
            mf = f.get(e, 0)
            delta_sum = 0
            for c1, _ in reaction.left:
                if e in compound_formula[c1.name]:
                    delta_sum += p.var(('delta', c1, c2, e))
            try:
                p.add_linear_constraints(delta_sum == v * mf)
            except ValueError:
                raise UnbalancedReactionError('Unable to add constraint')

    for c1, v1 in reaction.left:
        for c2, _ in reaction.right:
            f1 = dict(compound_formula[c1.name].items())
            f2 = dict(compound_formula[c2.name].items())
            for e in elements:
                mf = f1.get(e, 0)

                # Eq 4
                delta = p.var(('delta', c1, c2, e))
                omega = p.var(('omega', c1, c2))
                p.add_linear_constraints(delta <= float(v1) * mf * omega)

                # Eq 5
                if e in f2:
                    gamma = p.var(('gamma', c1, c2, e))
                    p.add_linear_constraints(delta <= float(v1) * mf * gamma)

            # Eq 7, 10
            m = p.var(('m', c1, c2))
            q = p.var(('q', c1, c2))
            x = p.var(('x', c1, c2))
            y = p.var(('y', c1, c2))
            p.add_linear_constraints(
                m <= float(v1) * zs[c1] * x, q <= float(v1) * zs[c1] * y)

    result = p.solve(lp.ObjectiveSense.Maximize)
    if not result:
        raise UnbalancedReactionError('Unable to solve')

    transfer = {}
    for c1, c2 in product((c for c, _ in reaction.left),
                          (c for c, _ in reaction.right)):
        f1 = compound_formula[c1.name]
        f2 = compound_formula[c2.name]

        items = {}
        for e in elements:
            v = result.get_value(('delta', c1, c2, e))
            if v % 1 == 0:
                v = int(v)
            if v != 0:
                items[e] = v

        if len(items) > 0:
            transfer[c1, c2] = Formula(items)

    return transfer
