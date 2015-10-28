
import logging
from itertools import product
from operator import itemgetter

from .formula import Formula, Atom

from six import iteritems


logger = logging.getLogger(__name__)


ATOM_WEIGHT = {
    Atom('H'): 0.0,
}


def shared_elements(f1, f2):
    """Calculate score of shared elements."""

    d1, d2 = {}, {}
    for a, v in f1.flattened().items():
        d1[a] = ATOM_WEIGHT.get(a, 1.0) * v
    for a, v in f2.flattened().items():
        d2[a] = ATOM_WEIGHT.get(a, 1.0) * v

    elements = set(element for element, _ in iteritems(d1))
    elements.update(element for element, _ in iteritems(d2))

    count, total = 0, 0
    for element in elements:
        count += min(d1.get(element, 0), d2.get(element, 0))
        total += max(d1.get(element, 0), d2.get(element, 0))

    if total == 0.0:
        return 0.0
    return count / float(total)


def formula_common(f1, f2):
    """Return intersection of two formulas."""

    f1 = dict(f1.flattened().items())
    f2 = dict(f2.flattened().items())
    elements = set(element for element, _ in iteritems(f1))
    elements.update(element for element, _ in iteritems(f2))

    new = {}
    for element in elements:
        new[element] = min(f1.get(element, 0), f2.get(element, 0))

    return Formula(new)


def formula_remove(f1, f2):
    """Return a copy of formula f1 with the elements from f2 removed."""

    f1 = dict(f1.flattened().items())
    f2 = dict(f2.flattened().items())
    elements = set(element for element, _ in iteritems(f1))
    elements.update(element for element, _ in iteritems(f2))

    new = {}
    for element in elements:
        count = f1[element] - f2.get(element, 0)
        if count < 0:
            raise ValueError('Unable to remove {} from formula'.format(
                element))
        elif count > 0:
            new[element] = count

    return Formula(new)


class _CompoundInstance(object):
    def __init__(self, compound, n, formula):
        self.compound = compound
        self.n = n
        self.formula = formula

    def __repr__(self):
        return '<{} {} of {}>'.format(
            self.__class__.__name__, self.n, self.compound)

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                self.compound == other.compound and
                self.n == other.n)

    def __hash__(self):
        return hash(self.compound) ^ hash(self.n)


def predict_rpair(reaction, compound_formula):
    """Predict reaction pairs."""

    transfer = {}

    def compound_instances(compounds):
        for compound, value in compounds:
            f = compound_formula[compound.name]
            for i in range(value):
                yield _CompoundInstance(compound, i, f)

    left = list(compound_instances(reaction.left))
    right = list(compound_instances(reaction.right))
    instances = left + right

    pairs = {}
    for inst1, inst2 in product(left, right):
        pairs[inst1, inst2] = shared_elements(inst1.formula, inst2.formula)

    while len(pairs) > 0:
        (inst1, inst2), _ = max(iteritems(pairs), key=itemgetter(1))
        common = formula_common(inst1.formula, inst2.formula)

        logger.debug('Pair: {}, {}'.format(inst1, inst2))
        logger.debug('Common: {}'.format(common))

        key = (inst1.compound, inst1.n), (inst2.compound, inst2.n)
        transfer[key] = common
        for inst in (inst1, inst2):
            inst.formula = formula_remove(inst.formula, common)

        to_update = {inst1, inst2}
        to_delete = set()
        for inst1, inst2 in pairs:
            if inst1 in to_update or inst2 in to_update:
                if (len(dict(inst1.formula.items())) > 0 and
                        len(dict(inst2.formula.items())) > 0):
                    shared = shared_elements(inst1.formula, inst2.formula)
                    if shared == 0.0:
                        to_delete.add((inst1, inst2))
                    else:
                        pairs[inst1, inst2] = shared
                else:
                    to_delete.add((inst1, inst2))

        for pair in to_delete:
            del pairs[pair]

    balance = {}
    for inst in instances:
        if len(dict(inst.formula.items())) > 0:
            key = inst.compound, inst.n
            balance[key] = inst.formula

    return transfer, balance
