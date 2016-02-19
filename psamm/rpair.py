
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

    d1 = dict(f1.flattened().items())
    d2 = dict(f2.flattened().items())

    elements = set(element for element, _ in iteritems(d1))
    elements.update(element for element, _ in iteritems(d2))

    count, total = 0, 0
    w_count, w_total = 0, 0
    for element in elements:
        mi = min(d1.get(element, 0), d2.get(element, 0))
        mx = max(d1.get(element, 0), d2.get(element, 0))
        count += mi
        total += mx
        w_count += ATOM_WEIGHT.get(element, 1.0) * mi
        w_total += ATOM_WEIGHT.get(element, 1.0) * mx

    score = 0.0 if total == 0.0 else count / float(total)
    w_score = 0.0 if w_total == 0.0 else w_count / float(w_total)

    return score, w_score


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

    for _, value in reaction.compounds:
        if not isinstance(value, int):
            logger.warning('Reaction pairs are not available for non-integer'
                           ' reaction {}'.format(reaction))
            return {}, {}

    uninstantiated_left = dict(reaction.left)
    uninstantiated_right = dict(reaction.right)

    def compound_instances(uninstantiated):
        instances = []
        for compound, value in iteritems(uninstantiated):
            if value > 0:
                f = compound_formula[compound.name]
                instances.append(_CompoundInstance(compound, value, f))

        for inst in instances:
            uninstantiated[inst.compound] -= 1

        return instances

    def instantiate(uninstantiated, compound):
        n = uninstantiated[compound]
        if n > 0:
            f = compound_formula[compound.name]
            inst = _CompoundInstance(compound, n, f)
            uninstantiated[compound] -= 1
            return inst

        return None

    left = compound_instances(uninstantiated_left)
    right = compound_instances(uninstantiated_right)
    instances = left + right

    pairs = {}
    for inst1, inst2 in product(left, right):
        score, w_score = shared_elements(inst1.formula, inst2.formula)
        if score > 0.0:
            pairs[inst1, inst2] = w_score

    logger.debug('Reaction: {}'.format(reaction))

    transfer = {}
    while len(pairs) > 0:
        logger.debug('Pairs: {}'.format(pairs))

        (inst1, inst2), _ = max(iteritems(pairs), key=itemgetter(1))
        common = formula_common(inst1.formula, inst2.formula)

        logger.debug('Pair: {}, {}'.format(inst1, inst2))
        logger.debug('Common: {}'.format(common))

        key = (inst1.compound, inst1.n), (inst2.compound, inst2.n)
        transfer[key] = common
        for inst in (inst1, inst2):
            inst.formula = formula_remove(inst.formula, common)

        to_insert = set()

        inst = instantiate(uninstantiated_left, inst1.compound)
        if inst is not None:
            left.append(inst)
            instances.append(inst)
            to_insert.add(inst)

        inst = instantiate(uninstantiated_right, inst2.compound)
        if inst is not None:
            right.append(inst)
            instances.append(inst)
            to_insert.add(inst)

        to_update = {inst1, inst2}
        logger.debug('To update: {}'.format(to_update))

        to_delete = set()
        for inst1, inst2 in pairs:
            if inst1 in to_update or inst2 in to_update:
                if (len(dict(inst1.formula.items())) > 0 and
                        len(dict(inst2.formula.items())) > 0):
                    score, w_score = shared_elements(
                        inst1.formula, inst2.formula)
                    if score == 0.0:
                        to_delete.add((inst1, inst2))
                    else:
                        pairs[inst1, inst2] = w_score
                else:
                    to_delete.add((inst1, inst2))

        logger.debug('To delete: {}'.format(to_delete))
        for pair in to_delete:
            del pairs[pair]

        logger.debug('To insert: {}'.format(to_insert))
        for inst1, inst2 in product(left, right):
            if inst1 in to_insert or inst2 in to_insert:
                score, w_score = shared_elements(inst1.formula, inst2.formula)
                if score > 0.0:
                    pairs[inst1, inst2] = w_score

    balance = {}
    for inst in instances:
        if len(dict(inst.formula.items())) > 0:
            key = inst.compound, inst.n
            balance[key] = inst.formula

    return transfer, balance
