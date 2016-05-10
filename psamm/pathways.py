# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

import re
import math
import logging
from collections import Counter
from itertools import product, chain

from .reaction import Reaction, Compound, Direction
from .formula import Formula, Atom
from .heap import Heap
from . import rpair

from six import iteritems
from six.moves import zip_longest


logger = logging.getLogger(__name__)


def shared_elements(f1, f2):
    """Calculate score of shared elements"""

    f1 = dict(f1.flattened().items())
    f2 = dict(f2.flattened().items())
    elements = set(element for element, _ in iteritems(f1))
    elements.update(element for element, _ in iteritems(f2))

    count, total = 0, 0
    for element in elements:
        count += min(f1.get(element, 0), f2.get(element, 0))
        total += max(f1.get(element, 0), f2.get(element, 0))
    return count / float(total)


def model_compound_formulas(model):
    formulas = {}
    for compound in model.parse_compounds():
        if compound.formula is not None:
            formulas[compound.id] = Formula.parse(compound.formula)
    return formulas


def complete_pathway(pathways, pathway, count):
    final_compound, _, cost, _ = pathway[-1]
    all_pathways = pathways.setdefault(final_compound, [])
    if len(all_pathways) < count:
        all_pathways.append((cost, list(pathway)))
        worst_cost = None
    else:
        all_pathways.sort()
        worst_cost, _ = all_pathways[-1]
        if cost < worst_cost:
            all_pathways[:] = all_pathways[:-1] + [(cost, list(pathway))]

    return worst_cost


class FormulaCostFunction(object):
    def __init__(self, model):
        self._formulas = model_compound_formulas(model)

    def _cost(self, score):
        if score == 0.0:
            return float('inf')
        return -math.log(score)

    def actual_cost(self, source, dest):
        f1, f2 = self._formulas[source.name], self._formulas[dest.name]
        return self._cost(0.999 * shared_elements(f1, f2))

    def admissible_cost(self, source, dest):
        """Admissible estimation of the cost from source to dest

        The cost of going from source to dest must never be overestimated by
        this function.
        """
        if source == dest:
            return 0.0
        return self.actual_cost(source, dest)


class JaccardCostFunction(object):
    def __init__(self, model):
        self._formulas = model_compound_formulas(model)

    def _cost(self, score):
        return 1 - score

    def actual_cost(self, source, dest):
        f1, f2 = self._formulas[source.name], self._formulas[dest.name]
        return self._cost(0.999 * shared_elements(f1, f2))

    def admissible_cost(self, source, dest):
        if source == dest:
            return 0.0
        return self.actual_cost(source, dest)


class UniformCostFunction(object):
    def actual_cost(self, source, dest):
        if source == dest:
            return 0
        return 1

    def admissible_cost(self, source, dest):
        return self.actual_cost(source, dest)


class AltFormulaCostFunction(object):
    def __init__(self, model):
        self._formulas = model_compound_formulas(model)

    def actual_cost(self, source, dest):
        f1 = self._formulas[source.name]
        f2 = self._formulas[dest.name]

        f1 = dict(f1.flattened().items())
        f2 = dict(f2.flattened().items())
        elements = set(element for element, _ in iteritems(f1))
        elements.update(element for element, _ in iteritems(f2))

        count = 0
        for element in elements:
            m1 = min(f1.get(element, 0), f2.get(element, 0))
            m2 = max(f1.get(element, 0), f2.get(element, 0))
            count += m2 - m1
        return count

    def admissible_cost(self, source, dest):
        return self.actual_cost(source, dest)


class ConnectivityCostFunction(object):
    def __init__(self, model):
        self._connectivity = Counter()
        for reaction in model.reactions:
            for compound, _ in model.get_reaction_values(reaction):
                self._connectivity[compound] += 1

    def actual_cost(self, source, dest):
        return self._connectivity[dest]

    def admissible_cost(self, source, dest):
        if source == dest:
            return 0.0
        return self.actual_cost(source, dest)


class Connector(object):
    def __init__(self, model, subset):
        self._model = model
        self._subset = subset
        self._forward = {}
        self._reverse = {}
        self._compounds = set()
        self._cache_pairs()

    def _cache_pairs(self):
        for reaction in self._model.parse_reactions():
            if self._subset is not None and reaction.id not in self._subset:
                continue

            rx = reaction.equation
            if rx.direction.forward:
                for metabolite, _ in rx.right:
                    self._compounds.add(metabolite)
                    known_reactants = self._reverse.setdefault(
                        metabolite, {})
                    for reactant, _ in rx.left:
                        if self._filter_reaction_pairs(
                                reaction, reactant, metabolite):
                            continue

                        entry = self._cache_state(
                            reaction, Direction.Forward, reactant, metabolite)
                        if entry is not None:
                            self._compounds.add(reactant)
                            reactions = known_reactants.setdefault(
                                reactant, {})
                            reactions[reaction.id, Direction.Forward] = entry
                            self._forward.setdefault(
                                reactant, {}).setdefault(metabolite, {})[
                                    reaction.id, Direction.Forward] = entry
            if rx.direction.reverse:
                for metabolite, _ in rx.left:
                    self._compounds.add(metabolite)
                    known_reactants = self._reverse.setdefault(
                        metabolite, {})
                    for reactant, _ in rx.right:
                        if self._filter_reaction_pairs(
                                reaction, metabolite, reactant):
                            continue

                        entry = self._cache_state(
                            reaction, Direction.Reverse, reactant, metabolite)
                        if entry is not None:
                            self._compounds.add(reactant)
                            reactions = known_reactants.setdefault(
                                reactant, {})
                            reactions[reaction.id, Direction.Reverse] = entry
                            self._forward.setdefault(
                                reactant, {}).setdefault(metabolite, {})[
                                    reaction.id, Direction.Forward] = entry

    def compounds(self):
        return iter(self._compounds)

    def iter_all(self, compound):
        return chain(self.iter_all_forward(compound),
                     self.iter_all_reverse(compound))

    def iter_all_forward(self, compound):
        return iteritems(self._forward.get(compound, {}))

    def iter_all_reverse(self, compound):
        return iteritems(self._reverse.get(compound, {}))

    def has(self, c1, c2):
        return self.has_forward(c1, c2) or self.has_reverse(c1, c2)

    def has_forward(self, c1, c2):
        return c1 in self._forward and c2 in self._forward[c1]

    def has_reverse(self, c1, c2):
        return c1 in self._reverse and c2 in self._reverse[c1]

    def get_forward(self, c1, c2):
        return self._forward.get(c1, {}).get(c2, None)

    def get_reverse(self, c1, c2):
        return self._reverse.get(c1, {}).get(c2, None)

    def _filter_reaction_pairs(self, reaction, c_left, c_right):
        return False


class CostConnector(Connector):
    def __init__(self, model, subset, cost_func):
        self._cost_func = cost_func
        super(CostConnector, self).__init__(model, subset)

    def _cache_state(self, reaction, direction, reactant, metabolite):
        return self.actual_cost(reaction.id, reactant, metabolite)

    def admissible_cost(self, source, dest):
        return self._cost_func.admissible_cost(source, dest)

    def actual_cost(self, reaction, source, dest):
        return self._cost_func.actual_cost(source, dest)


class UndirectedConnector(Connector):
    """Connector wrapper that makes every edge bidirectional."""
    def __init__(self, connector):
        self._connector = connector

    def compounds(self):
        return self._connector.compounds()

    def iter_all(self, compound):
        return self._connector.iter_all(compound)

    iter_all_forward = iter_all
    iter_all_reverse = iter_all

    def has(self, c1, c2):
        return self._connector.has(c1, c2)

    has_forward = has
    has_reverse = has

    def _get(self, c1, c2):
        d = {}
        forward = self._connector.get_forward(c1, c2)
        if forward is not None:
            d.update(forward)

        reverse = self._connector.get_reverse(c1, c2)
        if reverse is not None:
            d.update(reverse)

        return d if len(d) > 0 else None

    get_forward = _get
    get_reverse = _get


class RpairConnector(CostConnector):
    def __init__(self, model, subset, cost_func):
        self._rpairs = {}
        formulas = model_compound_formulas(model)
        for reaction in model.parse_reactions():
            if subset is not None and reaction.id not in subset:
                continue
            transfer, _ = rpair.predict_rpair(
                reaction.equation, formulas)

            rpairs = {}
            for ((c1, _), (c2, _)), form in iteritems(transfer):
                form_dict = dict(form.items())
                if (form_dict.get(Atom('C'), 0) <= 1 and
                        (Atom('C') in formulas[c1.name] or
                         Atom('C') in formulas[c2.name])):
                    continue
                rpairs[c1, c2] = form

            self._rpairs[reaction.id] = rpairs
            logger.info('{}: {}'.format(reaction.id, rpairs))

        super(RpairConnector, self).__init__(model, subset, cost_func)

    def _filter_reaction_pairs(self, reaction, c_left, c_right):
        if (c_left, c_right) not in self._rpairs.get(reaction.id, {}):
            logger.info('{}: Filtering {} -> {}'.format(
                reaction.id, c_left, c_right))
            return True

        logger.info('{}: Allowing {} -> {}'.format(
            reaction.id, c_left, c_right))
        return False


class RpairConnectorCommon(RpairConnector):
    def __init__(self, model, subset, cost_func):
        self._formulas = model_compound_formulas(model)
        self._backbone = None
        super(RpairConnectorCommon, self).__init__(model, subset, cost_func)

    def set_source_dest(self, source, dest):
        self._backbone = rpair.formula_common(
            self._formulas[source.name], self._formulas[dest.name])
        self._cache_pairs()

    def actual_cost(self, reaction, source, dest):
        if self._backbone is not None:
            transfer = self._rpairs.get(reaction, {}).get((source, dest))
            if transfer is None:
                return None

            score, w_score = rpair.shared_elements(
                rpair.formula_common(transfer, self._backbone),
                self._backbone)

            if w_score < 0.3:
                logger.info(
                    '{}: Filtering {} -> {} because of low backbone'
                    ' similarity: {}'.format(reaction, source, dest, w_score))
                return None

        return super(RpairConnectorCommon, self).actual_cost(
            reaction, source, dest)


class CompoundNode(object):
    def __init__(self, compound, reaction, f_score, g_score, previous):
        self.compound = compound
        self.reaction = reaction
        self.f_score = f_score
        self.g_score = g_score
        self.previous = previous

    def pathway(self):
        node = self
        while node is not None:
            yield node
            node = node.previous

    def __repr__(self):
        reaction = None if self.reaction is None else '{}[{}]'.format(
            self.reaction[0], self.reaction[1].symbol)
        return '<{} at {}<-{}, f={:.4f}, g={:.4f}>'.format(
            self.__class__.__name__, self.compound, reaction,
            self.f_score, self.g_score)


def find_pathways(connector, source, dest):
    initial_cost = connector.admissible_cost(dest, source)
    initial_node = CompoundNode(dest, None, initial_cost, 0.0, None)

    compound_nodes = {dest: initial_node}

    def path_priority(node):
        return node.f_score

    paths_heap = Heap([initial_node], key=path_priority)

    while len(paths_heap) > 0:
        node = paths_heap.pop()

        logger.debug(
            'Following reaction from {} to {} for new cost {} ({})...'.format(
                node.previous, node.compound, node.g_score, node.f_score))

        if node.compound == source:
            logger.info('At end of path: {}'.format(node))
            yield list(node.pathway()), node.g_score
            continue

        neighbors = connector.iter_all_reverse(node.compound)
        for next_compound, next_reactions in neighbors:
            skip_compound = False
            for n in node.pathway():
                if n.compound == next_compound:
                    skip_compound = True
                    break

            if skip_compound:
                continue

            for (next_reaction, direction), next_cost in iteritems(
                    next_reactions):
                skip_reaction = False
                for n in node.pathway():
                    rev_rx = next_reaction, direction.flipped()
                    if n.reaction == rev_rx:
                        logger.debug(
                            'Skipping reaction that is already'
                            ' used in the reverse direction: {}'.format(
                                next_reaction))
                        skip_reaction = True
                        break

                if skip_reaction:
                    continue

                g_score = node.g_score + next_cost

                try:
                    node_cost = connector.admissible_cost(
                        next_compound, source)
                except ValueError:
                    continue
                logger.debug(
                    'Allowing reaction from {} to {} through {}'
                    ' with admissible cost {}, step_cost={}'.format(
                        node.compound, next_compound,
                        next_reaction, g_score + node_cost, next_cost))

                next_node = CompoundNode(
                    next_compound, (next_reaction, direction),
                    g_score + node_cost, g_score, node)
                paths_heap.push(next_node)
                if (next_compound not in compound_nodes or
                        compound_nodes[next_compound].g_score > g_score):
                    compound_nodes[next_compound] = next_node


def check_cost_function(connector):
    """Check that the triangle inequalite holds for valid triplets."""

    for c1 in connector.compounds():
        logger.info('Checking compound {}...'.format(c1))
        for c2, c1_reactions in connector.iter_all_reverse(c1):
            for c3, c2_reactions in connector.iter_all_reverse(c2):
                if connector.has_reverse(c1, c3):
                    for (r1, cost1), (r2, cost2), (r3, cost3) in product(
                            iteritems(c1_reactions),
                            iteritems(c2_reactions),
                            iteritems(connector.get_reverse(c1, c3))):
                        if cost1 + cost2 < cost3:
                            logger.warning(
                                'Invalid distance:'
                                ' {} -> {} ({}): {}, {} -> {} ({}): {},'
                                ' {} -> {} ({}): {}'.format(
                                    c1, c2, r1, cost1, c2, c3, r2, cost2,
                                    c1, c3, r3, cost3))

                # Must not overestimate
                ad_cost3 = connector.admissible_cost(c1, c3)
                for (r1, cost1), (r2, cost2) in product(
                        iteritems(c1_reactions),
                        iteritems(c2_reactions)):
                    if cost1 + cost2 < ad_cost3:
                        logger.warning(
                            'Invalid admissible cost:'
                            ' {} -> {} ({}): {}, {} -> {} ({}): {},'
                            ' {} -> {}: ({})'.format(
                                c1, c2, r1, cost1, c2, c3, r2, cost2,
                                c1, c3, ad_cost3))
