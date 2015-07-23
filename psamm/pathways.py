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

from .reaction import Reaction, Compound
from .formula import Formula
from .heap import Heap

from six import iteritems


logger = logging.getLogger(__name__)


def parse_compound(s):
    m = re.match(r'(.*)\[(\w+)\]', s)
    if m:
        return Compound(m.group(1), compartment=m.group(2))
    return Compound(s, compartment=None)


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
        self._formulas = {}
        for compound in model.parse_compounds():
            self._formulas[compound.id] = Formula.parse(compound.formula)

    def _cost(self, score):
        if score == 0.0:
            return None
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


class AltFormulaCostFunction(object):
    def __init__(self, model):
        self._formulas = {}
        for compound in model.parse_compounds():
            self._formulas[compound.id] = Formula.parse(compound.formula)

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
    def __init__(self, model, cost_func, disconnect=None):
        self._connections = {}

        for reaction in model.reactions:
            if disconnect is not None and reaction in disconnect:
                continue

            rx = model.get_reaction(reaction)
            if (rx.direction == Reaction.Bidir or
                    rx.direction == Reaction.Right):
                for metabolite, _ in rx.right:
                    known_reactants = self._connections.setdefault(
                        metabolite, {})
                    for reactant, _ in rx.left:
                        if reactant not in known_reactants:
                            cost = cost_func.actual_cost(reactant, metabolite)
                            if cost is None:
                                continue
                            known_reactants[reactant] = cost, set([reaction])
                        else:
                            _, reaction_set = known_reactants[reactant]
                            reaction_set.add(reaction)
            if (rx.direction == Reaction.Bidir or
                    rx.direction == Reaction.Left):
                for metabolite, _ in rx.left:
                    known_reactants = self._connections.setdefault(
                        metabolite, {})
                    for reactant, _ in rx.right:
                        if reactant not in known_reactants:
                            cost = cost_func.actual_cost(reactant, metabolite)
                            if cost is None:
                                continue
                            known_reactants[reactant] = cost, set([reaction])
                        else:
                            _, reaction_set = known_reactants[reactant]
                            reaction_set.add(reaction)

    def get(self, compound):
        return iteritems(self._connections.get(compound, {}))


class CompoundNode(object):
    def __init__(self, compound, reactions, f_score, g_score, previous):
        self.compound = compound
        self.reactions = reactions
        self.f_score = f_score
        self.g_score = g_score
        self.previous = previous

    def pathway(self):
        node = self
        while node is not None:
            yield node
            node = node.previous

    def __repr__(self):
        return '<{} at {}, f={:.4f}, g={:.4f}>'.format(
            self.__class__.__name__, self.compound, self.f_score, self.g_score)


def find_pathways(connector, cost_func, source, dest):
    initial_cost = cost_func.admissible_cost(dest, source)
    initial_node = CompoundNode(dest, None, initial_cost, 0.0, None)

    compound_nodes = {dest: initial_node}

    def path_priority(node):
        return node.f_score

    paths_heap = Heap([initial_node], key=path_priority)

    while len(paths_heap) > 0:
        node = paths_heap.pop()

        logger.info(
            'Following reaction from {} to {} for new cost {} ({})...'.format(
                node.previous, node.compound, node.g_score, node.f_score))

        if node.compound == source:
            logger.info('At end of path: {}'.format(node))
            yield list(node.pathway()), node.g_score
            continue

        for next_compound, (next_cost, next_reactions) in connector.get(
                node.compound):
            if any(n.compound == next_compound for n in node.pathway()):
                logger.info('Skipping neighbor already in path {}'.format(
                    next_compound))
                continue

            g_score = node.g_score + next_cost

            try:
                node_cost = cost_func.admissible_cost(next_compound, source)
            except ValueError:
                continue
            logger.info('Allowing reaction from {} to {}'
                        ' with admissible cost {}, step_cost={}'.format(
                            node.compound, next_compound,
                            g_score + node_cost, next_cost))

            next_node = CompoundNode(
                next_compound, next_reactions, g_score + node_cost,
                g_score, node)
            paths_heap.push(next_node)
            if (next_compound not in compound_nodes or
                    compound_nodes[next_compound].g_score > g_score):
                compound_nodes[next_compound] = next_node
