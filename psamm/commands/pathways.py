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

import logging
from itertools import product
from collections import defaultdict

from ..command import Command, MetabolicMixin, CommandError
from ..reaction import Compound, Direction
from .. import pathways
from ..datasource.reaction import parse_compound

from six import iteritems, itervalues

logger = logging.getLogger(__name__)

_REACTION_COLOR = '#ccebc5'
_COMPOUND_COLOR = '#b3cde3'
_ACTIVE_COLOR = '#fbb4ae'


class PathwaysCommand(MetabolicMixin, Command):
    """Find shortest paths between two compounds."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--source', type=parse_compound,
            action='append', default=None, help='Source compound')
        parser.add_argument(
            '--dest', type=parse_compound,
            action='append', default=None, help='Destination compound')
        parser.add_argument('-n', type=int, default=None,
                            help='Number of pathways to find')
        super(PathwaysCommand, cls).init_parser(parser)

    def run(self):
        biomass_reaction = self._model.get_biomass_reaction()

        if self._model.has_model_definition():
            subset = set(self._model.parse_model())
        else:
            subset = set(r.id for r in self._model.parse_reactions())

        if biomass_reaction is not None:
            subset.discard(biomass_reaction)

        #cost_func = pathways.FormulaCostFunction(self._model)
        cost_func = pathways.JaccardCostFunction(self._model)
        #cost_func = pathways.UniformCostFunction()
        # cost_func = pathways.AltFormulaCostFunction(self._model)
        # cost_func = pathways.ConnectivityCostFunction(self._mm)
        #connector = pathways.Connector(self._model, cost_func, disconnect)
        connector = pathways.RpairConnector(self._model, subset, cost_func)

        with open('reactions.tsv', 'w') as f:
            self.write_reaction_matrix(f, self._mm)

        with open('reaction_compounds.tsv', 'w') as f:
            self.write_reaction_compounds_matrix(f, self._mm)

        with open('connector_compounds.tsv', 'w') as f:
            self.write_connector_compounds_matrix(f, self._mm, connector)

        with open('connector.dot', 'w') as f:
            self.write_connector_graph(f, connector)

        with open('reactions.dot', 'w') as f:
            self.write_reaction_graph(f, self._mm, subset)

        if self._args.source is not None:
            sources = set(self._args.source)
        else:
            sources = set()
            for reaction_id in self._mm.reactions:
                if self._mm.is_exchange(reaction_id):
                    compound, _ = next(
                        self._mm.get_reaction_values(reaction_id))
                    sources.add(compound)

        if self._args.dest is not None:
            dests = set(self._args.dest)
        else:
            biomass = self._mm.get_reaction(biomass_reaction)
            dests = {c for c, _ in biomass.left}

        for source, dest in product(sources, dests):
            logger.info('Searching for {} -> {}...'.format(source, dest))
            paths = []
            stats = []
            for pathway, cost in pathways.find_pathways(
                    connector, source, dest):
                pathway = [(node.compound, node.reaction, node.g_score)
                           for node in reversed(pathway)]
                stats.append((len(pathway), cost))
                paths.append(pathway)
                if len(paths) == self._args.n:
                    break

            logger.info('Found {} paths'.format(len(paths)))

            if len(paths) > 0:
                with open('{}_{}.tsv'.format(source, dest), 'w') as f:
                    for length, cost in stats:
                        f.write('{}\t{}\n'.format(length, cost))

                with open('{}_{}.dot'.format(source, dest), 'w') as f:
                    self.write_graph(f, paths, source, dest)

    def write_reaction_matrix(self, f, model):
        compounds = sorted(model.compounds)
        reactions = sorted(model.reactions)

        connections = {}
        for reaction in model.reactions:
            rx = model.get_reaction(reaction)
            for compound, value in rx.compounds:
                connections.setdefault(compound, {})[reaction] = value

        f.write('\t'.join(str(r) for r in reactions) + '\n')
        for compound in compounds:
            values = [connections.get(compound, {}).get(r, 0)
                      for r in reactions]
            f.write('{}\t{}\n'.format(
                compound, '\t'.join(str(x) for x in values)))

    def write_reaction_compounds_matrix(self, f, model):
        compounds = sorted(model.compounds)
        compound_index = {c: i for i, c in enumerate(compounds)}

        connections = {}
        for reaction in model.reactions:
            rx = model.get_reaction(reaction)
            for (c1, _), (c2, _) in product(rx.left, rx.right):
                if rx.direction.forward:
                    connections.setdefault(c1, set()).add(c2)
                if rx.direction.reverse:
                    connections.setdefault(c2, set()).add(c1)

        f.write('\t'.join(str(c) for c in compounds) + '\n')
        for compound in compounds:
            values = (int(c in connections.get(compound, set()))
                      for c in compounds)
            f.write('{}\t{}\n'.format(
                compound, '\t'.join(str(x) for x in values)))

    def write_reaction_graph(self, f, model, subset):
        f.write('digraph pathways {\n')

        for compound in model.compounds:
            f.write('  "{}" [style=filled,fillcolor="#b3cde3"];\n'.format(
                compound))

        def dir_value(direction):
            if rx.direction == Direction.Forward:
                return 'forward'
            elif rx.direction == Direction.Reverse:
                return 'back'
            return 'both'

        for reaction in model.reactions:
            if reaction not in subset:
                continue
            rx = model.get_reaction(reaction)
            f.write('  "{}" [shape=box,style=filled,'
                    'fillcolor="#ccebc5"];\n'.format(reaction))
            for c, _ in rx.left:
                f.write('  "{}" -> "{}" [dir={}];\n'.format(
                    c, reaction, dir_value(rx.direction)))
            for c, _ in rx.right:
                f.write('  "{}" -> "{}" [dir={}];\n'.format(
                    reaction, c, dir_value(rx.direction)))

        f.write('}\n')

    def write_connector_compounds_matrix(self, f, model, connector):
        compounds = sorted(model.compounds)
        compound_index = {c: i for i, c in enumerate(compounds)}

        connections = {}
        for compound in compounds:
            for other, _ in connector.iter_all(compound):
                connections.setdefault(other, set()).add(compound)

        f.write('\t'.join(str(c) for c in compounds) + '\n')
        for compound in compounds:
            values = (int(c in connections.get(compound, set()))
                      for c in compounds)
            f.write('{}\t{}\n'.format(
                compound, '\t'.join(str(x) for x in values)))

    def write_connector_graph(self, f, connector):
        f.write('digraph pathways {\n')

        def ids(prefix):
            i = 0
            while True:
                yield '{}{}'.format(prefix, i)
                i += 1

        reaction_ids = ids('r')

        compound_set = set()
        compound_reaction = {}
        reaction_compound = {}
        reaction_labels = {}
        for compound in connector.compounds():
            compound_set.add(compound)
            for other, reactions in connector.iter_all(compound):
                compound_set.add(other)
                for (reaction, direction), cost in iteritems(reactions):
                    if cost is None:
                        continue

                    logger.info('Compounds: {}, {}, reaction: {}'.format(
                        compound, other, reaction))

                    key1 = compound, reaction
                    key2 = other, reaction
                    if (key1 not in compound_reaction and
                            key2 not in compound_reaction):
                        rid = next(reaction_ids)
                        compound_reaction[key1] = rid
                        compound_reaction[key2] = rid
                        reaction_compound[rid] = set([compound, other])
                        logger.info('Putting {}, {} in {}'.format(compound, other, rid))
                    elif (key1 in compound_reaction and
                            key2 in compound_reaction):
                        rid = compound_reaction[key1]
                        other_rid = compound_reaction[key2]
                        if rid != other_rid:
                            logger.info('Removing {}...'.format(other_rid))
                            for c in reaction_compound[other_rid]:
                                compound_reaction[c, reaction] = rid
                                reaction_compound[rid].add(c)
                                logger.info('Moving {} to {}'.format(c, rid))
                            del reaction_compound[other_rid]
                            del reaction_labels[other_rid]
                    elif key1 in compound_reaction:
                        rid = compound_reaction[key1]
                        compound_reaction[key2] = rid
                        reaction_compound[rid].add(other)
                        logger.info('Putting {} in {}'.format(other, rid))
                    elif key2 in compound_reaction:
                        rid = compound_reaction[key2]
                        compound_reaction[key1] = rid
                        reaction_compound[rid].add(compound)
                        logger.info('Putting {} in {}'.format(compound, rid))

                    reaction_labels[rid] = reaction

        inbound_reaction = {}
        outbound_reaction = {}
        for compound in connector.compounds():
            for other, reactions in connector.iter_all(compound):
                for (reaction, direction), cost in iteritems(reactions):
                    if cost is None:
                        continue

                    pen_width = 3 * (1 - cost)
                    reaction_id = compound_reaction[compound, reaction]

                    key = other, reaction_id
                    if key not in inbound_reaction:
                        inbound_reaction[key] = 0
                    inbound_reaction[key] += pen_width

                    key = reaction_id, compound
                    if key not in outbound_reaction:
                        outbound_reaction[key] = 0
                    outbound_reaction[key] += pen_width

        for (compound, reaction), pen_width in iteritems(inbound_reaction):
            f.write('  "{}" -> "{}" [penwidth={}];\n'.format(
                compound, reaction, pen_width))

        for (reaction, compound), pen_width in iteritems(outbound_reaction):
            f.write('  "{}" -> "{}" [penwidth={}];\n'.format(
                reaction, compound, pen_width))

        for compound in compound_set:
            f.write('  "{}" [style=filled,fillcolor="{}"];\n'.format(
                compound, _COMPOUND_COLOR))

        for reaction, label in iteritems(reaction_labels):
            f.write('  "{}" [shape=box,label="{}",style=filled,'
                    'fillcolor="{}"];\n'.format(
                        reaction, label, _REACTION_COLOR))

        f.write('}\n')

    def write_graph(self, f, paths, source, dest):
        if len(paths) == 0:
            return

        node_compounds = set()
        node_reactions = set()
        edges = set()
        edge_use_count = defaultdict(int)
        for pathway in paths:
            first_compound, _, cost = pathway[-1]
            logger.info('{}: {}'.format(first_compound, cost))
            prev_compound, _, _ = pathway[0]
            for compound, reaction, cost in pathway[1:]:
                node_compounds.add(compound)
                node_reactions.add(reaction)
                for edge in ((reaction[0], prev_compound),
                             (compound, reaction[0])):
                    edges.add(edge)
                    edge_use_count[edge] += 1
                logger.info(
                    '- {} <- {}, cost={}, ({})'.format(
                        prev_compound, compound, cost,
                        '{}[{}]'.format(reaction[0], reaction[1].symbol)))
                prev_compound = compound

        max_use_count = max(itervalues(edge_use_count))

        width, height = 50, 50
        margin = 0.5

        edge_props = defaultdict(dict)
        reaction_props = defaultdict(dict)
        compound_props = defaultdict(dict)

        for compound in node_compounds:
            compound_props[compound]['style'] = 'filled'
            compound_props[compound]['fillcolor'] = '"' + _COMPOUND_COLOR + '"'

        for reaction, direction in node_reactions:
            reaction_props[reaction]['shape'] = 'box'
            reaction_props[reaction]['style'] = 'filled'
            reaction_props[reaction]['fillcolor'] = '"' + _REACTION_COLOR + '"'

        for edge in edges:
            use_count = edge_use_count[edge]
            if max_use_count <= 1:
                pen_width = 1
            else:
                pen_width = (10.0 * (use_count - 1) / (max_use_count - 1)) + 1
            edge_props[edge]['penwidth'] = pen_width

        step = (height - 2*margin) / len(paths[0])
        prev_compound = None
        for i, (compound, reaction, _) in enumerate(paths[0]):
            compound_props[compound]['color'] = '"' + _ACTIVE_COLOR + '"'
            compound_props[compound]['pos'] = '"{},{}!"'.format(
                width / 2, margin + i*step)

            if reaction is not None:
                r, d = reaction
                edge_props[compound, r]['color'] = '"' + _ACTIVE_COLOR + '"'
                edge_props[r, prev_compound]['color'] = (
                    '"' + _ACTIVE_COLOR + '"')

            prev_compound = compound

        f.write('digraph pathways {\n')
        f.write('  size="{},{}!";\n'.format(width, height))
        f.write('  dpi=300;\n')

        for reaction, direction in node_reactions:
            if reaction in reaction_props:
                props = reaction_props[reaction]
                f.write('  "{}" [{}];\n'.format(
                    reaction, ','.join('{}={}'.format(k, v)
                                       for k, v in iteritems(props))))

        for compound in node_compounds:
            if compound in compound_props:
                props = compound_props[compound]
                f.write('  "{}" [{}];\n'.format(
                    compound, ','.join('{}={}'.format(k, v)
                                       for k, v in iteritems(props))))

        for (source, dest) in edges:
            props = edge_props.get((source, dest), {})
            if len(props) > 0:
                prop_string = ' [{}]'.format(
                    ','.join('{}={}'.format(k, v)
                             for k, v in iteritems(props)))
            else:
                prop_string = ''
            f.write('  "{}" -> "{}"{};\n'.format(
                source, dest, prop_string))

        f.write('}\n')
