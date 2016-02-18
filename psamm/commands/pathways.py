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
from ..reaction import Compound
from .. import pathways
from ..datasource.reaction import parse_compound

from six import iteritems, itervalues

logger = logging.getLogger(__name__)


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
        parser.add_argument('-n', type=int, default=5,
                            help='Number of pathways to find')
        super(PathwaysCommand, cls).init_parser(parser)

    def run(self):
        biomass_reaction = self._model.get_biomass_reaction()
        disconnect = set()
        if biomass_reaction is not None:
            disconnect.add(biomass_reaction)

        if self._model.has_model_definition():
            subset = set(self._model.parse_model())
        else:
            subset = None

        #cost_func = pathways.FormulaCostFunction(self._model)
        cost_func = pathways.JaccardCostFunction(self._model)
        #cost_func = pathways.UniformCostFunction()
        # cost_func = pathways.AltFormulaCostFunction(self._model)
        # cost_func = pathways.ConnectivityCostFunction(self._mm)
        #connector = pathways.Connector(self._model, cost_func, disconnect)
        connector = pathways.RpairConnector(
            self._model, subset, cost_func, disconnect)

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
                    connector, cost_func, source, dest):
                pathway = [(node.compound, node.reactions, node.g_score)
                           for node in reversed(pathway)]
                stats.append((len(pathway), cost))
                paths.append(pathway)
                if len(paths) == self._args.n:
                    break

            if len(paths) > 0:
                with open('{}_{}.tsv'.format(source, dest), 'w') as f:
                    for length, cost in stats:
                        f.write('{}\t{}\n'.format(length, cost))

                with open('{}_{}.dot'.format(source, dest), 'w') as f:
                    self.write_graph(f, paths, source, dest)

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
            for compound, reactions, cost in pathway[1:]:
                node_compounds.add(compound)
                for reaction in reactions:
                    node_reactions.add(reaction)
                    for edge in ((reaction, prev_compound),
                                 (compound, reaction)):
                        edges.add(edge)
                        edge_use_count[edge] += 1
                logger.info(
                    '- {} <- {}, cost={}, ({})'.format(
                        prev_compound, compound, cost,
                        ','.join(sorted(reactions))))
                prev_compound = compound

        max_use_count = max(itervalues(edge_use_count))

        width, height = 50, 50
        margin = 0.5

        edge_props = defaultdict(dict)
        reaction_props = defaultdict(dict)
        compound_props = defaultdict(dict)

        for reaction in node_reactions:
            reaction_props[reaction]['shape'] = 'box'

        for edge in edges:
            use_count = edge_use_count[edge]
            if max_use_count <= 1:
                width = 1
            else:
                width = (10.0 * (use_count - 1) / (max_use_count - 1)) + 1
            edge_props[edge]['penwidth'] = width

        step = (height - 2*margin) / len(paths[0])
        prev_compound = None
        for i, (compound, reactions, _) in enumerate(paths[0]):
            compound_props[compound]['color'] = 'red'
            compound_props[compound]['pos'] = '"{},{}!"'.format(
                width / 2, margin + i*step)

            if reactions is not None:
                for reaction in reactions:
                    reaction_props[reaction]['color'] = 'red'
                    edge_props[compound, reaction]['color'] = 'red'
                    edge_props[reaction, prev_compound]['color'] = 'red'

            prev_compound = compound

        f.write('digraph pathways {\n')
        f.write('  size="{},{}!";\n'.format(width, height))
        f.write('  dpi=300;\n')

        for reaction in node_reactions:
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
