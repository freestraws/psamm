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

from ..command import Command, CommandError
from ..reaction import Compound
from .. import pathways

from six import iteritems, itervalues

logger = logging.getLogger(__name__)


class PathwaysCommand(Command):
    """Find shortest paths between two compounds."""

    name = 'pathways'
    title = 'Find shortest paths between two compounds'

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument('source', type=pathways.parse_compound,
                            help='Source compound')
        parser.add_argument('dest', type=pathways.parse_compound,
                            help='Destination compound')
        parser.add_argument('-n', type=int, default=5,
                            help='Number of pathways to find')
        super(PathwaysCommand, cls).init_parser(parser)

    def run(self):
        source = self._args.source
        dest = self._args.dest

        biomass_reaction = self._model.get_biomass_reaction()
        disconnect = set()
        if biomass_reaction is not None:
            disconnect.add(biomass_reaction)

        cost_func = pathways.FormulaCostFunction(self._model)
        # cost_func = pathways.AltFormulaCostFunction(self._model)
        # cost_func = pathways.ConnectivityCostFunction(self._mm)
        #connector = pathways.Connector(self._model, cost_func, disconnect)
        connector = pathways.RpairConnector(
            self._model, cost_func, disconnect)

        dests = {
            Compound('M_ala_DASH_L_c', 'c'),
            Compound('M_arg_DASH_L_c', 'c'),
            Compound('M_asn_DASH_L_c', 'c'),
            Compound('M_asp_DASH_L_c', 'c'),
            Compound('M_atp_c', 'c'),
            Compound('M_coa_c', 'c'),
            Compound('M_cys_DASH_L_c', 'c'),
            Compound('M_gln_DASH_L_c', 'c'),
            Compound('M_leu_DASH_L_c', 'c'),
            Compound('M_nad_c', 'c')
        }

        for dest in dests:
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

            with open('{}_{}.tsv'.format(source, dest), 'w') as f:
                for length, cost in stats:
                    f.write('{}\t{}\n'.format(length, cost))

            with open('{}_{}.dot'.format(source, dest), 'w') as f:
                self.write_graph(f, paths, source, dest)

    def write_graph(self, f, paths, source, dest):
        node_compounds = set()
        node_reactions = set()
        edge_score = {}

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
                        if edge not in edge_score:
                            edge_score[edge] = 0
                        edge_score[edge] += 1
                logger.info(
                    '- {} <- {}, cost={}, ({})'.format(
                        prev_compound, compound, cost,
                        ','.join(sorted(reactions))))
                prev_compound = compound

        max_score = max(itervalues(edge_score))

        f.write('digraph pathways {\n')
        for reaction in node_reactions:
            f.write('  "{}" [shape=box];\n'.format(reaction))
        for compound in (source, dest):
            f.write('  "{}" [color=red];\n'.format(compound))
        for edge, score in iteritems(edge_score):
            source, dest = edge
            width = (10.0 * (score - 1) / (max_score - 1)) + 1
            f.write('  "{}" -> "{}" [penwidth={}];\n'.format(
                source, dest, width))
        f.write('}\n')
