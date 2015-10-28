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
from .. import pathways

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

        paths = []
        for pathway, cost in pathways.find_pathways(
                connector, cost_func, source, dest):
            pathway = [(node.compound, node.reactions, node.g_score)
                       for node in reversed(pathway)]
            paths.append(pathway)
            if len(paths) == self._args.n:
                break

        for pathway in paths:
            first_compound, _, cost = pathway[-1]
            logger.info('{}: {}'.format(first_compound, cost))
            prev_compound, _, _ = pathway[0]
            for compound, reactions, cost in pathway[1:]:
                logger.info(
                    '- {} <- {}, cost={}, ({})'.format(
                        prev_compound, compound, cost,
                        ','.join(sorted(reactions))))
                prev_compound = compound
