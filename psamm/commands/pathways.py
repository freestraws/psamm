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

import math
import csv
import logging
from itertools import product, count
from collections import defaultdict, Counter

from ..command import Command, MetabolicMixin, CommandError
from ..reaction import Compound, Direction
from .. import pathways, graph
from ..datasource.reaction import parse_compound
from ..heap import Heap

from six import iteritems, itervalues

logger = logging.getLogger(__name__)

_REACTION_COLOR = '#ccebc5'
_COMPOUND_COLOR = '#b3cde3'
_ACTIVE_COLOR = '#fbb4ae'
_ALT_COLOR = '#ccb460'

_INVIS_NODE = {'shape': 'point', 'width': 0}


def path_use_values(paths):
    values = {}
    for pathway in paths:
        prev_compound, _, _ = pathway[0]
        for compound, reaction, cost in pathway[1:]:
            for edge in ((reaction[0], prev_compound),
                         (compound, reaction[0])):
                if edge not in values:
                    values[edge] = 0
                values[edge] += 1


def dijkstra_shortest(connector, breaks, source):
    """Dijkstra's shortest paths from source."""

    prev_node = {source: []}
    dist = {source: 0}
    path_count = Counter({source: 1})
    open_nodes = Heap([source], key=lambda x: dist[x])
    closed_nodes = set()

    while len(open_nodes) > 0:
        current = open_nodes.pop()
        if current not in dist:
            break

        closed_nodes.add(current)

        for other, reactions in connector.iter_all_forward(current):
            cpair = tuple(sorted([current, other]))
            if other in closed_nodes:
                continue

            reaction_set = set()
            for (reaction, direction), cost in iteritems(reactions):
                if cost is None or (reaction, cpair) in breaks:
                    continue
                reaction_set.add((reaction, cpair))

            if len(reaction_set) == 0:
                continue

            alt_dist = dist[current] + 1
            if other not in dist or alt_dist < dist[other]:
                dist[other] = alt_dist
                if other not in open_nodes:
                    open_nodes.push(other)
                else:
                    open_nodes.update(other)
                path_count[other] = path_count[current] * len(reaction_set)
                prev_node[other] = [(current, reaction_set)]
            elif other in dist and alt_dist == dist[other]:
                path_count[other] += path_count[current] * len(reaction_set)
                prev_node[other].append((current, reaction_set))

    return dist, prev_node, path_count


def reaction_centrality(connector, breaks):
    """Calculate reaction centrality."""
    centrality = Counter()
    for initial in connector.compounds_forward():
        dependency = Counter()
        reaction_dependency = Counter()
        dist, prev_node, path_count = dijkstra_shortest(
            connector, breaks, initial)

        for compound, d in sorted(
                iteritems(dist), key=lambda x: x[1], reverse=True):
            for other, reaction_set in prev_node[compound]:
                ratio = path_count[other] / float(path_count[compound])
                dep = ratio * (1 + dependency[compound])
                dependency[other] += dep

                for reaction in reaction_set:
                    reaction_dependency[reaction] += dep

        for reaction, value in iteritems(reaction_dependency):
            centrality[reaction] += value

    return centrality


def compound_centrality(connector, breaks):
    """Calculate compound centrality."""
    centrality = Counter()
    for initial in connector.compounds_forward():
        dependency = Counter()
        dist, prev_node, path_count = dijkstra_shortest(
            connector, breaks, initial)

        for compound, d in sorted(
                iteritems(dist), key=lambda x: x[1], reverse=True):
            for other, _ in prev_node[compound]:
                ratio = path_count[other] / float(path_count[compound])
                dependency[other] += ratio * (1 + dependency[compound])

            if compound != initial:
                centrality[compound] += dependency[compound]

    return centrality


def find_ebc_breaks(connector, n=10):
    breaks = set()
    break_order = []

    for i in range(n):
        new_breaks = set()
        break_compounds = set()
        c = sorted(iteritems(reaction_centrality(connector, breaks)),
                   key=lambda x: x[1], reverse=True)
        max_value = c[0][1]
        for (reaction, cpair), value in c:
            if value < max_value:
                logger.info('No more breaks: {}, {}'.format(
                    reaction, value))
                break
            c1, c2 = cpair
            if c1 in break_compounds or c2 in break_compounds:
                logger.info('Skipping {} because compound'
                            ' already broken'.format(reaction))
                continue
            logger.info('Break at {} ({}<->{}), {}'.format(
                reaction, c1, c2, value))
            new_breaks.add((reaction, cpair))
            breaks.add((reaction, cpair))
            break_compounds.add(c1)
            break_compounds.add(c2)

        with open('ebc-{}.tsv'.format(i), 'w') as f:
            for (reaction, cpair), value in c:
                f.write('{}\t{}\t{}\t{}\n'.format(
                    reaction, cpair[0], cpair[1], value))

        break_order.append(new_breaks)

    logger.info('Breaks: {}'.format(break_order))
    return breaks


def tarjan_components(connector, breaks):
    """Find strongly connected components."""
    next_index = count()
    compound_index = {}
    compound_lowlink = {}

    stack = []
    open_compounds = set()
    components = []

    def strong_connect(compound):
        compound_index[compound] = next(next_index)
        compound_lowlink[compound] = compound_index[compound]
        logger.info('Strong connect: {}, {}'.format(
            compound, compound_index[compound]))

        stack.append(compound)
        open_compounds.add(compound)

        for other, reactions in connector.iter_all_forward(compound):
            cpair = tuple(sorted([compound, other]))
            reaction_set = set()
            for (reaction, direction), cost in iteritems(reactions):
                if cost is None or (reaction, cpair) in breaks:
                    continue
                reaction_set.add((reaction, cpair))

            if len(reaction_set) == 0:
                continue

            if other not in compound_index:
                strong_connect(other)
                compound_lowlink[compound] = min(
                    compound_lowlink[compound], compound_lowlink[other])
                logger.info('Update lowlink: {}={}'.format(
                    compound, compound_lowlink[compound]))
            elif other in open_compounds:
                compound_lowlink[compound] = min(
                    compound_lowlink[compound], compound_index[other])
                logger.info('Update lowlink: {}={}'.format(
                    compound, compound_lowlink[compound]))

        if compound_index[compound] == compound_lowlink[compound]:
            component = []
            other = None
            while other != compound:
                other = stack.pop()
                open_compounds.remove(other)
                component.append(other)

            logger.info('New SC component: {}: {}'.format(compound, component))

            if len(component) > 2:
                components.append(component)

        logger.info('Done with {}'.format(compound))

    for compound in connector.compounds_forward():
        if compound not in compound_index:
            strong_connect(compound)

    return components


def calculate_compound_degree(connector):
    """Calculate the in-degree and out-degree for each compound."""
    indegree = Counter()
    outdegree = Counter()

    for compound in connector.compounds_forward():
        for other, reactions in connector.iter_all_forward(compound):
            outdegree[compound] += 1
            indegree[other] += 1

    return indegree, outdegree


def calculate_clustering_coefficients(connector):
    """Calculate the local clustering for each compounds.

    Return neighborhood size and raw edge saturation count.
    """
    values = {}
    for compound in connector.compounds_forward():
        # Neighborhood is the set of compounds with out-edges to compound or
        # in-edges from compound.
        neighborhood = set(
            other for other, _ in connector.iter_all_forward(compound))
        neighborhood.update(
            other for other, _ in connector.iter_all(compound))

        n_size = len(neighborhood)
        if n_size < 2:
            continue

        count = 0
        for c1, c2 in product(neighborhood, neighborhood):
            if c1 == c2:
                continue

            if connector.has_forward(c1, c2):
                count += 1

        values[compound] = n_size, count

    return values


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
        parser.add_argument(
            '--edge-values', type=str, default=None, help='Values for edges')
        parser.add_argument(
            '--clusters', type=str, default=None, help='Reaction clusters')
        parser.add_argument(
            '--breaks', type=int, default=0, help='EBC breaks')
        super(PathwaysCommand, cls).init_parser(parser)

    def run(self):
        biomass_reaction = self._model.get_biomass_reaction()
        subset = set(self._mm.reactions)

        if biomass_reaction is not None:
            subset.discard(biomass_reaction)

        edge_values = None
        if self._args.edge_values is not None:
            raw_values = {}
            with open(self._args.edge_values, 'r') as f:
                for row in csv.reader(f, delimiter='\t'):
                    raw_values[row[0]] = float(row[1])

            edge_values = {}
            for reaction in self._mm.reactions:
                rx = self._mm.get_reaction(reaction)
                if reaction in raw_values:
                    flux = raw_values[reaction]
                    if abs(flux) < 1e-9:
                        continue

                    if flux > 0:
                        for compound, value in rx.right:
                            edge_values[reaction, compound] = (
                                flux * float(value))
                        for compound, value in rx.left:
                            edge_values[compound, reaction] = (
                                flux * float(value))
                    else:
                        for compound, value in rx.left:
                            edge_values[reaction, compound] = (
                                -flux * float(value))
                        for compound, value in rx.right:
                            edge_values[compound, reaction] = (
                                -flux * float(value))

        clusters = {}
        if self._args.clusters is not None:
            with open(self._args.clusters, 'r') as f:
                f.readline()  # skip header
                for row in csv.reader(f, delimiter='\t'):
                    reaction_id = row[0]
                    cluster_id = row[1]
                    clusters[reaction_id] = cluster_id

        #cost_func = pathways.FormulaCostFunction(self._model)
        cost_func = pathways.JaccardCostFunction(self._model)
        #cost_func = pathways.UniformCostFunction()
        # cost_func = pathways.AltFormulaCostFunction(self._model)
        # cost_func = pathways.ConnectivityCostFunction(self._mm)
        #connector = pathways.Connector(self._model, cost_func, disconnect)
        connector = pathways.RpairConnector(self._model, subset, cost_func)

        breaks = find_ebc_breaks(connector, self._args.breaks)

        components = tarjan_components(connector, breaks)
        for component in components:
            logger.info('{}'.format(component))
        logger.info('{} SC components'.format(len(components)))

        indegree, outdegree = calculate_compound_degree(connector)
        degree = indegree + outdegree
        for compound, count in degree.most_common(50):
            logger.info('{}: {}'.format(compound, count))

        def connectivity_score(n_size, count):
            max_edges = n_size * (n_size-1)
            return float(max_edges - count) / (max_edges + count)

        coeffs = calculate_clustering_coefficients(connector)
        for compound, value in sorted(
                iteritems(coeffs), key=lambda x: connectivity_score(*x[1])):
            n_size, count = value
            logger.info('CC: {}: {}, {}'.format(
                compound, value, connectivity_score(n_size, count)))

        centrality = reaction_centrality(connector, breaks)

        with open('reactions.tsv', 'w') as f:
            self.write_reaction_matrix(f, self._mm)

        with open('reaction_compounds.tsv', 'w') as f:
            self.write_reaction_compounds_matrix(f, self._mm)

        with open('connector_compounds.tsv', 'w') as f:
            self.write_connector_compounds_matrix(f, self._mm, connector)

        def extended_compound_label(compound):
            label = str(compound)
            if compound in coeffs:
                label += '\n{:.2f}'.format(
                    connectivity_score(*coeffs[compound]))
            return label

        def compound_label(compound):
            return str(compound)

        def extended_reaction_label(reaction, pairs):
            label = str(reaction)
            for pair in pairs:
                spair = tuple(sorted(pair))
                if (reaction, spair) in centrality:
                    label += '\n{}<=>{}: {:.2f}'.format(
                        pair[0], pair[1], centrality[reaction, spair])
            return label

        def reaction_label(reaction, pairs):
            return str(reaction)

        with open('connector.dot', 'w') as f:
            self.write_connector_graph(
                f, connector, self._mm, biomass_reaction, edge_values,
                clusters, breaks, compound_label, reaction_label)

        with open('reactions.dot', 'w') as f:
            self.write_reaction_graph(
                f, self._mm, self._mm.reactions, edge_values)

        for r, v in iteritems(centrality):
            logger.info('{}: {}'.format(r, v))
        return

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
                if edge_values is None:
                    edge_values = path_use_values(paths)

                with open('{}_{}.tsv'.format(source, dest), 'w') as f:
                    for length, cost in stats:
                        f.write('{}\t{}\n'.format(length, cost))

                with open('{}_{}.dot'.format(source, dest), 'w') as f:
                    self.write_graph(f, paths, source, dest, edge_values)

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

    def write_reaction_graph(self, f, model, subset, edge_values):
        g = graph.Graph()

        compound_nodes = {}
        for compound in model.compounds:
            node = graph.Node({
                'style': 'filled',
                'fillcolor': _COMPOUND_COLOR})
            g.add_node(node)
            compound_nodes[compound] = node

        min_edge_value = math.log(min(itervalues(edge_values)))
        max_edge_value = math.log(max(itervalues(edge_values)))
        edge_value_span = max_edge_value - min_edge_value

        def dir_value(direction):
            if rx.direction == Direction.Forward:
                return 'forward'
            elif rx.direction == Direction.Reverse:
                return 'back'
            return 'both'

        for reaction in model.reactions:
            if reaction not in subset:
                logger.info('Skip {}'.format(reaction))
                continue
            rx = model.get_reaction(reaction)

            color = _REACTION_COLOR
            if model.is_exchange(reaction):
                color = _ACTIVE_COLOR

            node = graph.Node({
                'shape': 'box',
                'style': 'filled',
                'fillcolor': color
            })
            g.add_node(node)

            edge_props = {'dir': dir_value(rx.direction)}

            def pen_width(value):
                return (9 * (math.log(value) - min_edge_value) /
                        edge_value_span) + 1

            def final_props(edge1, edge2):
                if edge_values is not None:
                    p = {}
                    if edge1 in edge_values:
                        p['dir'] = 'forward'
                        p['penwidth'] = pen_width(edge_values[edge1])
                    elif edge2 in edge_values:
                        p['dir'] = 'back'
                        p['penwidth'] = pen_width(edge_values[edge2])
                    else:
                        p['style'] = 'dotted'
                        p['dir'] = dir_value(rx.direction)
                else:
                    p = props
                return p

            for c, _ in rx.left:
                edge1 = c, reaction
                edge2 = reaction, c
                g.add_edge(graph.Edge(
                    compound_nodes[c], node, final_props(edge1, edge2)))

            for c, _ in rx.right:
                edge1 = reaction, c
                edge2 = c, reaction
                g.add_edge(graph.Edge(
                    node, compound_nodes[c], final_props(edge1, edge2)))

        g.write_graphviz(f)

    def write_connector_compounds_matrix(self, f, model, connector):
        compounds = sorted(model.compounds)
        compound_index = {c: i for i, c in enumerate(compounds)}

        connections = {}
        for compound in compounds:
            for other, _ in connector.iter_all_forward(compound):
                connections.setdefault(other, set()).add(compound)

        f.write('\t'.join(str(c) for c in compounds) + '\n')
        for compound in compounds:
            values = (int(c in connections.get(compound, set()))
                      for c in compounds)
            f.write('{}\t{}\n'.format(
                compound, '\t'.join(str(x) for x in values)))

    def write_connector_graph(
            self, f, connector, model, biomass, edge_values, clusters,
            breaks, compound_label=None, reaction_label=None):
        g = graph.Graph()

        if compound_label is None:
            compound_label = str

        if reaction_label is None:
            reaction_label = lambda x: str(x[0])

        def ids(prefix):
            i = 0
            while True:
                yield '{}{}'.format(prefix, i)
                i += 1

        reaction_ids = ids('r')
        compound_ids = ids('c')

        # Split reaction nodes into nodes that do not share compounds.
        compound_reaction = {}
        reaction_name = {}
        reaction_compound = {}
        reaction_pairs = {}
        reaction_props = {}
        split_reaction = set()
        for compound in connector.compounds_forward():
            for other, reactions in connector.iter_all_forward(compound):
                for (reaction, direction), cost in iteritems(reactions):
                    if cost is None:
                        continue

                    logger.info('Compounds: {}, {}, reaction: {}'.format(
                        compound, other, reaction))

                    key1 = compound, reaction
                    key2 = other, reaction
                    cpair = tuple(sorted([compound, other]))
                    if (reaction, cpair) in breaks:
                        if key1 not in compound_reaction:
                            rid = next(reaction_ids)
                            compound_reaction[key1] = rid
                            reaction_compound[rid] = {compound}
                            reaction_name[rid] = reaction
                            split_reaction.add(rid)
                        else:
                            split_reaction.add(compound_reaction[key1])

                        if key2 not in compound_reaction:
                            rid = next(reaction_ids)
                            compound_reaction[key2] = rid
                            reaction_compound[rid] = {other}
                            reaction_name[rid] = reaction
                            split_reaction.add(rid)
                        else:
                            split_reaction.add(compound_reaction[key2])

                        logger.info('Keeping {}, {} separate'.format(
                            compound, other))
                    elif (key1 not in compound_reaction and
                            key2 not in compound_reaction):
                        rid = next(reaction_ids)
                        compound_reaction[key1] = rid
                        compound_reaction[key2] = rid
                        reaction_compound[rid] = {compound, other}
                        reaction_name[rid] = reaction
                        reaction_pairs[rid] = {(compound, other)}
                        logger.info('Putting {}, {} in {}'.format(
                            compound, other, rid))
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

                            # Add new pair and transfer from other reaction
                            reaction_pairs.setdefault(rid, set()).add(
                                (compound, other))
                            reaction_pairs[rid].update(
                                reaction_pairs.get(other_rid, set()))

                            del reaction_name[other_rid]
                            del reaction_compound[other_rid]
                            reaction_pairs.pop(other_rid, None)
                    elif key1 in compound_reaction:
                        rid = compound_reaction[key1]
                        compound_reaction[key2] = rid
                        reaction_compound[rid].add(other)
                        reaction_name[rid] = reaction
                        reaction_pairs.setdefault(rid, set()).add(
                            (compound, other))
                        logger.info('Putting {} in {}'.format(other, rid))
                    elif key2 in compound_reaction:
                        rid = compound_reaction[key2]
                        compound_reaction[key1] = rid
                        reaction_compound[rid].add(compound)
                        reaction_name[rid] = reaction
                        reaction_pairs.setdefault(rid, set()).add(
                            (compound, other))
                        logger.info('Putting {} in {}'.format(compound, rid))

        for rid, reaction in iteritems(reaction_name):
            color = _REACTION_COLOR
            if rid in split_reaction:
                color = _ALT_COLOR

            logger.info('{}, {}: {}'.format(
                rid, reaction, reaction_pairs.get(rid, set())))

            label = reaction_label(reaction, reaction_pairs.get(rid, set()))

            reaction_props[rid] = {
                'label': label,
                'fillcolor': color
            }

        min_edge_value = math.log(min(itervalues(edge_values)))
        max_edge_value = math.log(max(itervalues(edge_values)))
        edge_value_span = max_edge_value - min_edge_value

        # Create dicts of inbound/outbound edges for reactions.
        compound_props = {}
        compound_id_map = {}
        inbound_reaction = {}
        outbound_reaction = {}
        for compound in connector.compounds_forward():
            compound_connected = False
            for other, reactions in connector.iter_all_forward(compound):
                for (reaction, direction), cost in iteritems(reactions):
                    if cost is None:
                        continue

                    cluster = clusters.get(reaction)

                    # Create compound nodes for cluster
                    if (compound, cluster) not in compound_id_map:
                        logger.info('Creating {}, {}'.format(compound, cluster))
                        compound_id = next(compound_ids)
                        compound_id_map[compound, cluster] = compound_id
                        compound_props[compound_id] = {
                            'label': compound_label(compound),
                            'fillcolor': _COMPOUND_COLOR
                        }
                    else:
                        compound_id = compound_id_map[compound, cluster]

                    if (other, cluster) not in compound_id_map:
                        logger.info('Creating {}, {}'.format(other, cluster))
                        other_id = next(compound_ids)
                        compound_id_map[other, cluster] = other_id
                        compound_props[other_id] = {
                            'label': compound_label(other),
                            'fillcolor': _COMPOUND_COLOR
                        }
                    else:
                        other_id = compound_id_map[other, cluster]

                    # Create inbound/outbound connections
                    key = compound_reaction[other, reaction], other_id
                    if key not in outbound_reaction:
                        outbound_reaction[key] = 0
                    outbound_reaction[key] += edge_values.get(
                        (reaction, other), 0)

                    key = compound_id, compound_reaction[compound, reaction]
                    if key not in inbound_reaction:
                        inbound_reaction[key] = 0
                    inbound_reaction[key] += edge_values.get(
                        (compound, reaction), 0)

        # Add exchange reaction nodes to graph
        for reaction in model.reactions:
            if not model.is_exchange(reaction):
                continue

            cluster = clusters.get(reaction)
            rx = model.get_reaction(reaction)
            for c, _ in rx.compounds:
                if (c, cluster) in compound_id_map:
                    reaction_id = next(reaction_ids)
                    reaction_name[reaction_id] = reaction
                    reaction_props[reaction_id] = {
                        'label': reaction_label(reaction, set()),
                        'fillcolor': _ACTIVE_COLOR
                    }

                    compound_id = compound_id_map[c, cluster]
                    inbound_reaction[compound_id, reaction_id] = (
                        edge_values.get((c, reaction), 0))
                    outbound_reaction[reaction_id, compound_id] = (
                        edge_values.get((reaction, c), 0))

        # Add biomass reaction node to graph
        if biomass is not None:
            cluster = clusters.get(biomass)
            rx = model.get_reaction(biomass)
            for c, _ in rx.left:
                if (c, cluster) in compound_id_map:
                    reaction_id = next(reaction_ids)
                    reaction_name[reaction_id] = biomass
                    reaction_props[reaction_id] = {
                        'label': reaction_label(biomass, set()),
                        'fillcolor': _ACTIVE_COLOR
                    }

                    compound_id = compound_id_map[c, cluster]
                    inbound_reaction[compound_id, reaction_id] = (
                        edge_values.get((c, biomass), 0))

        def edge_props(flux):
            props = {}
            if flux == 0:
                props['style'] = 'dotted'
            else:
                props['penwidth'] = (
                    9 * ((math.log(flux) - min_edge_value) /
                         edge_value_span) + 1)
            return props

        compound_nodes = {}
        for compound, props in iteritems(compound_props):
            node = graph.Node({'style': 'filled'})
            node.props.update(props)
            g.add_node(node)
            compound_nodes[compound] = node

        reaction_nodes = {}
        for reaction, props in iteritems(reaction_props):
            node = graph.Node({
                'shape': 'box',
                'style': 'filled',
            })
            node.props.update(props)
            g.add_node(node)

            input_node = graph.Node(_INVIS_NODE)
            g.add_node(input_node)

            output_node = graph.Node(_INVIS_NODE)
            g.add_node(output_node)

            g.add_edge(graph.Edge(input_node, node))
            g.add_edge(graph.Edge(node, output_node))

            reaction_nodes[reaction] = input_node, output_node

            if reaction in reaction_compound:
                for c, v in model.get_reaction_values(reaction_name[reaction]):
                    if c not in reaction_compound[reaction]:
                        c_node = graph.Node({
                            'xlabel': str(c), 'fontsize': 10,
                            'shape': 'point', 'width': 0, 'margin': 0})
                        g.add_node(c_node)
                        if v < 0:
                            g.add_edge(graph.Edge(c_node, input_node))
                        else:
                            g.add_edge(graph.Edge(output_node, c_node))

        for (compound, reaction), flux in iteritems(inbound_reaction):
            g.add_edge(graph.Edge(
                compound_nodes[compound], reaction_nodes[reaction][0],
                edge_props(flux)))

        for (reaction, compound), flux in iteritems(outbound_reaction):
            g.add_edge(graph.Edge(
                reaction_nodes[reaction][1], compound_nodes[compound],
                edge_props(flux)))

        g.write_graphviz(f)

    def write_graph(self, f, paths, source, dest, edge_values):
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
            compound_props[compound]['fillcolor'] = _COMPOUND_COLOR

        for reaction, direction in node_reactions:
            reaction_props[reaction]['shape'] = 'box'
            reaction_props[reaction]['style'] = 'filled'
            reaction_props[reaction]['fillcolor'] = _REACTION_COLOR

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
            compound_props[compound]['color'] = _ACTIVE_COLOR
            compound_props[compound]['pos'] = '{},{}!'.format(
                width / 2, margin + i*step)

            if reaction is not None:
                r, d = reaction
                edge_props[compound, r]['color'] = _ACTIVE_COLOR
                edge_props[r, prev_compound]['color'] = _ACTIVE_COLOR

            prev_compound = compound

        f.write('digraph pathways {\n')
        f.write('  size="{},{}!";\n'.format(width, height))
        f.write('  dpi=300;\n')

        for reaction, direction in node_reactions:
            if reaction in reaction_props:
                props = reaction_props[reaction]
                f.write('  "{}" [{}];\n'.format(reaction, gv_props(props)))

        for compound in node_compounds:
            if compound in compound_props:
                props = compound_props[compound]
                f.write('  "{}" [{}];\n'.format(compound, gv_props(props)))

        for (source, dest) in edges:
            props = edge_props.get((source, dest), {})
            if len(props) > 0:
                prop_string = ' [{}]'.format(gv_props(props))
            else:
                prop_string = ''
            f.write('  "{}" -> "{}"{};\n'.format(
                source, dest, prop_string))

        f.write('}\n')
