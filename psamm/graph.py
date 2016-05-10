
from itertools import count

from six import iteritems


def _graphviz_prop_string(d):
    return ','.join('{}="{}"'.format(k, v) for k, v in iteritems(d))


class Entity(object):
    """Base class for graph entities."""
    def __init__(self, props={}):
        self._props = dict(props)

    @property
    def props(self):
        return self._props


class Graph(Entity):
    """Graph entity representing a collection of nodes and edges."""
    def __init__(self, props={}):
        super(Graph, self).__init__(props)
        self._nodes = set()
        self._edges = set()

    def add_node(self, node):
        self._nodes.add(node)

    def add_edge(self, edge):
        if edge.source not in self._nodes or edge.dest not in self._nodes:
            raise ValueError('Edge nodes not in graph')

        self._edges.add(edge)

    def write_graphviz(self, f):
        f.write('digraph {\n')

        for k, v in iteritems(self.props):
            f.write(' {}="{}";\n'.format(k, v))

        next_id = count(0)
        for node in self._nodes:
            if 'id' not in node.props:
                node.props['id'] = 'n{}'.format(next(next_id))

            f.write(' "{}"[{}]\n'.format(
                node.props['id'], _graphviz_prop_string(node.props)))

        for edge in self._edges:
            f.write(' "{}" -> "{}"[{}]\n'.format(
                edge.source.props['id'], edge.dest.props['id'],
                _graphviz_prop_string(edge.props)))

        f.write('}\n')


class Node(Entity):
    """Node entity represents a vertex in the graph."""
    def __init__(self, props={}):
        super(Node, self).__init__(props)


class Edge(Entity):
    """Edge entity represents a connection between nodes."""
    def __init__(self, source, dest, props={}):
        super(Edge, self).__init__(props)
        self.source = source
        self.dest = dest
