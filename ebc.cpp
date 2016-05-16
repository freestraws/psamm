
#include <iostream>
#include <sstream>
#include <algorithm>
#include <array>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>


typedef std::pair<std::string, std::string> node_pair;

namespace std {
    template <>
    struct hash<node_pair> {
        size_t operator()(const node_pair& x) const throw() {
            std::hash<std::string> hasher;
            return 31 * hasher(x.first) + hasher(x.second);
        }
    };
}


struct graph_type {
    typedef std::string node_type;
    typedef std::unordered_multimap<node_type, node_type>::const_iterator
        edge_iterator;
    typedef std::pair<edge_iterator, edge_iterator> edge_range;

    std::unordered_set<node_type> nodes;
    std::unordered_set<std::pair<node_type, node_type>> edges;
    std::unordered_multimap<node_type, node_type> forward;

    void insert_edge(const node_type& source, const node_type& dest) {
        auto pair = std::make_pair(source, dest);
        auto got = edges.find(pair);
        if (got == edges.end()) {
            nodes.insert(source);
            nodes.insert(dest);
            edges.insert(pair);
            forward.insert(pair);
        }
    }

    void erase_edge(const node_type& source, const node_type& dest) {
        auto got = edges.find(std::make_pair(source, dest));
        if (got != edges.end()) {
            auto range = forward.equal_range(source);
            for (auto it = range.first; it != range.second; ++it) {
                if (it->second == dest) {
                    forward.erase(it);
                    break;
                }
            }

            edges.erase(got);
        }
    }

    size_t count_edges() const {
        return edges.size();
    }

    edge_range follow_edges(const node_type& source) const {
        return forward.equal_range(source);
    }
};


void print_graph(const graph_type& graph)
{
    for (auto it = graph.nodes.begin(); it != graph.nodes.end(); ++it) {
        std::cout << *it << std::endl;
        auto range = graph.follow_edges(*it);
        for (auto it2 = range.first; it2 != range.second; ++it2) {
            std::cout << "  " << it2->second << std::endl;
        }
    }
}


struct dijkstra_counts {
    int path_count;
    int dist;
};

struct dijkstra_result {
    std::unordered_multimap<std::string, std::string> prev_node;
    std::unordered_map<std::string, dijkstra_counts> counts;
};


void dijkstra_shortest(const graph_type& graph, const std::string& source,
                       dijkstra_result& result) {
    result.prev_node.clear();
    result.counts.clear();

    result.counts.insert(std::make_pair(source, dijkstra_counts{1, 0}));

    std::queue<std::string> open;
    std::unordered_set<std::string> closed;

    open.push(source);

    while (!open.empty()) {
        std::string current = open.front();
        open.pop();

        dijkstra_counts& counts = result.counts.at(current);

        closed.insert(current);

        auto range = graph.follow_edges(current);
        for (auto it = range.first; it != range.second; ++it) {
            const std::string& other = it->second;
            auto got_closed = closed.find(other);
            if (got_closed != closed.end()) continue;

            int dist = counts.dist + 1;
            int path_count = counts.path_count;

            auto got = result.counts.find(other);
            if (got == result.counts.end()) {
                open.push(other);
                result.counts.insert(std::make_pair(
                    other, dijkstra_counts{path_count, dist}));
                result.prev_node.insert(std::make_pair(other, current));
            } else {
                dijkstra_counts& other_counts = result.counts[other];
                if (other_counts.dist == dist) {
                    other_counts.path_count += path_count;
                    result.prev_node.insert(std::make_pair(other, current));
                }
            }
        }
    }
}


typedef std::unordered_map<node_pair, double> centrality_result;

struct node_dist_comparator {
    typedef std::pair<std::string, dijkstra_counts> node_pair;
    bool operator()(const node_pair& first,
                    const node_pair& second) const throw() {
        return first.second.dist < second.second.dist;
    }
};


void compound_pair_centrality(const graph_type& graph,
                              centrality_result& result)
{
    node_dist_comparator comparator;

    std::vector<std::string> nodes(graph.nodes.cbegin(), graph.nodes.cend());

    #pragma omp parallel
    {
        // Local result for the thread. These will be summed to the global
        // result at the end.
        std::unordered_map<node_pair, double> local_result;

        #pragma omp for schedule(dynamic)
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            const std::string& node = nodes[i];

            std::unordered_map<std::string, double> dependency;
            std::vector<std::pair<std::string, dijkstra_counts>> ordered_nodes;
            dijkstra_result dij_result;
            dijkstra_shortest(graph, node, dij_result);

            // Sort nodes in vector
            ordered_nodes.assign(
                dij_result.counts.begin(), dij_result.counts.end());
            std::sort(ordered_nodes.begin(), ordered_nodes.end(), comparator);

            for (auto it2 = ordered_nodes.rbegin();
                 it2 != ordered_nodes.rend(); ++it2) {
                const std::string& compound = it2->first;
                const dijkstra_counts& counts = it2->second;

                auto range = dij_result.prev_node.equal_range(compound);
                for (auto it3 = range.first; it3 != range.second; ++it3) {
                    const std::string& other = it3->second;
                    const dijkstra_counts& other_counts =
                        dij_result.counts.at(other);

                    double compound_dep = 0;
                    auto got_dep = dependency.find(compound);
                    if (got_dep != dependency.end()) {
                        compound_dep = got_dep->second;
                    }

                    double ratio = static_cast<double>(
                        other_counts.path_count) / counts.path_count;
                    double dep = ratio * (1 + compound_dep);

                    auto got = dependency.find(other);
                    if (got == dependency.end()) {
                        dependency.insert(std::make_pair(other, dep));
                    } else {
                        got->second += dep;
                    }

                    node_pair c_pair;
                    if (compound < other) {
                        c_pair = std::make_pair(compound, other);
                    } else {
                        c_pair = std::make_pair(other, compound);
                    }

                    // Add to local centrality result
                    {
                        auto got = local_result.find(c_pair);
                        if (got == local_result.end()) {
                            local_result.insert(std::make_pair(c_pair, dep));
                        } else {
                            got->second += dep;
                        }
                    }
                }
            }
        }

        // Add up thread local centrality to final result
        #pragma omp critical
        for (auto it = local_result.begin(); it != local_result.end(); ++it) {
            auto got = result.find(it->first);
            if (got == result.end()) {
                result.insert(*it);
            } else {
                got->second += it->second;
            }
        }
    }
}


struct centrality_comparator {
    typedef std::pair<node_pair, double> centrality_value;
    bool operator()(const centrality_value& first,
                    const centrality_value& second) const throw() {
        return first.second < second.second;
    }
};


typedef std::vector<node_pair> ebc_break_result;


void find_ebc_breaks(const graph_type& graph, int count,
                     ebc_break_result& result)
{
    graph_type g = graph;
    centrality_result cent_result;
    centrality_comparator comparator;
    std::unordered_set<std::string> nodes_cut;

    for (int i = 1; i <= count; ++i) {
        if (g.count_edges() == 0) break;

        cent_result.clear();
        compound_pair_centrality(g, cent_result);

        std::cerr << "Break " << i << std::endl;

        // Sort nodes in vector
        std::vector<std::pair<node_pair, double>> centrality_vec(
            cent_result.begin(), cent_result.end());
        std::sort(centrality_vec.begin(), centrality_vec.end(), comparator);

        auto max_it = centrality_vec.rbegin();
        double max_value = max_it->second;

        nodes_cut.clear();

        for (auto it = centrality_vec.rbegin();
             it != centrality_vec.rend(); ++it) {
            const node_pair& c_pair = it->first;

            if (it->second < max_value) {
                std::cerr << "No more breaks" << std::endl;
                break;
            }

            if (nodes_cut.find(c_pair.first) != nodes_cut.end()) {
                std::cerr << "Already broken " << c_pair.first <<
                    "; Skipping." << std::endl;
                continue;
            }

            if (nodes_cut.find(c_pair.second) != nodes_cut.end()) {
                std::cerr << "Already broken " << c_pair.second <<
                    "; Skipping." << std::endl;
                continue;
            }

            std::cerr << "Break at " << c_pair.first << ", " <<
                c_pair.second << ": " << it->second << std::endl;

            g.erase_edge(c_pair.first, c_pair.second);
            g.erase_edge(c_pair.second, c_pair.first);

            result.push_back(c_pair);

            nodes_cut.insert(c_pair.first);
            nodes_cut.insert(c_pair.second);
        }
    }
}


int main(int argc, char* argv[])
{
    graph_type graph;

    std::string line, value;
    std::istringstream line_values;
    while (std::getline(std::cin, line)) {
        std::string source;
        line_values.str(line);
        for (int i = 0; i < 2; i++) {
            std::getline(line_values, value, '\t');
            if (i == 0) {
                source = value;
            } else {
                graph.insert_edge(source, value);
                graph.insert_edge(value, source);
            }
        }
        line_values.clear();
    }

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " break-count" << std::endl;
        return 1;
    }

    int count = std::atoi(argv[1]);

    ebc_break_result result;
    find_ebc_breaks(graph, count, result);

    for (auto it = result.begin(); it != result.end(); ++it) {
        std::cout << it->first << "\t" << it->second << std::endl;
    }
}
