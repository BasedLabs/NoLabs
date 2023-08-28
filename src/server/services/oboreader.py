import os

import obonet
import numbers


def read_obo():
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    f = os.path.join(__location__, 'go.obo')
    graph = obonet.read_obo(f)

    ids = {
        "GO:0009274",
        "GO:0000724",
        "GO:0009314",
        "GO:0009432",
        "GO:0005524",
        "GO:0006281",
        "GO:0006310",
        "GO:0005524",
        "GO:0006281",
        "GO:0006974",
        "GO:0000166"
    }

    def get_attributes(edge_data):
        for key in edge_data:
            if isinstance(key, numbers.Number):
                return get_attributes(edge_data[key])
            return [k for k in edge_data.keys()]

    # cleanup a graph from the nodes that we don't have in list
    nodes_to_remove = []
    for node in graph.nodes:
        if node in ids:
            continue
        for in_edge in graph.in_edges(node):
            for out_edge in graph.out_edges(node):
                if not graph.has_edge(in_edge[0], out_edge[1]):
                    attributes = graph.get_edge_data(out_edge[0], out_edge[1], 0)
                    if not attributes:
                        attributes = graph.get_edge_data(out_edge[0], out_edge[1])
                    graph.add_edge(in_edge[0], out_edge[1], **attributes)
        nodes_to_remove.append(node)

    for node in nodes_to_remove:
        graph.remove_node(node)

    shaped_graph = {} # {nodeId: {name: '', edges: { nodeId: {linkType}, }}}
    for node in graph.nodes:
        if node not in shaped_graph:
            shaped_graph[node] = {'name': graph.nodes[node]['name'], 'namespace': graph.nodes[node]['namespace'],
                                  'edges': {}}
        for out_edge in graph.out_edges(node):
            out_node = out_edge[1]
            edge_data = get_attributes(graph.get_edge_data(node, out_edge[1]))
            shaped_graph[node]['edges'][out_node] = edge_data

    return shaped_graph


if __name__ == "__main__":
    obo = read_obo()