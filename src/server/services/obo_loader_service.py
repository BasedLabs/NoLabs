import copy
import json

import networkx as nx
import obonet


def load_obo(ids):
    url = '../../media/go.obo'
    original_graph = obonet.read_obo(url)
    graph = copy.deepcopy(original_graph)

    nodes_to_remove = []
    new_edges = []
    for node in graph.nodes:
        if node in ids:
            continue  # 'GO:1900396'
        for in_edge in graph.in_edges(node):
            for out_edge in graph.out_edges(node):
                attributes = graph.get_edge_data(out_edge[0], out_edge[1])
                new_edges.append((in_edge[0], out_edge[1], attributes))
        nodes_to_remove.append(node)

    for (a, b, attributes) in new_edges:
        graph.add_edge(a, b, **{k: v for k, v in attributes.items() if isinstance(k, str)})
    for node in nodes_to_remove:
        graph.remove_node(node)

    data = nx.node_link_data(graph)

    # Dump data to JSON format
    graph_json = json.dumps(data, indent=4)

    data = nx.node_link_data(original_graph)
    original_graph_json = json.dumps(data, indent=4)
    return {
        'graph': graph_json,
        'original_graph_json': original_graph_json
    }