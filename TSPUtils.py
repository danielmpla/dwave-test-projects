import networkx as nx

class TSPUtils():
    def __init__(self):
        pass

    @classmethod
    def make_complete(cls, graph: nx.Graph):
        mapping = {}
        edges_to_remove = []

        for i in graph.nodes():
            mapping[i] = i - 1
            edges_to_remove.append((i - 1, i - 1))
        
        graph = nx.relabel_nodes(graph, mapping)
        graph.remove_edges_from(edges_to_remove)

        return graph

