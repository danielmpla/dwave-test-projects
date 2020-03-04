from pprint import pprint

import matplotlib.pyplot as plt

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

import networkx as nx
import dwave_networkx as dnx

from dimod import BinaryQuadraticModel

from pyqubo import Binary

import minorminer


def get_max_chain_length(embedded):
    max_chain_length = 0
    for _, chain in embedded.items():
        if len(chain) > max_chain_length:
            max_chain_length = len(chain)
    return max_chain_length


def read_integers(filename):
    with open(filename) as f:
        return [int(elem) for elem in f.read().split()]


def test_embeddings(graph: nx.Graph, tries: int):
    shortest_chain = {}
    longest_chain = {}

    shortest = 100
    longest = 0

    for i in range(0, tries):
        embedded_graph = minorminer.find_embedding(G.edges(), dnx.chimera_graph(16, 16, 4))
        chain = get_max_chain_length(embedded_graph)

        if chain < shortest:
            shortest = chain
            shortest_chain = embedded_graph
        if chain > longest:
            longest = chain
            longest_chain = embedded_graph

    pprint((shortest, longest))

    plt.figure(figsize=(40, 40))
    dnx.draw_chimera_embedding(dnx.chimera_graph(16, 16, 4), shortest_chain)
    plt.savefig('embedding/max_cut_60.8_short.png')

    plt.figure(figsize=(40, 40))
    dnx.draw_chimera_embedding(dnx.chimera_graph(16, 16, 4), longest_chain)
    plt.savefig('embedding/max_cut_60.8_long.png')

    return shortest_chain


file_it = iter(read_integers('data/g05_60.8'))
# Number of vertices
n = next(file_it)
# Number of edges
m = next(file_it)

# Origin of each edge
origin = [None]*m
# Destination of each edge
dest = [None]*m
# Weight of each edge
w = [None]*m

weighted_edges = []

for e in range(m):
    origin[e] = next(file_it)
    dest[e] = next(file_it)
    w[e] = next(file_it)

    weighted_edges.append((origin[e], dest[e], float(w[e])))

G = nx.Graph()

G.add_weighted_edges_from(weighted_edges)

h = {v: 0. for v in G}
J = {(u, v): 1 for u, v in G.edges}

sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_ising(h, J, chain_strength=1, num_reads=1000, annealing_time=20, return_embedding=True)

embedding = response.info['embedding']

plt.figure(figsize=(40, 40))
dnx.draw_chimera_embedding(dnx.chimera_graph(16, 16, 4), embedding)
plt.savefig('embedding/embedding_g05_60.8_sampler.png')

pprint(response.info)
pprint(get_max_chain_length(embedding))

occ = max(list(response.data(['num_occurrences']))).num_occurrences
min_e = min(list(response.data(['energy']))).energy

for datum in response.data(['sample', 'energy', 'num_occurrences']):
    if datum.energy == min_e:
        print(datum.sample, datum.energy, "Occurrences: ", datum.num_occurrences)
