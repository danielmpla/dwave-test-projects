import matplotlib.pyplot as plt

from TSPUtils import TSPUtils

import tsplib95

import dwave_networkx as dnx
import networkx as nx

from dimod.reference.samplers import ExactSolver
from dwave.system.samplers import DWaveSampler

from dwave.system.composites import EmbeddingComposite

from dwave_qbsolv import QBSolv


problem = tsplib95.load_problem('./burma14.tsp')
problemG = TSPUtils.make_complete(problem.get_graph())

sampler = ExactSolver()
d_sampler = EmbeddingComposite(DWaveSampler())

# TODO: visualize qubo
qubo = dnx.traveling_salesperson_qubo(problemG)

response = QBSolv().sample_qubo(qubo, num_reads=1)

# TODO: validate result
# TODO: visualize result

# TODO: do measurements
# TODO: visualize measurements

print("samples=" + str(list(response.samples())))
print("energies=" + str(list(response.data_vectors['energy'])))
print(response)

# d_answer = dnx.traveling_salesman(G, d_sampler)

# print(answer)
# print(answerP)
# print(d_answer)

# TODO: how we get the measurements?

# TODO: measure

# time for the compilation (dependent on local machine)
# time for Embedding
# time waiting for computation at dwave
# time running at dwave

# total computation