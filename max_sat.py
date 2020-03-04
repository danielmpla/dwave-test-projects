import dimod
import matplotlib.pyplot as plt
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

quadratic = {(1,2): -1, (1,4): 1, (2,3): 10, (2,4): -20, (3,4): -20}
linear = {1: 0, 2: 0, 3: -1, 4: 30}


model = dimod.BinaryQuadraticModel(linear, quadratic, 0, dimod.BINARY)

s_sampler = dimod.SimulatedAnnealingSampler()
sampler = EmbeddingComposite(DWaveSampler())
# response = s_sampler.sample(model, num_reads=5)

response = sampler.sample(model, num_reads=1000, annealing_time=20, chain_strength=10, return_embedding=True)

embedding = response.info['embedding']
#
# plt.figure(figsize=(40,40))
# dnx.draw_chimera_embedding(dnx.chimera_graph(16, 16, 4), embedding)
# plt.savefig('embedding/embedding_max_sat.png')

min_e = min(list(response.data(['energy']))).energy

for datum in response.data(['sample', 'energy', 'num_occurrences', 'chain_break_fraction']):
    if datum.energy == min_e:
        print(datum.sample, datum.energy, "Occurrences: ", datum.num_occurrences, 'Chain Break: ' + str(100 * datum.chain_break_fraction))