from pprint import pprint

import matplotlib.pyplot as plt

from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

import dwave_networkx as dnx

from pyqubo import Binary

product = 42

x0, x1, x2, y0, y1, y2 = Binary('x0'), Binary('x1'), Binary('x2'), Binary('y0'), Binary('y1'), Binary('y2')
H = (product - (x0 + 2*x1 + 4*x2)*(y0 + 2*y1 + 4*y2))**2

m = H.compile((product**2) * 100).to_dimod_bqm()

max_val = 1.0

for k, v in m.linear.items():
    if v > max_val:
        max_val = v

sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample(m, chain_strength=max_val, num_reads=1000, annealing_time=20, return_embedding=True)

embedding = response.info['embedding']

plt.figure(figsize=(40,40))
dnx.draw_chimera_embedding(dnx.chimera_graph(16, 16, 4), embedding)
plt.savefig('embedding/embedding15.png')

pprint(response.info)

occ = max(list(response.data(['num_occurrences']))).num_occurrences
min_e = min(list(response.data(['energy']))).energy

for datum in response.data(['sample', 'energy', 'num_occurrences']):
    if datum.energy == min_e:
        print(datum.sample, datum.energy, "Occurrences: ", datum.num_occurrences)
        print('x=' + str(datum.sample['x0']) + str(datum.sample['x1']) + str(datum.sample['x2']) + ' y=' + str(datum.sample['y0']) + str(datum.sample['y1']) + str(datum.sample['y2']))
        print('x=' + str(datum.sample['x0'] + 2*datum.sample['x1'] + 4*datum.sample['x2']) + ' y=' + str(datum.sample['y0'] + 2*datum.sample['y1'] + 4*datum.sample['y2']))
