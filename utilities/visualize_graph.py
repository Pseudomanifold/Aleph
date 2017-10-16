#!/usr/bin/env python3

import networkx          as nx
import matplotlib.pyplot as plt
import sys

G         = nx.read_gml( sys.argv[1] )
weights   = nx.get_edge_attributes(G, 'value')
positions = nx.spring_layout(G,scale=2)

for weight in sorted(set( weights.values() ), reverse=True):
    edges = [ (u,v) for (u,v,d) in G.edges(data=True) if d['value'] >= weight ]

    plt.figure(figsize=[5,5])
    plt.axis('off')

    nx.draw_networkx_edges(G, positions, edge_color='white', arrows=False)
    nx.draw_networkx_edges(G, positions, edgelist=edges, edge_color='black', arrows=False)
    nx.draw_networkx_nodes(G, positions, node_size=50, node_color='gray')

    plt.savefig("/tmp/%02d.svg" % weight)
