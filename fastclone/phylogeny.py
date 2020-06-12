"""Subclonal phylogeny inference"""

import networkx as nx
import random
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

def infer(subclones, score, output):
    G = nx.Graph()

    if len(subclones) > 1:
        root = 0
        edge_ls = [(root, 1)]

        if len(subclones) > 2:
            subclone_snv_number = Counter([np.nanargmax(row) for row in subclones])
            i = 2
            while i != len(subclones):
                if subclone_snv_number[i] < subclone_snv_number[i - 1]:
                    edge_ls.append((root, i))
                else:
                    edge_ls.append((i - 1, i))
                    root = i - 1
                i += 1

        G.add_edges_from(edge_ls)

    else:
        G.add_nodes_from([0])

    pos = hierarchy_pos(G, root=0)
    nx.draw(G, pos=pos, with_labels=True, node_color='orange', alpha=0.7)
    plt.savefig(output + '/phylogeny.png')


def hierarchy_pos(G, root=None, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5):
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5, pos = None, parent = None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''

        if pos is None:
            pos = {root:(xcenter,vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)
        if len(children)!=0:
            dx = width/len(children)
            nextx = xcenter - width/2 - dx/2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G,child, width = dx, vert_gap = vert_gap,
                                    vert_loc = vert_loc-vert_gap, xcenter=nextx,
                                    pos=pos, parent = root)
        return pos


    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)
