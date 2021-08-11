"""Subclonal phylogeny inference"""

import networkx as nx
import random
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

def infer(subclones, score, output):
    G = nx.Graph()

    if len(subclones) > 1:
        edge_ls = []
        root = len(subclones) - 1
        
        snv_assignments = score.values.argmax(axis=1)
        subclone_snv_number = Counter(snv_assignments)
        
        child_parent = {root:root}
        tmp_root = root
        for i in range(len(subclones) - 1):
            i_inverse = root - i - 1
            while True:
                if subclone_snv_number[i_inverse] < subclone_snv_number[tmp_root]:
                    edge_ls.append([tmp_root, i_inverse])
                    child_parent[i_inverse] = tmp_root
                    tmp_root = i_inverse
                    break
                else:
                    if tmp_root == root:
                        edge_ls.append([tmp_root, i_inverse])
                        child_parent[i_inverse] = tmp_root
                        break
                    tmp_root = child_parent[tmp_root]
                 
                    

        G.add_edges_from(edge_ls)
        
        child_in_parent_prop = list()
        for edge in edge_ls:
            p = subclones.iloc[edge[1], 0]/subclones.iloc[edge[0], 0]
            edge.append(p)
            child_in_parent_prop.append(edge)

        pd.DataFrame(child_in_parent_prop, columns=['Parent', 'Child', 'The proportion of child clone']).to_csv(output + '/phylogeny_proportion.csv')

    else:
        G.add_nodes_from([0])

    pos = hierarchy_pos(G, root=len(subclones) - 1)
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
