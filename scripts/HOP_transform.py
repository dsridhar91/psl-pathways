#!/usr/bin/env python


import pandas
import numpy
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
import itertools

hop_scores=pandas.read_csv("HOP_scores.txt", sep="\t", index_col=0)
sub_matrix=hop_scores[~numpy.isnan(hop_scores).any(axis=1)]
score_matrix=sub_matrix.transpose()[ [not a.endswith('z-score') for a in hop_scores.axes[1]] ].transpose()
score_matrix.to_csv("HOP_scores.txt.fixed", sep="\t")

distxy = squareform(pdist(score_matrix, metric='correlation'))

link = linkage(distxy, method='complete')


def load_gmt(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            tmp = line.rstrip().split("\t")
            out[tmp[0]] = tmp[1:]
    return out


out = load_gmt("yeast_pathways.gmt")

def intersect(a,b):
    out = {}
    for i in a:
       out[i] = True
    for i in b:
        if i in out:
            yield i




# <codecell>

def cluster_to_matrix(labels, clust):
    df = pandas.DataFrame(index=labels, columns=labels).fillna(0)
    for i in numpy.unique(clust):
        for x,y in itertools.combinations(labels[clust == i], 2):
            df[x][y] = 1
            df[y][x] = 1
    for i in labels:
        df[i][i] = numpy.nan
    return df

# <codecell>

def gmt_to_matrix(gmt):
    key_map = {}
    for a in gmt.values():
        for i in a:
            key_map[i] = True
    keys = key_map.keys()
    df = pandas.DataFrame(index=keys, columns=keys).fillna(0)
    for a in gmt.values():
        for x,y in itertools.combinations(a, 2):
            df[x][y] = 1
            df[y][x] = 1
    for a in keys:
        df[a][a] = numpy.nan
    return df


# <codecell>

def accuracy(data, pred):
    #true positive: 1 in both matrices
    tp = float(sum(((data+pred) == 2).values))
    #false negative: 0 in prediction 1 in data
    fn = float(sum(((pred-data) == -1).values) )
    #true negative: 0 in both matrices
    tn = float(sum((data+pred == 0).values))
    #false positives: 1 in prediction 0 in data
    fp = float(sum( (pred-data == 1).values))
    out = {
        'tp' : tp,
        'fn' : fn,
        'tn' : tn,
        'fp' : fp
    }
    #print out
    return out


# <codecell>

out_df = gmt_to_matrix(out)

# <codecell>


# <codecell>

for i in numpy.arange(0.5, 1.2, 0.05):
    l = fcluster(link,i)
    df = cluster_to_matrix(score_matrix.axes[0], l)
    a = accuracy(out_df, df)
    print i, a
    sens = a['tp'] / (a['tp']+a['fn'])
    spes = a['tn'] / (a['fp']+a['tn'])
    fdr =  a['fp'] / (a['tp']+a['fp'])

    print i, sens, spes, fdr
