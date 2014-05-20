__author__ = 'dhanyasridhar'

import pandas
import numpy
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
import scipy
from numpy import argmax
import itertools
import csv
import numexpr

score_types = ['protein1', 'protein2', 'neighborhood', 'fusion',  'cooccurence',  'coexpression', 'experimental', 'database', 'textmining']


def main():
    directory_name = '/Users/dhanyasridhar/Documents/pathways-data/'
    string_db = 'string-db.txt'

    string_db_matrix = pandas.read_csv(directory_name + string_db, sep='\s')
    scores = pandas.DataFrame(string_db_matrix, columns=score_types)

    scores['protein1'] = scores['protein1'].apply(lambda x: x[5:])
    scores['protein2'] = scores['protein2'].apply(lambda x: x[5:])

    uniqueness = scores['protein1'] != scores['protein2']
    gene_pairs = scores[uniqueness][['protein1', 'protein2']]
    gene_pairs.to_csv(directory_name + 'pathwayPairs.txt', sep='\t')

    for score in (set(score_types) - set(['protein1', 'protein2'])):
        column = scores[score]
        max_value = column.max(axis=1)

        scores[score] = column.apply(lambda x: float(x)/float(max_value))

        non_zero = scores[score] != 0
        similarityCol = scores[non_zero][['protein1', 'protein2', score]]
        # print similarityCol
        similarityCol.to_csv(directory_name + score + '.txt', sep='\t')


main()