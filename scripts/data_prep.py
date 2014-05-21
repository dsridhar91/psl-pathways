__author__ = 'dhanyasridhar'

import pandas
import numpy
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
import scipy
from numpy import argmax
import itertools
import csv
import sys
import numexpr

score_types = ['protein1', 'protein2', 'neighborhood', 'fusion',  'cooccurence',  'coexpression', 'experimental', 'database', 'textmining']
directory_name = '/Users/dhanyasridhar/Documents/pathways-data/'
output_directory_name = '/Users/dhanyasridhar/Documents/psl-pathways/edu.ucsc.cs/data/'


def main():
    # create_string_similarity_scores()
    # create_ground_truth()
    create_drug_sensitivity_similarity()

def create_string_similarity_scores():
    string_db = 'string-db.txt'

    string_db_matrix = pandas.read_csv(directory_name + string_db, sep='\s')
    scores = pandas.DataFrame(string_db_matrix, columns=score_types)

    scores['protein1'] = scores['protein1'].apply(lambda x: x[5:])
    scores['protein2'] = scores['protein2'].apply(lambda x: x[5:])

    uniqueness = scores['protein1'] != scores['protein2']
    gene_pairs = scores[uniqueness][['protein1', 'protein2']]
    gene_pairs = gene_pairs.drop_duplicates(['protein1', 'protein2'])
    gene_pairs.to_csv(output_directory_name + 'pathwayPairs.txt', sep='\t', header=False, index=False)

    for score in (set(score_types) - set(['protein1', 'protein2'])):
        column = scores[score]
        max_value = column.max(axis=1)

        scores[score] = column.apply(lambda x: float(x)/float(max_value))

        non_zero = scores[score] != 0
        similarityCol = scores[non_zero][['protein1', 'protein2', score]]
        similarityCol.to_csv(output_directory_name + score + '.txt', sep='\t', header=False, index=False)


def create_ground_truth():
    pathways_file = 'yeast_pathways.txt'
    pathways = csv.reader(open(directory_name + pathways_file, 'rb'), delimiter='\t')

    known_pathway_pairs = csv.writer(open(output_directory_name + 'knownPathwayPairs.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    for pathway in pathways:
        for i in xrange(1, len(pathway)):
            for j in xrange(i + 1, len(pathway)):
                known_pathway_pairs.writerow([pathway[i], pathway[j]])


def create_drug_sensitivity_similarity():
    matrix = pandas.read_csv('/Users/dhanyasridhar/Documents/yeast_data/HOP_scores.txt', sep="\t", index_col=0)
    #remove any rows where there are nans
    sub_matrix=matrix[~numpy.isnan(matrix).any(axis=1)]
    #get the correlation coefficients of rows
    co=numpy.corrcoef(sub_matrix)
    nz=numpy.nonzero(co>0)

    drug_similarity = csv.writer(open(output_directory_name + 'drugSensitivity.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    for x, y in zip(nz[0], nz[1]):
        if x != y:
            drug_similarity.writerow([sub_matrix.axes[0][x], sub_matrix.axes[0][y], co[x][y]])
            # print sub_matrix.axes[0][x], sub_matrix.axes[0][y], co[x][y]


main()