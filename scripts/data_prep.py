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
from collections import defaultdict
import random
import dill

score_types = ['protein1', 'protein2', 'neighborhood', 'fusion',  'cooccurence',  'coexpression', 'experimental', 'database', 'textmining']
directory_name = '/Users/dhanyasridhar/Documents/pathways-data/'
output_directory_name = '/Users/dhanyasridhar/Documents/psl-pathways/edu.ucsc.cs/data/'
train = 'train/'
test = 'test/'
training_pairs = set()
testing_pairs = set()
pathways = defaultdict(lambda: set())


def main():
#    collect_pathways()
#    create_splits()
#    create_string_similarity_scores()
     initSplits()
     create_drug_sensitivity_similarity()
#    create_ground_truth()


def initSplits():
    global training_pairs, testing_pairs
    training_pairs = dill.load(open('training_pairs.dat', 'rb'))
    testing_pairs = dill.load(open('testing_pairs.dat', 'rb'))


def create_string_similarity_scores():
    string_db = 'string-db.txt'

    string_db_matrix = pandas.read_csv(directory_name + string_db, sep='\s')
    scores = pandas.DataFrame(string_db_matrix, columns=score_types)

    scores['protein1'] = scores['protein1'].apply(lambda x: x[5:])
    scores['protein2'] = scores['protein2'].apply(lambda x: x[5:])

    # queries for different conditions we want
    uniqueness = scores['protein1'] != scores['protein2']

    # random variables for the test split
    gene_pairs = scores[uniqueness][['protein1', 'protein2']]
    gene_dict = gene_pairs.to_dict()
    genes = csv.writer(open(output_directory_name + test + 'pathwayPairs.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    for index, value in gene_dict['protein1'].items():
            gene_pair = (value, gene_dict['protein2'][index])
            if gene_pair not in training_pairs:
                genes.writerow([gene_pair[0], gene_pair[1]])

    for score in (set(score_types) - set(['protein1', 'protein2'])):
        column = scores[score]
        max_value = column.max(axis=1)

        scores[score] = column.apply(lambda x: float(x)/float(max_value))

        non_zero = scores[score] != 0

        # add to training split
        similarityCol = scores[non_zero][['protein1', 'protein2', score]]
        similarityDict = similarityCol.to_dict()
        similarityTuples = set()

        for index, value in similarityDict['protein1'].items():
            gene_pair = (value, similarityDict['protein2'][index])
            if gene_pair in training_pairs:
                score_tuple = (gene_pair[0], gene_pair[1], similarityDict[score][index])
                similarityTuples.add(score_tuple)


        training = csv.writer(open(output_directory_name + train + score + '.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)
        for tuple in similarityTuples:
            training.writerow([tuple[0], tuple[1], tuple[2]])

        # add to testing split
        similarityCol = scores[non_zero][['protein1', 'protein2', score]]
        similarityDict = similarityCol.to_dict()
        similarityTuples = set()

        for index, value in similarityDict['protein1'].items():
            gene_pair = (value, similarityDict['protein2'][index])
            if gene_pair not in training_pairs:
                score_tuple = (gene_pair[0], gene_pair[1], similarityDict[score][index])
                similarityTuples.add(score_tuple)


        testing = csv.writer(open(output_directory_name + test + score + '.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)
        for tuple in similarityTuples:
            testing.writerow([tuple[0], tuple[1], tuple[2]])


def create_ground_truth():
    train_truth = csv.writer(open(output_directory_name + train + 'knownPathwayPairs.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)
    test_truth = csv.writer(open(output_directory_name + test + 'knownPathwayPairs.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    for gene_pair in training_pairs:
        train_truth.writerow([gene_pair[0], gene_pair[1]])

    for gene_pair in testing_pairs:
        test_truth.writerow([gene_pair[0], gene_pair[1]])


def create_drug_sensitivity_similarity():
    matrix = pandas.read_csv('/Users/dhanyasridhar/Documents/yeast_data/HOP_scores.txt', sep="\t", index_col=0)
    #remove any rows where there are nans
    sub_matrix=matrix[~numpy.isnan(matrix).any(axis=1)]
    score_matrix=sub_matrix.transpose()[ [not a.endswith('z-score') for a in matrix.axes[1]] ].transpose()
    #get the correlation coefficients of rows
    co=numpy.corrcoef(score_matrix)
    nz=numpy.nonzero(co>0)

    # training split
    training_drug_similarity = csv.writer(open(output_directory_name + train + 'drugSensitivity.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)
    testing_drug_similarity = csv.writer(open(output_directory_name + test + 'drugSensitivity.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    for x, y in zip(nz[0], nz[1]):
        if x != y and (x, y) in training_pairs:
            training_drug_similarity.writerow([sub_matrix.axes[0][x], sub_matrix.axes[0][y], co[x][y]])
            # print sub_matrix.axes[0][x], sub_matrix.axes[0][y], co[x][y]
        elif x != y and not (x, y) in training_pairs:
            testing_drug_similarity.writerow([sub_matrix.axes[0][x], sub_matrix.axes[0][y], co[x][y]])


def collect_pathways():
    pathways_file = 'yeast_pathways.txt'
    pathways_reader = csv.reader(open(directory_name + pathways_file, 'rb'), delimiter='\t')

    for pathway in pathways_reader:
        for i in xrange(1, len(pathway)):
            for j in xrange(i + 1, len(pathway)):
                if pathway[i] != pathway[j]:
                    pathways[pathway[0]].add((pathway[i], pathway[j]))


def create_splits():
    is_training = True

    while len(pathways.keys()) > 0:
        if is_training:
            p_name = random.choice(pathways.keys())
            random_pathway = pathways.pop(p_name)
            for gene_pair in random_pathway:
                training_pairs.add(gene_pair)
                is_training = False
        else:
            p_name = random.choice(pathways.keys())
            random_pathway = pathways.pop(p_name)
            for gene_pair in random_pathway:
                testing_pairs.add(gene_pair)
                is_training = True

    overlap = training_pairs & testing_pairs
    for item in overlap:
        testing_pairs.remove(item)

    dill.dump(training_pairs, open('training_pairs.dat', 'wb'))
    dill.dump(testing_pairs, open('testing_pairs.dat', 'wb'))


main()