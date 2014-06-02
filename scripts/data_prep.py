__author__ = 'dhanyasridhar'

#TODO: use command line args for input, output directory names and train/test split percentage

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
import shutil
import os

score_types = ['protein1', 'protein2', 'neighborhood', 'fusion',  'cooccurence',  'coexpression', 'experimental', 'database', 'textmining']
directory_name = '/Users/dhanyasridhar/Documents/pathways-data/'
output_directory_name = '/Users/dhanyasridhar/Documents/psl-pathways/edu.ucsc.cs/data/'
train = 'train/'
test = 'test/'
dev = 'dev/'

known_training_pairs = set()
known_testing_pairs = set()
unknown_pairs = set()

negative_training_examples = set()
negative_testing_examples = set()

pathways = defaultdict(lambda: set())
similarity_scores = defaultdict(lambda: defaultdict())

drug_corr_threshold = 0.207721581179
threshold_dict ={}

def main():
    initSplits()
    init_negative_examples()

    print "Reading from string db file.."
    rv_pairs = set()
    string_db = 'string-db.txt'

    string_db_matrix = pandas.read_csv(directory_name + string_db, sep='\s')
    scores = pandas.DataFrame(string_db_matrix, columns=score_types)

    # strips numeric prefix from gene names
    scores['protein1'] = scores['protein1'].apply(lambda x: x[5:])
    scores['protein2'] = scores['protein2'].apply(lambda x: x[5:])

    # queries to make sure gene1 is not gene2
    uniqueness = scores['protein1'] != scores['protein2']

    # random variables (gene pairs) for the test split
    gene_pairs = scores[uniqueness][['protein1', 'protein2']]
    gene_dict = gene_pairs.to_dict()

    print "setting up random variable pairs.."

    all_known_pairs = (known_training_pairs | negative_training_examples | known_testing_pairs | negative_testing_examples)

    for index, value in gene_dict['protein1'].items():
            gene_pair = (value, gene_dict['protein2'][index])
            inverse_pair = (gene_dict['protein2'][index], value)
            if gene_pair not in all_known_pairs:
                if inverse_pair not in rv_pairs:
                    rv_pairs.add(gene_pair)

    one_hop_pairs = get_one_hop_unknown_pathways(rv_pairs, known_testing_pairs)

    print "one hop from only known pathway genes region, ", len(one_hop_pairs)

    print 'length of all unknown pathways', len(rv_pairs)

    print 'length of known testing pathways, ', len(known_testing_pairs)

    print 'no of negative example pairs in test, ', len(negative_testing_examples)

    # setup_data()
    # create_dev_output()
    # create_dev_ground_truth()

    # co_sub = numpy.tril(co, -1)
    # co_sub_1d = co_sub.flatten()
    # co_sub_1d[co_sub_1d <= 0] = numpy.NaN
    # noNan = co_sub_1d[~numpy.isnan(co_sub_1d)]
    #
    # sorted_co = numpy.sort(noNan)
    # print sorted_co[4398988]


    # threshold_dict= {}
    # print "Reading from string db file.."
    # string_db = 'string-db.txt'
    #
    # string_db_matrix = pandas.read_csv(directory_name + string_db, sep='\s')
    # scores = pandas.DataFrame(string_db_matrix, columns=score_types)
    #
    # # set up similarity data
    # for score in (set(score_types) - set(['protein1', 'protein2'])):
    #     print "reading values for " + score + " column.."
    #     column = scores[score]
    #     max_value = column.max(axis=1)
    #     scores[score] = column.apply(lambda x: float(x)/float(max_value))
    #     non_zero = scores[score] != 0
    #
    #     similarityCol = scores[non_zero][['protein1', 'protein2', score]]
    #
    #     sorted_col = numpy.sort(similarityCol[score])
    #
    #     threshold_index = 0.96 * len(sorted_col)
    #     threshold_value = sorted_col[threshold_index]
    #
    #     threshold_dict[score] = threshold_value
    #
    #     dill.dump(threshold_dict, open('threshold_dict.dat', 'wb'))


def initSplits():
    global known_training_pairs, known_testing_pairs, threshold_dict
    known_training_pairs = dill.load(open('training_pairs.dat', 'rb'))
    known_testing_pairs = dill.load(open('testing_pairs.dat', 'rb'))

    threshold_dict = dill.load(open('threshold_dict.dat', 'rb'))


def init_negative_examples():
    global negative_testing_examples, negative_training_examples

    negative_testing_examples = dill.load(open('negative_testing_examples.dat', 'rb'))
    negative_training_examples = dill.load(open('negative_training_examples.dat', 'rb'))


def setup_data_for_known_pathways_only():
    print "Reading from string db file.."
    string_db = 'string-db.txt'

    string_db_matrix = pandas.read_csv(directory_name + string_db, sep='\s')
    scores = pandas.DataFrame(string_db_matrix, columns=score_types)

    print "setting up similarity dict.."
    # set up similarity data
    for score in (set(score_types) - set(['protein1', 'protein2'])):
        column = scores[score]
        max_value = column.max(axis=1)

        scores[score] = column.apply(lambda x: float(x)/float(max_value))
        non_zero = scores[score] != 0
        similarityCol = scores[non_zero][['protein1', 'protein2', score]]

        above_threshold = similarityCol[score] > float(threshold_dict[score])
        blocked_col = similarityCol[above_threshold][['protein1', 'protein2', score]]

        blocked_col_dict = blocked_col.to_dict()

        for index, value in blocked_col_dict['protein1'].items():
            gene_pair = (value, blocked_col_dict['protein2'][index])
            similarity_scores[score][gene_pair] = blocked_col_dict[score][index]

    print "finished data set up.."


def create_dev_output():
    print "writing to files.."
    genes = csv.writer(open(output_directory_name + dev + train + 'pathwayPairs.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    for pair in known_training_pairs:
        genes.writerow([pair[0], pair[1]])
    for pair in negative_training_examples:
        genes.writerow([pair[0], pair[1]])

    genes = csv.writer(open(output_directory_name + dev + test + 'pathwayPairs.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    for pair in known_testing_pairs:
        genes.writerow([pair[0], pair[1]])
    for pair in negative_testing_examples:
        genes.writerow([pair[0], pair[1]])

    for score, similarity in similarity_scores.items():
        training = csv.writer(open(output_directory_name + dev + train + score + '.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

        for pair in known_training_pairs:
            if pair in similarity:
                training.writerow([pair[0], pair[1], similarity[pair]])
        for pair in negative_training_examples:
            if pair in similarity:
                training.writerow([pair[0], pair[1], similarity[pair]])

        testing = csv.writer(open(output_directory_name + dev + test + score + '.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

        for pair in known_testing_pairs:
            if pair in similarity:
                testing.writerow([pair[0], pair[1], similarity[pair]])
        for pair in negative_testing_examples:
            if pair in similarity:
                testing.writerow([pair[0], pair[1], similarity[pair]])


def create_dev_ground_truth():
    train_truth = csv.writer(open(output_directory_name + dev + train + 'knownPathwayPairs.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)
    test_truth = csv.writer(open(output_directory_name + dev + test + 'knownPathwayPairs.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    for gene_pair in known_training_pairs:
        train_truth.writerow([gene_pair[0], gene_pair[1], 1])
    for gene_pair in negative_training_examples:
        train_truth.writerow([gene_pair[0], gene_pair[1], 0])
        
    for gene_pair in known_testing_pairs:
        test_truth.writerow([gene_pair[0], gene_pair[1], 1])
    for gene_pair in negative_testing_examples:
        test_truth.writerow([gene_pair[0], gene_pair[1], 0])


def collect_pathways():
    pathways_file = 'yeast_pathways.txt'
    pathways_reader = csv.reader(open(directory_name + pathways_file, 'rb'), delimiter='\t')

    for pathway in pathways_reader:
        for i in xrange(1, len(pathway)):
            for j in xrange(i + 1, len(pathway)):
                if pathway[i] != pathway[j]:
                    pathways[pathway[0]].add((pathway[i], pathway[j]))

    dill.dump(pathways, open('pathways.dat', 'wb'))


def create_splits():
    is_training = True

    while len(pathways.keys()) > 0:
        if is_training:
            p_name = random.choice(pathways.keys())
            random_pathway = pathways.pop(p_name)
            for gene_pair in random_pathway:
                known_training_pairs.add(gene_pair)
                is_training = False
        else:
            p_name = random.choice(pathways.keys())
            random_pathway = pathways.pop(p_name)
            for gene_pair in random_pathway:
                known_testing_pairs.add(gene_pair)
                is_training = True

    overlap = known_training_pairs & known_testing_pairs
    for item in overlap:
        known_testing_pairs.remove(item)

    dill.dump(known_training_pairs, open('training_pairs.dat', 'wb'))
    dill.dump(known_testing_pairs, open('testing_pairs.dat', 'wb'))


def collect_negative_examples():
    all_training_genes = set()
    all_testing_genes = set()
    all_known_training_pairs = set()
    all_known_testing_pairs = set()

    for gene_pair in known_training_pairs:
        all_training_genes.add(gene_pair[0])
        all_training_genes.add(gene_pair[1])
    
#    for g in all_training_genes:
#        for g2 in all_training_genes:
#            if g != g2:
#                all_known_training_pairs.add((g, g2))

    for i in xrange(1, len(all_training_genes)):
        for j in xrange(i+1, len(all_training_genes)):
            all_known_training_pairs.add((list(all_training_genes)[i], list(all_training_genes)[j]))

    for gene_pair in known_testing_pairs:
        all_testing_genes.add(gene_pair[0])
        all_testing_genes.add(gene_pair[1])
        
#    for g in all_testing_genes:
#        for g2 in all_testing_genes:
#            if g != g2:
#                all_known_testing_pairs.add((g, g2))

    for i in xrange(1, len(all_testing_genes)):
        for j in xrange(i+1, len(all_testing_genes)):
            all_known_testing_pairs.add((list(all_testing_genes)[i], list(all_testing_genes)[j]))

    global negative_training_examples, negative_testing_examples

    negative_training_examples = all_known_training_pairs - known_training_pairs
    negative_testing_examples = all_known_testing_pairs - known_testing_pairs

    dill.dump(negative_training_examples, open('negative_training_examples.dat', 'wb'))
    dill.dump(negative_testing_examples, open('negative_testing_examples.dat', 'wb'))


main()