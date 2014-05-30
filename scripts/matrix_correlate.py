#!/usr/bin/env python

import sys
import pandas
import numpy
import dill
import csv

if __name__ == "__main__":
    known_training_pairs = dill.load(open('training_pairs.dat', 'rb'))
    known_testing_pairs = dill.load(open('testing_pairs.dat', 'rb'))
    negative_testing_examples = dill.load(open('negative_testing_examples.dat', 'rb'))
    negative_training_examples = dill.load(open('negative_training_examples.dat', 'rb'))

    all_train_pairs = known_training_pairs | negative_training_examples
    all_test_pairs = known_testing_pairs | negative_testing_examples

    writer1 = csv.writer(open('/Users/dhanyasridhar/Documents/psl-pathways/edu.ucsc.cs/data/dev/train/drug.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)
    writer2 = csv.writer(open('/Users/dhanyasridhar/Documents/psl-pathways/edu.ucsc.cs/data/dev/test/drug.txt', 'wb'), delimiter='\t', quoting=csv.QUOTE_NONE)

    matrix = pandas.read_csv('/Users/dhanyasridhar/Documents/yeast_data/HOP_scores.txt', sep="\t", index_col=0)
    #remove any rows where there are nans
    sub_matrix=matrix[~numpy.isnan(matrix).any(axis=1)]
    #get the correlation coefficients of rows
    co=numpy.corrcoef(sub_matrix)
    nz=numpy.nonzero(co>0.207721581179)

    for x, y in zip(nz[0], nz[1]):
        if x != y and (sub_matrix.axes[0][x],sub_matrix.axes[0][y]) in all_train_pairs:
            writer1.writerow([sub_matrix.axes[0][x],sub_matrix.axes[0][y],co[x][y]])
        elif x != y and (sub_matrix.axes[0][x],sub_matrix.axes[0][y]) in all_test_pairs:
            writer2.writerow([sub_matrix.axes[0][x],sub_matrix.axes[0][y],co[x][y]])


