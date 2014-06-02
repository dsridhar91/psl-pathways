__author__ = 'dhanyasridhar'
import sys
import pandas
import numpy
import dill
import csv

score_types = ['protein1', 'protein2', 'neighborhood', 'fusion',  'cooccurence',  'coexpression', 'experimental', 'database', 'textmining']
directory_name = '/Users/dhanyasridhar/Documents/pathways-data/'
output_directory_name = '/Users/dhanyasridhar/Documents/psl-pathways/edu.ucsc.cs/data/'
train = 'train/'
test = 'test/'


def main():

    matrix = pandas.read_csv('/Users/dhanyasridhar/Documents/yeast_data/HOP_scores.txt', sep="\t", index_col=0)
    #remove any rows where there are nans
    sub_matrix=matrix[~numpy.isnan(matrix).any(axis=1)]
    #get the correlation coefficients of rows
    co=numpy.corrcoef(sub_matrix)

    co_sub = numpy.tril(co, -1)
    co_sub_1d = co_sub.flatten()
    co_sub_1d[co_sub_1d <= 0] = numpy.NaN
    noNan = co_sub_1d[~numpy.isnan(co_sub_1d)]

    sorted_co = numpy.sort(noNan)
    threshold_index = 0.96 * len(sorted_co)
    print sorted_co[threshold_index]

    threshold_dict= {}
    print "Reading from string db file.."
    string_db = 'string-db.txt'

    string_db_matrix = pandas.read_csv(directory_name + string_db, sep='\s')
    scores = pandas.DataFrame(string_db_matrix, columns=score_types)

    for score in (set(score_types) - set(['protein1', 'protein2'])):
        print "reading values for " + score + " column.."
        column = scores[score]
        max_value = column.max(axis=1)
        scores[score] = column.apply(lambda x: float(x)/float(max_value))
        non_zero = scores[score] != 0

        similarityCol = scores[non_zero][['protein1', 'protein2', score]]

        sorted_col = numpy.sort(similarityCol[score])

        threshold_index = 0.96 * len(sorted_col)
        threshold_value = sorted_col[threshold_index]

        threshold_dict[score] = threshold_value

        dill.dump(threshold_dict, open('threshold_dict.dat', 'wb'))

main()