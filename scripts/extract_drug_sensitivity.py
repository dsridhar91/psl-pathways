__author__ = 'dhanyasridhar'

import csv
import sys
import string
import numpy
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

yeast_file_name = "/Users/dhanyasridhar/Documents/yeast_data/HOP_scores.txt"
sensitivity_file_name = "/Users/dhanyasridhar/Documents/yeast_data/sensitivity_scores.txt"
prefix = ['Ad.', 'MADL']
suffix = 'z-score'
drug_index_dict = defaultdict(lambda: defaultdict())
gene_index_dict = defaultdict()
sensitivity_matrix = numpy.array


def main():
    initialize_data()


def create_sensitivity_matrix():

    yeast_file = open(yeast_file_name, 'rb')
    HOP_reader = csv.reader(yeast_file, delimiter='\t', quoting=csv.QUOTE_NONE)

    matrix_file = open(sensitivity_file_name, 'wb')
    scores_writer = csv.writer(matrix_file, delimiter='\t', quoting=csv.QUOTE_NONE)

    header_row = list()

    for index, dictionary in drug_index_dict.items():
        header_row.append(dictionary['cmd_id'] + '_' + dictionary['exp'] + '_' + dictionary['score_type'])

    scores_writer.writerow(header_row)

    for row in HOP_reader:
        gene_name = row[0]
        gene_name = gene_name[1:-1]

        matrix_row = list()
        matrix_row.append(gene_name)

        for index, dictionary in drug_index_dict.items():
            matrix_row.append(row[index])

        scores_writer.writerow(matrix_row)

    yeast_file.close()
    matrix_file.close()


def initialize_data():
    yeast_file = open(yeast_file_name, 'rb')
    HOP_reader = csv.reader(yeast_file, delimiter='\t', quoting=csv.QUOTE_NONE)

    header = HOP_reader.next()
    for i in xrange(len(header)):
        colName = header[i][1:-1]
        colName = colName.split(" ")
        type_of_score = colName[0]

        if type_of_score in prefix and colName[-1] != suffix:
            identifier = colName[4]
            identifier = identifier.split("_")
            drug_compound_id = identifier[0]
            drug_concentration = identifier[1]
            drug_experiment_no = identifier[3]

            drug_dict = defaultdict()
            drug_dict['cmd_id'] = drug_compound_id
            drug_dict['conc'] = drug_concentration
            drug_dict['exp'] = drug_experiment_no
            drug_dict['score_type'] = type_of_score

            drug_index_dict[i] = drug_dict

    for i, row in enumerate(HOP_reader):
        gene_name = row[0]
        gene_name = gene_name[1:-1]

        gene_index_dict[i] = gene_name

        for j, dictionary in drug_index_dict.items():
            sensitivity_matrix[i][j] = row[j]

    yeast_file.close()

def output_sensitivity_predicate_data():
    sensitivity_file = open(sensitivity_file_name, 'rb')


main()