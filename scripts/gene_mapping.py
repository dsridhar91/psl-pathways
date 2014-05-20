__author__ = 'dhanyasridhar'

import csv
import sys
from collections import defaultdict

dir_name = "/Users/dhanyasridhar/Documents/pathways-data/"
file_name = "yeast.txt"

gene_name_map = defaultdict(list)

if __name__ == "__main__":
    map_file = open(dir_name + file_name, 'rb')
    map_reader = csv.reader(map_file, delimiter='\t', quoting=csv.QUOTE_NONE)

    for row in map_reader:
        row = row[0].split(" ")
        print row