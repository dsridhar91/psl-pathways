#!/usr/bin/env python

import sys
import pandas
import numpy

if __name__ == "__main__":
	matrix = pandas.read_csv(sys.argv[1], sep="\t", index_col=0)
	#remove any rows where there are nans
	sub_matrix=matrix[~numpy.isnan(matrix).any(axis=1)]
	#get the correlation coefficients of rows
	co=numpy.corrcoef(sub_matrix)
	nz=numpy.nonzero(co>0.7)

	for x, y in zip(nz[0], nz[1]):
		if x != y:
			print "%s\t%s\%s" % (sub_matrix.axes[0][x],sub_matrix.axes[0][y],co[x][y])

