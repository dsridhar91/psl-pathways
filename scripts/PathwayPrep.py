# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import urllib
import csv

# <codecell>

urllib.urlretrieve("http://downloads.yeastgenome.org/curation/chromosomal_feature/dbxref.tab", "dbxref.tab")

# <codecell>

symbol_map = {}
with open("dbxref.tab") as handle:
    reader = csv.reader(handle, delimiter="\t")
    for line in reader:
        symbol_map[line[5]] = line[3]

# <codecell>

urllib.urlretrieve("http://downloads.yeastgenome.org/curation/literature/biochemical_pathways.tab", "biochemical_pathways.tab")

# <codecell>

path_map = {}
with open("biochemical_pathways.tab") as handle:
    reader = csv.reader(handle, delimiter="\t")
    for line in reader:
        if line[3] in symbol_map:
            gene = symbol_map[line[3]]
            path_map[line[0]] = path_map.get(line[0], []) + [gene]

# <codecell>

with open("pathway.gmt", "w") as handle:
    for a in path_map:
        handle.write("%s\t%s\n" % (a, "\t".join(path_map[a])))

# <codecell>


