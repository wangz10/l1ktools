'''
This script contains examples for reading .gctx files in Python.
'''

import cmap.io.gct as gct

# give input file
path_to_gctx_file = '/xchip/cogs/l1ktools/data/modzs_n272x978.gctx'

# read the full data file
GCTObject = gct.GCT(path_to_gctx_file)
GCTObject.read()
print(GCTObject.matrix)

# read the first 100 rows and 10 columns of the data
GCTObject = gct.GCT(path_to_gctx_file)
GCTObject.read(row_inds=range(100),col_inds=range(10))
print(GCTObject.matrix)

# get the available meta data headers for data columns and row
column_headers = GCTObject.get_chd()
row_headers = GCTObject.get_rhd()

# get the perturbagen description meta data field from the column data
inames = GCTObject.get_column_meta('pert_iname')

# get the gene symbol meta data field from the row data
symbols = GCTObject.get_row_meta('pr_gene_symbol')

GCTObject.write('/xchip/cogs/l1ktools/data/python_example.gctx')
