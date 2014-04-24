'''
This script contains examples for reading .gctx files in Python, and making
calls to the LINCS API annotations via Python.
'''

import cmap.io.gct as gct
import cmap.util.api_utils as apiu

# give input file
path_to_gctx_file = '../data/modzs_n272x978.gctx'

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

# get annotations by calling CMAP web API
# this requires a LINCS API key. If you do not have one, contact:
# lincs@broadinstitute.org
ac = apiu.APIContainer()
row_annots = ac.geneinfo.find({'pr_id' : {'$in' : GCTObject.get_rids()}},
                              limit = 100, toDataFrame = True,
                              fields = ['pr_id', 'pr_gene_symbol',
                              			'pr_gene_title', 'is_lm'])
print row_annots.head()

col_annots = ac.siginfo.find({'sig_id' : {'$in' : GCTObject.get_cids()}},
                             limit = 10, toDataFrame = True,
                             fields = ['sig_id', 'pert_id', 'pert_iname',
                             		   'pert_itime', 'pert_idose'])
print col_annots
