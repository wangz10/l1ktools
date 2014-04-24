'''
This script contains examples for interacting with the LINCS API programatically
through Python
'''

import cmap.util.api_utils as apiu

# set up a class to handle API calls
ac = apiu.APIContainer()

# show collections available to the container class
print ac

# search for the first 20 signatures treated with Vemurafenib
# return their sig_id, pert_id, pert_itime, pert_iname, and pert_idose
# return result as Pandas DataFrame
sigs = ac.siginfo.find({'pert_iname' : 'vemurafenib'},
                       fields = ['sig_id', 'pert_id', 'pert_itime', 
                      			 'pert_idose', 'pert_iname'],
                       toDataFrame = True, limit = 20)

# get information on all HDAC genes
genes = ac.geneinfo.find({'pr_gene_symbol' : {'$regex' : 'HDAC.*'},
                          'is_l1000' : True},
                         fields = ['pr_id', 'pr_gene_symbol', 'pr_gene_title'],
                         toDataFrame = True, limit = 20)

# get information on some Statin drugs
# return all available information as a list of dictionaries
drugnames = ['atorvastatin', 'simvastatin', 'fluvastatin']
drugs = ac.pertinfo.find({'pert_iname' : {'$in' : drugnames}})

# or, return a few select fields as a DataFrame
drugs_alt = ac.pertinfo.find({'pert_iname' : {'$in' : drugnames}},
                             fields = ['num_sig', 'pert_icollection',
                             		   'pert_id', 'pert_summary'],
                             toDataFrame = True)