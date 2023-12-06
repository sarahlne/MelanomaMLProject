""" Module for extracting network/pathway data using rpy2 package """



# import required package
import pandas as pd
import numpy as np
import json

import rpy2.robjects.packages as rpackages
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.conversion import localconverter
import rpy2.robjects as robj


# path to datafiles directory
DATAFILES = '../../source_code/datafiles/'


# Curated Melanoma Network data from the Virtual Cell

def extract_melanoma_vcells(path=DATAFILES+'MelanomaMap_vcells.cyjs'):
    ''' Extract Melanoma network from Virtual cell. '''
    # Load data
    network_data = json.load(open(path))
    # select nodes which are either protein or gene
    nodes = pd.DataFrame(pd.DataFrame(network_data['elements']['nodes'])['data'].to_list())
    nodes = nodes[nodes['class'].isin(['GENE', 'PROTEIN'])][['id', 'hgnc_symbol']]
    nodes.dropna(inplace=True)
    # get edge data for protein and genes
    edges = pd.DataFrame(pd.DataFrame(network_data['elements']['edges'])['data'].to_list())
    edges = edges[['source', 'target']]
    edges = edges[(edges.target.isin(nodes.id.unique())) & (edges.source.isin(nodes.id.unique()))]
    # map hgnc symbols to edge data
    network = pd.DataFrame(columns=['node1', 'node2'])
    network['node1'] = edges['source'].apply(lambda x: nodes[nodes.id == x]['hgnc_symbol'].values[0])
    network['node2'] = edges['target'].apply(lambda x: nodes[nodes.id == x]['hgnc_symbol'].values[0])
    network = network[network.node1 != network.node2].drop_duplicates()
    # processed melanoma network
    # network_rev = pd.DataFrame({'node1' : network['node2'], 'node2': network['node1']})
    # network_processed = pd.concat([network, network_rev], ignore_index=True).drop_duplicates()
    network.reset_index(drop=True, inplace=True)
    
    return network



### Pathway data from KEGG using graphite package

# import graphite package
graphite = rpackages.importr('graphite', lib_loc='C:/Users/kriti/AppData/Local/R/win-library/4.2')

# KEGG pathways
kegg_pathways = graphite.pathways('hsapiens', 'kegg')


# Melanoma Network data from the KEGG database

def extract_melanoma_kegg():
    """ Extract melanoma network from KEGG database. """

    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    # find index of 'Melanoma' pathway
    entry = [i for i, j in enumerate(kegg_pathways.names) if 'Melanoma' in j][0]
    # get melanoma pathway object
    melanoma_pathway = kegg_pathways.slots['entries'][entry]
    # convert all node identifiers from entrez ids to  HGNC symbols
    melanoma_pathway = graphite.convertIdentifiers(melanoma_pathway, 'SYMBOL')
    # Conver R dataframe to pandas dataframe to extract pathway info
    edgedata_list = []
    with localconverter(robj.default_converter + pandas2ri.converter):
        for edge_type in ['protEdges', 'protPropEdges']:
            edge_df = robj.conversion.rpy2py(melanoma_pathway.slots[edge_type])
            edgedata_list.append(edge_df)
    melanoma = pd.concat(edgedata_list, ignore_index=True)

    return melanoma


# Network data from all pathways in KEGG

def extract_pathways_kegg():
    """ Extract gene/proteins in networks from all pathways of KEGG. """

    all_pathways = []
    for entry in kegg_pathways.slots['entries']:
        name = entry.slots['title'][0]
        # print(name)
        pathway = graphite.convertIdentifiers(entry, 'SYMBOL')
        # extract unique nodes in the pathway network
        unique_nodes = set()
        for edge_type in ['protPropEdges', 'protEdges']:
            edge_df = pathway.slots[edge_type]
            src = edge_df[1]
            dest = edge_df[3]
            unique_nodes.update(src, dest)
        # store pathway info
        pathway_nodes = pd.Series([name, list(unique_nodes)], index=['pathway', 'nodes'])
        all_pathways.append(pathway_nodes)
    # dataframe of all pathways 
    all_pathways = pd.DataFrame(all_pathways)
    all_pathways = all_pathways.explode(column='nodes', ignore_index=True).dropna()

    return all_pathways



## Wrapper function to extract all data at once

def extract_data():
    """ Extract network data from resources. """

    # extract processed Melanoma network from Virtual cell
    melanoma_vcells = extract_melanoma_vcells()
    # save extracted data
    melanoma_vcells.to_csv(DATAFILES+'Melanoma_vcells.csv', index=False)

    # extract Melanoma network from KEGG database
    melanoma_kegg = extract_melanoma_kegg()
    # save extracted data
    melanoma_kegg.to_csv(DATAFILES+'Melanoma_kegg.csv', index=False)

    # extract gene/proteins data from all pathways of KEGG.
    pathways = extract_pathways_kegg()
    # save extracted data
    pathways.to_csv(DATAFILES+'Pathways_kegg.csv', index=False)



############ EXTRA FUNCTIONS ############

# functions to convert r objects to python objects and vice versa

from collections import OrderedDict
from rpy2.robjects.vectors import DataFrame, FloatVector, IntVector, StrVector, ListVector, Matrix, FloatArray, FloatMatrix, BoolVector

# Ref :
# https://stackoverflow.com/questions/24152160/converting-an-rpy2-listvector-to-a-python-dictionary
# https://stackoverflow.com/questions/55636932/convert-null-from-python-to-r-using-rpy2
def recurse_r_tree(data):
    """
    step through an R object recursively and convert the types to python types as appropriate. 
    Leaves will be converted to e.g. numpy arrays or lists as appropriate and the whole tree to a dictionary.
    """
    r_dict_types = [DataFrame, ListVector]
    r_array_types = [FloatVector, IntVector, Matrix, FloatArray, FloatMatrix, BoolVector]
    r_list_types = [StrVector]
    if type(data) == robj.rinterface.sexp.NULLType:
        return np.array(None)
    elif type(data) in r_list_types:
        return [recurse_r_tree(elt) for elt in data]
    elif type(data) in r_array_types:
        return np.array(data)
    elif type(data) in r_dict_types:
        if type(data.names) == robj.rinterface.sexp.NULLType:
            pos = list(range(len(data)))
            return OrderedDict(zip(pos, [recurse_r_tree(elt) for elt in data]))
        else:
            return OrderedDict(zip(data.names, [recurse_r_tree(elt) for elt in data]))
    else:
        if hasattr(data, "rclass"):  # An unsupported r class
            raise KeyError('Could not proceed, type {} is not defined'
                           'to add support for this type, just add it to the imports '
                           'and to the appropriate type list above'.format(type(data)))
        else:
            return data  # We reached the end of recursion


# Ref :
# https://stackoverflow.com/questions/14463600/converting-numpy-array-to-rpy2-matrix-forecast-package-xreg-parameter
def nparray2rmatrix(x):
  nr, nc = x.shape
  xvec = robj.FloatVector(x.transpose().reshape((x.size)))
  xr = robj.r.matrix(xvec, nrow=nr, ncol=nc, )
  return xr

