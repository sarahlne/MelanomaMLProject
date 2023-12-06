""" Module to integrate the extracted pathways data with clinical data """



# import required packages
import os
import numpy as np
import pandas as pd

from scipy.stats import hypergeom

from load_data import load_snps, load_clinical_raw



############ CHECK REQUIRED FILES ############

# datafiles directory
DATAFILES = '../../source_code/datafiles/'
# files with extracted pathways information
REQUIRED_FILES = [
    'Melanoma_vcells.csv',
    'Melanoma_kegg.csv',
    'Pathways_kegg.csv'
]

# files present in the directory
PRESENT_FILES = os.listdir(DATAFILES)
# check for required files in the directory
if len(np.intersect1d(REQUIRED_FILES, PRESENT_FILES)) < len(REQUIRED_FILES):
    print('Extracting data from resources.')
    from extract_pathway_data import extract_data
    extract_data()
    print(f'Extracted data stored in "{DATAFILES}" directory')



############ MUTATION DATA ############

# function to get list of all mutated genes in patient from clinical and snps data
# used to integrate network/pathways data with patients data

def patients_mut_genes(type = None):
    """ Returns mapping of all mutated genes to respective patients. """
    if (type=='validation'): # Load mut genes for Blateau & Louveau patients
        snps_data = load_snps(type='validation').copy()
        sources_to_keep = ['Pauline Blateau, Jerome Solassol', 'Baptiste Louveau, Samia Mourah'] #
        snps_data = snps_data.loc[snps_data['source'].isin(sources_to_keep)]
        snps_data = snps_data[['patient_id', 'HGNC', 'mutated']]
        snps_data = snps_data[snps_data.mutated == 'yes'].drop(columns='mutated').drop_duplicates()
    # # BRAF mutation from clinical data
    # patients_braf = load_clinical_raw()[['patient_ID', 'BRAF_mut']].copy()
    # patients_braf.dropna(subset=['BRAF_mut'], inplace=True) # mutation info in snps data
    # patients_braf['BRAF_mut'] = patients_braf['BRAF_mut'].map({'MND' : 0}).fillna(1)
    # patients_braf = patients_braf[patients_braf.BRAF_mut != 0] # patients with no mutation info removed
    # patients_braf.loc[:, 'HGNC'] = 'BRAF'
    # patients_braf.drop(columns=['BRAF_mut'], inplace=True)
    # patients_braf.rename(columns={'patient_ID' : 'patient_id'}, inplace=True)
    # all mutated genes in patients from snps data
    else:
        snps_data = load_snps().copy()
        snps_data = snps_data[snps_data.source != 'Pauline Blateau, Jerome Solassol']
        snps_data = snps_data[['patient_id', 'HGNC', 'mutated']]
        snps_data = snps_data[snps_data.mutated == 'yes'].drop(columns='mutated').drop_duplicates()

    # dataframe with list of all mutated genes of patients
    # mutation_data = pd.concat([patients_braf, snps_data], ignore_index=True)
    # mutation_data.drop_duplicates(inplace=True)
    #mutation_data = mutation_data.groupby(['patient_id'], as_index=False).agg({'HGNC': lambda x: x.tolist()})

    # # For Cscape analysis - consider only cscape genes for analysis
    # cscape_data = pd.read_csv(DATAFILES+'Binary_cscape_oncogenic.csv', index_col=0)
    # cscape_genes_df = pd.DataFrame()
    # cscape_genes_df['HGNC'] = cscape_data[cscape_data.columns].astype('bool').apply(lambda row: list(cscape_data.columns[row]), axis=1)
    # cscape_genes_df = cscape_genes_df.explode('HGNC')
    # cscape_genes_df.reset_index(inplace=True)
    # cscape_genes_df.rename(columns={'index':'patient_id'}, inplace=True)

    return snps_data # cscape_genes_df

############ MELANOMA NETWORK ############

# create binary table of mutated genes respective to patients
def mut_genes_table(mutation_data, network_genes):
    """
    Function to create binary table of mutated genes in patient.

    arguments:
    df_snps : dataframe with patient id, HGNC symbols and mutation status
    network_genes : genes common between the patients and preferred pathway

    returns:
    pd.DataFrame : binary table of mutated genes

    """
    # get list of mutated genes for each patient
    mutated_df = mutation_data[mutation_data['HGNC'].isin(network_genes)].copy(deep=True)
    mutated_df.drop_duplicates(inplace=True)
    mutated_df = mutated_df.groupby(['patient_id'], as_index=False).agg({'HGNC': lambda x: x.tolist()})
    # create binary table of mutated genes for patients
    mutation_dict = dict()
    for _, (id, genes) in mutated_df.iterrows():
        if id not in mutation_dict:
            mutation_dict[id] = pd.Series([1 for i in genes], index=genes)
    binary_table = pd.DataFrame(mutation_dict).T#.fillna(0)
    binary_table = pd.DataFrame(binary_table, columns=network_genes)
    binary_table = binary_table[sorted(binary_table.columns)]
    return binary_table


# map patients snp data to melanoma network from kegg
def map_melanoma_kegg(filename = 'Binary_kegg.csv', type = 'None'):
    if type=='validation':
        filename = 'Binary_kegg_validation.csv'
        mutation_data = patients_mut_genes(type='validation').copy()
    else:
        mutation_data = patients_mut_genes().copy()

    if filename in PRESENT_FILES:
        return pd.read_csv(DATAFILES+filename, index_col=0)
    else:
        melanoma_kegg = pd.read_csv(DATAFILES+'Melanoma_kegg.csv')
        nodes_kegg = pd.concat([melanoma_kegg['src'], melanoma_kegg['dest']]).unique()
        binary_kegg = mut_genes_table(mutation_data, nodes_kegg)
        binary_kegg.to_csv(DATAFILES+filename, index=True)
        return binary_kegg
    
#map validation patients snp data to melanoma from kegg
# def map_melanoma_kegg_validation(filename = 'Binary_kegg_validation.csv'):
    
#     if filename in PRESENT_FILES:
#         return pd.read_csv(DATAFILES+filename, index_col=0)
#     else:
#         mutation_data = patients_mut_genes(type='validation').copy()
#         melanoma_kegg = pd.read_csv(DATAFILES+'Melanoma_kegg.csv')
#         nodes_kegg = pd.concat([melanoma_kegg['src'], melanoma_kegg['dest']]).unique()
#         binary_kegg = mut_genes_table(mutation_data, nodes_kegg)
#         binary_kegg.to_csv(DATAFILES+filename, index=True)
#         return binary_kegg

# map patients snp data to curated melanoma network from virtual cell
def map_melanoma_vcells(filename = 'Binary_vcells.csv', type='None'):
    if type=='validation':
        filename = 'Binary_vcells_validation.csv'
        mutation_data = patients_mut_genes(type='validation').copy()
    else:
        mutation_data = patients_mut_genes().copy()

    if filename in PRESENT_FILES:
        return pd.read_csv(DATAFILES+filename, index_col=0)
    else:
        melanoma_vcells = pd.read_csv(DATAFILES+'Melanoma_vcells.csv')
        nodes_vcells = pd.concat([melanoma_vcells['node1'], melanoma_vcells['node2']]).unique()
        binary_vcells = mut_genes_table(mutation_data, nodes_vcells)
        binary_vcells.to_csv(DATAFILES+filename, index=True)
        return binary_vcells
    

# # map validation patients snp data to curated melanoma network from virtual cell
# def map_melanoma_vcells_validation(filename = 'Binary_vcells_validation.csv'):
    
#     if filename in PRESENT_FILES:
#         return pd.read_csv(DATAFILES+filename, index_col=0)
#     else:
#         mutation_data = patients_mut_genes(type='validation').copy()
#         melanoma_vcells = pd.read_csv(DATAFILES+'Melanoma_vcells.csv')
#         nodes_vcells = pd.concat([melanoma_vcells['node1'], melanoma_vcells['node2']]).unique()
#         binary_vcells = mut_genes_table(mutation_data, nodes_vcells)
#         binary_vcells.to_csv(DATAFILES+filename, index=True)
#         return binary_vcells
    

############ HYPERGEOMETRIC TEST ############

# Perform hypergeometric test of mutated genes in all pathways
def perform_hypergeometric_test(patient_mutations, pathway_mapping, genes_population):
    # population of protein coding genes
    M = len(genes_population)

    patients_pvals = {}  # store info
    for i, (patient, pat_genes) in patient_mutations.iterrows():
        # mutated genes in patient
        #pat_genes = df_snps_comp[df_snps_comp.patient_id == patient].HGNC.unique()
        n = np.intersect1d(pat_genes, genes_population)
        
        pathway_pvals = {} # store pvals
        for pathway in pathway_mapping.pathway.unique():
            # total no. of genes in the pathway
            path_nodes = pathway_mapping[pathway_mapping.pathway == pathway].nodes.unique()
            N = np.intersect1d(path_nodes, genes_population)
            # no. of genes mutated in pathway
            #path_mut_genes = np.intersect1d(pat_genes, path_nodes)
            k = len(np.intersect1d(n, N))
            # perform hypergeometric test
            pval = hypergeom.pmf(k=k, M=M, n=len(n), N=len(N))
            value = -np.log10(pval)
            #lens = (len(pat_genes), len(n), len(path_nodes), len(N), len(path_mut_genes), len(k))
            #print(lens, pval, value)
            pathway_pvals[pathway] = value
        
        patients_pvals[patient] = pathway_pvals

    hg_results = pd.DataFrame(patients_pvals).T
    #print(hg_results.shape)
    return hg_results


# Calculate the pvalues for all pathways
def compute_pathway_pvalues(type=None):
    """ Perform hypergeometric test on mutated gene for all pathways and return pvalues """

    ## All protein coding genes fron HGNC database
    genes_population = pd.read_csv(DATAFILES+'hgnc_complete_set_2022-09-01.txt', sep='\t')
    genes_population = genes_population[genes_population.locus_group == 'protein-coding gene'].symbol
    #genes_population = genes_population[['symbol', 'alias_symbol', 'prev_symbol',]]

    ## Genes mapped to pathways
    # Load file with genes mapped to all pathways in kegg
    pathways_kegg = pd.read_csv(DATAFILES+'Pathways_kegg.csv')
    # Load file with genes mapped to Vcells melanoma pathway
    edges_vcells = pd.read_csv(DATAFILES+'Melanoma_vcells.csv')
    nodes_vcells = pd.concat([edges_vcells['node1'],edges_vcells['node2']]).unique()
    pathway_vcells = pd.DataFrame({'pathway': 'Melanoma Map (Curated)', 'nodes': nodes_vcells})
    # dataframe with genes mapped to all pathways (kegg + vcells)
    pathway_mapping = pd.concat([pathways_kegg, pathway_vcells], ignore_index=True)

    # Patients mutation data
    if type=='validation':
        patient_mutations = load_snps(type='validation')
        # only patients with complete sequencing information
        sources_to_keep = ['Pauline Blateau, Jerome Solassol', 'Baptiste Louveau, Samia Mourah'] 
        patient_mutations = patient_mutations.loc[patient_mutations['source'].isin(sources_to_keep)]
        # Remove unmutated HGNC
        patient_mutations = patient_mutations[patient_mutations['mutated']=='yes']
        patient_mutations = patient_mutations[['patient_id', 'HGNC']].drop_duplicates()
        patient_mutations = patient_mutations.groupby(['patient_id'], as_index=False).agg({'HGNC': lambda x: x.tolist()})
    else:
        patient_mutations = load_snps()
        # only patients with complete sequencing information
        patient_mutations = patient_mutations[patient_mutations['mutated']=='yes']
        patient_mutations = patient_mutations[patient_mutations.source != 'Pauline Blateau, Jerome Solassol']
        patient_mutations = patient_mutations[['patient_id', 'HGNC']].drop_duplicates()
        patient_mutations = patient_mutations.groupby(['patient_id'], as_index=False).agg({'HGNC': lambda x: x.tolist()})

    # # For Cscape analysis - considering only cscape genes for analysis
    # cscape_data = pd.read_csv(DATAFILES+'Binary_cscape_oncogenic.csv', index_col=0)
    # cscape_genes_df = pd.DataFrame()
    # cscape_genes_df['HGNC'] = cscape_data[cscape_data.columns].astype('bool').apply(lambda row: list(cscape_data.columns[row]), axis=1)
    # #cscape_genes_df = cscape_genes_df.explode('HGNC')
    # cscape_genes_df.reset_index(inplace=True)
    # cscape_genes_df.rename(columns={'index':'patient_id'}, inplace=True)

    ## Perform hypergeometric test for all pathways
    test_results = perform_hypergeometric_test(patient_mutations, pathway_mapping, genes_population)
    # test_results = perform_hypergeometric_test(cscape_genes_df, pathway_mapping, genes_population) # Cscape analysis

    ## Scale the values
    from sklearn.preprocessing import MinMaxScaler
    scaled_results = MinMaxScaler().fit_transform(test_results)
    scaled_results = pd.DataFrame(scaled_results, columns=test_results.columns, index=test_results.index)

    return scaled_results


# Load computed pathway pvalues
def get_pathway_pvalues(filename = 'Pathway_pvalues.csv', type = None):
    if type=='validation':
        filename='Pathway_pvalues_validation.csv'
        if filename in PRESENT_FILES:
            return pd.read_csv(DATAFILES+filename, index_col=0)
        else:
            print('computing pvalues ....')
            pval_df = compute_pathway_pvalues(type='validation')
            pval_df.to_csv(DATAFILES+filename, index=True)
    else:
    
        if filename in PRESENT_FILES:
            return pd.read_csv(DATAFILES+filename, index_col=0)
        else:
            pval_df = compute_pathway_pvalues()
            pval_df.to_csv(DATAFILES+filename, index=True)
        return pval_df



############ ADJACENCY MATRIX ############

# create adjacency matrix of a network
def network_adjacency(edgelist, source, target):
    """ Return adjacency matrix of network. """
    # process network
    edges_rev = pd.DataFrame({source : edgelist[target], target: edgelist[source]})
    edges = pd.concat([edgelist, edges_rev], ignore_index=True).drop_duplicates()
    # adjacency matrix
    network_dict = dict()
    for node in edges[source].unique():
        targets = edges[edges[source] == node][target].unique()
        if node not in network_dict:
            network_dict[node] = pd.Series([1 for _ in targets], index=targets)
    adj_mat = pd.DataFrame(network_dict).T.fillna(0)
    adj_mat = adj_mat.reindex(sorted(adj_mat.index))
    return adj_mat


# adjacency matrix of melanoma network from kegg
def adj_melanoma_kegg(filename = 'Melanoma_kegg_adj.csv'):
    
    if filename in PRESENT_FILES:
        return pd.read_csv(DATAFILES+filename, index_col=0)
    else:
        melanoma_kegg = pd.read_csv(DATAFILES+'Melanoma_kegg.csv')
        adj_mat_kegg = network_adjacency(melanoma_kegg, 'src', 'dest')
        adj_mat_kegg.to_csv(DATAFILES+filename, index=True)
        return adj_mat_kegg


# adjacency matrix of curated melanoma network from virtual cell
def adj_melanoma_vcells(filename = 'Melanoma_vcells_adj.csv'):
    
    if filename in PRESENT_FILES:
        return pd.read_csv(DATAFILES+filename, index_col=0)
    else:
        melanoma_vcells = pd.read_csv(DATAFILES+'Melanoma_vcells.csv')
        adj_mat_vcells = network_adjacency(melanoma_vcells, 'node1', 'node2')
        adj_mat_vcells.to_csv(DATAFILES+filename, index=True)
        return adj_mat_vcells



############ MUTATION LOAD ############

# function to find mutation load of given mutated genes
def pathway_mutation_load(mutated_genes: list, pathway_mapping: pd.DataFrame) -> dict:
    
    # store all unique pathways where either or both mutated genes are present
    pathways = set()
    for gene in mutated_genes:
        paths = pathway_mapping[pathway_mapping.nodes == gene].pathway
        pathways.update(paths)
    #print(len(pathways))

    # mutation load for each unique pathway
    mutation_load = {}
    for pathway in pathways:
        pathway_df = pathway_mapping[pathway_mapping.pathway == pathway]
        total_nodes = pathway_df.shape[0]
        total_mut_genes = np.intersect1d(mutated_genes, pathway_df.nodes.unique())
        mutation_rate = len(total_mut_genes) / total_nodes
        #print(f'{pathway}\t{total_nodes}\t{total_mut_genes}\t{len(total_mut_genes)}\t{mutation_rate}')
        mutation_load[pathway] = mutation_rate # store mutation rate
    
    return mutation_load


# calculate mutational load in patients for all pathways
def calculate_mutation_load():
    """ Calculates mutational load of mutated genes in patients for all pathways (kegg + vcells). """

    ## Genes mapped to pathways
    # Load file with genes mapped to all pathways in kegg
    pathways_kegg = pd.read_csv(DATAFILES+'Pathways_kegg.csv')
    # Load file with genes mapped to Vcells melanoma pathway
    edges_vcells = pd.read_csv(DATAFILES+'Melanoma_vcells.csv')
    nodes_vcells = pd.concat([edges_vcells['node1'],edges_vcells['node2']]).unique()
    pathway_vcells = pd.DataFrame({'pathway': 'Melanoma Map (Curated)', 'nodes': nodes_vcells})
    # dataframe with genes mapped to all pathways (kegg + vcells)
    pathways_genes = pd.concat([pathways_kegg, pathway_vcells], ignore_index=True)

    ## Mutation data of patients
    mutation_data = patients_mut_genes().copy()
    mutation_data = mutation_data.groupby(['patient_id'], as_index=False).agg({'HGNC': lambda x: x.tolist()})

    ## Mutational load for all pathways
    mutation_load = dict()
    load_copy = dict() # copy of mutation load of repeated set of mutated genes
    for i, (patient, mut_genes) in mutation_data.iterrows():
        if tuple(mut_genes) in load_copy:
            load = load_copy[tuple(mut_genes)]
        else:
            load = pathway_mutation_load(mut_genes, pathways_genes)
            load_copy[tuple(mut_genes)] = load
        #print(i)
        mutation_load[patient] = pd.Series(list(load.values()), index=list(load.keys()))
    
    # dataframe with patients mutation load in all pathways
    return pd.DataFrame(mutation_load).T


# load mutational load of patients
def get_mutational_load():
    filename = 'Mutational_load2.csv'
    if filename in PRESENT_FILES:
        return pd.read_csv(DATAFILES+filename, index_col=0)
    else:
        ml_df = calculate_mutation_load()
        ml_df.to_csv(DATAFILES+filename, index=True)
        return ml_df

