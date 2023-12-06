### Module to create and load different datasets
## additional features are added to clinical dataset



# import required packages
import pandas as pd
import numpy as np
from load_data import load_clinical, load_clinical_validation
import integrate_pathway_data as ipd



## Threshold to drop column with many missing values
def threshold_level(dataset, columns, percentage=5):
    """ Returns columns with missing values greater than threshold. """
    threshold = dataset.shape[0] * (percentage / 100)
    # number of filled values
    not_null = pd.Series(dataset[columns].notna().sum())
    #not_null = pd.Series((dataset[columns] == 1).sum())
    # columns to drop ( below threshold )
    drop_cols = not_null[not_null < threshold].index
    return drop_cols


## Prepare FATHMM and CSCAPE dataset

def dataset_from_binary(dataset_from_binary_file, thr):
    clinical_x_dataset = clinical_data.join(dataset_from_binary_file, on='patient_ID')
    # Drop empty columns
    to_drop= threshold_level(clinical_x_dataset, dataset_from_binary_file.columns, thr)
    clinical_x_dataset.drop(columns=to_drop, inplace=True)
    # fill na with 0
    cols_added = list(set(clinical_x_dataset.columns) - set(clinical_data.columns))
    clinical_x_dataset[cols_added] = clinical_x_dataset[cols_added].fillna(value=0)
    
    return(clinical_x_dataset)

########## CLINICAL DATASET ##########

# cleaned clinical data
clinical_data = load_clinical(binary_BRAF=False)
#clinical_data = load_clinical(binary_BRAF=False, drop_dcr=True)
#clinical_data = load_clinical(binary_BRAF=False, censor_val=12.0, drop_dcr=True)


########## ADDITIONAL FEATURES ##########

# # Melanoma network (KEGG)
# binary_kegg = ipd.map_melanoma_kegg()
# # Melanoma network (Virtual cell)
# binary_vcells = ipd.map_melanoma_vcells()
# # Pathway pvalues
# pathway_pvalues = ipd.get_pathway_pvalues()


# #### For Cscape Analysis (CA)
# # Melanoma network (KEGG)
# binary_kegg = ipd.map_melanoma_kegg('CA_Binary_kegg.csv')
# # Melanoma network (Virtual cell)
# binary_vcells = ipd.map_melanoma_vcells('CA_Binary_vcells.csv')
# # Pathway pvalues
# pathway_pvalues = ipd.get_pathway_pvalues('CA_Pathway_pvalues.csv')


# #### For Catalanotti Analysis (CAT)
# # Melanoma network (KEGG)
# binary_kegg = ipd.map_melanoma_kegg('CAT_Binary_kegg.csv')
# # Melanoma network (Virtual cell)
# binary_vcells = ipd.map_melanoma_vcells('CAT_Binary_vcells.csv')
# # Pathway pvalues
# pathway_pvalues = ipd.get_pathway_pvalues('CAT_Pathway_pvalues.csv')

#full_genes='ALL'
#full_genes="CA"
full_genes="CAT"

if full_genes == 'ALL':
    # Melanoma network (KEGG)
    binary_kegg = ipd.map_melanoma_kegg()
    # Melanoma network (Virtual cell)
    binary_vcells = ipd.map_melanoma_vcells()
    # Pathway pvalues
    pathway_pvalues = ipd.get_pathway_pvalues()
elif full_genes == 'CAT':
    #### For Catalanotti Analysis (CAT)
    # Melanoma network (KEGG)
    binary_kegg = ipd.map_melanoma_kegg('CAT_Binary_kegg.csv')
    # Melanoma network (Virtual cell)
    binary_vcells = ipd.map_melanoma_vcells('CAT_Binary_vcells.csv')
    # Pathway pvalues
    pathway_pvalues = ipd.get_pathway_pvalues('CAT_Pathway_pvalues.csv')
elif full_genes == 'CA':
    # #### For Cscape Analysis (CA)
    # Melanoma network (KEGG)
    binary_kegg = ipd.map_melanoma_kegg('CA_Binary_kegg.csv')
    # Melanoma network (Virtual cell)
    binary_vcells = ipd.map_melanoma_vcells('CA_Binary_vcells.csv')
    # Pathway pvalues
    pathway_pvalues = ipd.get_pathway_pvalues('CA_Pathway_pvalues.csv')

########## MELANOMA NETWORK (KEGG) ##########

# join clinical and network data on patient id
clinical_kegg = clinical_data.join(binary_kegg, on='patient_ID')
# # Drop empty columns
# #clinical_kegg.dropna(how='all', axis=1, inplace=True)
# drop_kegg = threshold_level(clinical_kegg, binary_kegg.columns)
# clinical_kegg.drop(columns=drop_kegg, inplace=True)
# # fill na with 0
# cols_added = list(set(clinical_kegg) - set(clinical_data))
clinical_kegg[binary_kegg.columns] = clinical_kegg[binary_kegg.columns].fillna(value=0)



########## MELANOMA NETWORK (VCELLS) ##########

# join clinical and network data on patient id
clinical_vcells = clinical_data.join(binary_vcells, on='patient_ID')
# # Drop empty columns
# #clinical_vcells.dropna(how='all', axis=1, inplace=True)
# drop_vcells = threshold_level(clinical_vcells, binary_vcells.columns)
# clinical_vcells.drop(columns=drop_vcells, inplace=True)
# # fill na with 0
# cols_added = list(set(clinical_vcells) - set(clinical_data))
clinical_vcells[binary_vcells.columns] = clinical_vcells[binary_vcells.columns].fillna(value=0)



########## PATHWAYS P-VALUE ##########

# join clinical and pathways pvalues on patient id
clinical_pp = clinical_data.join(pathway_pvalues, on='patient_ID')
# # fill na with 0
# cols_added = list(set(clinical_pp) - set(clinical_data))
clinical_pp[pathway_pvalues.columns] = clinical_pp[pathway_pvalues.columns].fillna(value=0)



########## MELANOMA PATHWAY P-VALUE (KEGG) ##########

# add pvalues of melanoma pathway from kegg
clinical_pp_kegg = clinical_data.join(pathway_pvalues['Melanoma'], on='patient_ID')
clinical_pp_kegg['Melanoma'] = clinical_pp_kegg['Melanoma'].fillna(value=0)



########## MELANOMA PATHWAY P-VALUE (VCELLS) ##########

# add pvalues of melanoma pathway from vcells
clinical_pp_vcells = clinical_data.join(pathway_pvalues['Melanoma Map (Curated)'], on='patient_ID')
clinical_pp_vcells['Melanoma Map (Curated)'] = clinical_pp_vcells['Melanoma Map (Curated)'].fillna(value=0)



# ########## MUTATIONAL LOAD ##########

# # join clinical and mutation load on patient id
# clinical_ml = clinical_data.join(mutation_load, on='patient_ID')
# # # Drop empty columns
# # #clinical_ml.dropna(how='all', axis=1, inplace=True)
# # drop_ml = threshold_level(clinical_ml, mutation_load.columns)
# # clinical_ml.drop(columns=drop_ml, inplace=True)
# # # fill na with 0
# # cols_added = list(set(clinical_ml) - set(clinical_data))
# clinical_ml[mutation_load.columns] = clinical_ml[mutation_load.columns].fillna(value=0)



# ########## MELANOMA MUTATIONAL LOAD (KEGG) ##########

# # add mutation load of melanoma pathway from kegg
# clinical_ml_kegg = clinical_data.join(mutation_load['Melanoma'], on='patient_ID')
# clinical_ml_kegg['Melanoma'] = clinical_ml_kegg['Melanoma'].fillna(value=0)



# ########## MELANOMA MUTATIONAL LOAD (VCELLS) ##########

# # add mutation load of melanoma pathway from vcells
# clinical_ml_vcells = clinical_data.join(mutation_load['Melanoma Map (Curated)'], on='patient_ID')
# clinical_ml_vcells['Melanoma Map (Curated)'] = clinical_ml_vcells['Melanoma Map (Curated)'].fillna(value=0)



########## COMBINED DATASET ##########

## Combine all additional features
# common genes between kegg and vcells network
common_genes = np.intersect1d(binary_kegg.columns, binary_vcells.columns)
# add kegg network data
combined_c_k = clinical_data.join(binary_kegg, on='patient_ID')
# add vcells network data ( dropped common genes )
combined_c_k_v = combined_c_k.join(binary_vcells.drop(columns=common_genes), on='patient_ID')
# add pathway pvalues
combined_c_k_v_p = combined_c_k_v.join(pathway_pvalues, on='patient_ID')
# # Drop empty columns
# #combined_c_k_v_p.dropna(how='all', axis=1, inplace=True)
# drop_cols = threshold_level(combined_c_k_v_p, combined_c_k_v_p.drop(columns=clinical_data.columns).columns)
# combined_c_k_v_p.drop(columns=drop_cols, inplace=True)
# fill na with 0
cols_added = list(set(combined_c_k_v_p) - set(clinical_data))
combined_c_k_v_p[cols_added] = combined_c_k_v_p[cols_added].fillna(value=0)



########## FATHMM DATA ##########

fathmm_data = pd.read_csv('../datafiles/Binary_fathmm_sarah.csv', index_col=0)
fathmm_data = fathmm_data.replace({0 : np.nan})
# join clinical and fathmm data on patient id
clinical_fathmm = clinical_data.join(fathmm_data, on='patient_ID')
# Drop empty columns
#clinical_fathmm.dropna(how='all', axis=1, inplace=True)
drop_fathmm = threshold_level(clinical_fathmm, fathmm_data.columns, 3)
clinical_fathmm.drop(columns=drop_fathmm, inplace=True)
# fill na with 0
cols_added = list(set(clinical_fathmm.columns) - set(clinical_data.columns))
clinical_fathmm[cols_added] = clinical_fathmm[cols_added].fillna(value=0)



########## CSCAPE HIGH DATA ##########

cscape_high_data = pd.read_csv('../datafiles/Binary_cscape_HIGH_oncogenic.csv', index_col=0)
cscape_high_data = cscape_high_data.replace({0 : np.nan})
# join clinical and cscape_high data on patient id
clinical_cscape_high = clinical_data.join(cscape_high_data, on='patient_ID')
# # Drop empty columns
# #clinical_cscape_high.dropna(how='all', axis=1, inplace=True)
# # drop 0%
# clinical_cscape_high0 = clinical_cscape_high.copy()
# # initial_drop
# drop_cscape_high = threshold_level(clinical_cscape_high, cscape_high_data.columns, 5)
# clinical_cscape_high.drop(columns=drop_cscape_high, inplace=True)
# # fill na with 0
# cols_added0 = list(set(clinical_cscape_high0.columns) - set(clinical_data.columns))
# clinical_cscape_high0[cols_added0] = clinical_cscape_high0[cols_added0].fillna(value=0)

cols_added = list(set(clinical_cscape_high.columns) - set(clinical_data.columns))
clinical_cscape_high[cols_added] = clinical_cscape_high[cols_added].fillna(value=0)


########## CSCAPE DATA ##########

cscape_data = pd.read_csv('../datafiles/Binary_cscape_oncogenic.csv', index_col=0)
cscape_data = cscape_data.replace({0 : np.nan})
# join clinical and cscape data on patient id
clinical_cscape = clinical_data.join(cscape_data, on='patient_ID')
# Drop empty columns
#clinical_cscape.dropna(how='all', axis=1, inplace=True)
#drop 3%
drop_cscape3 = threshold_level(clinical_cscape, cscape_data.columns, 3)
clinical_cscape3 = clinical_cscape.copy()
clinical_cscape3.drop(columns=drop_cscape3, inplace=True)
#drop 2%
drop_cscape2 = threshold_level(clinical_cscape, cscape_data.columns, 2)
clinical_cscape2 = clinical_cscape.copy()
clinical_cscape2.drop(columns=drop_cscape2, inplace=True)
#initial drop
drop_cscape = threshold_level(clinical_cscape, cscape_data.columns, 5)
clinical_cscape.drop(columns=drop_cscape, inplace=True)
# Fill na with 0 
cols_added3 = list(set(clinical_cscape3.columns) - set(clinical_data.columns))
clinical_cscape3[cols_added3] = clinical_cscape3[cols_added3].fillna(value=0)

cols_added2 = list(set(clinical_cscape2.columns) - set(clinical_data.columns))
clinical_cscape2[cols_added2] = clinical_cscape2[cols_added2].fillna(value=0)

cols_added = list(set(clinical_cscape.columns) - set(clinical_data.columns))
clinical_cscape[cols_added] = clinical_cscape[cols_added].fillna(value=0)


#############################################################
#################### VALIDATION DATASETS ####################
#############################################################


# ########## CLINICAL DATASET ##########

# cleaned clinical data
clinical_data_validation = load_clinical_validation(binary_BRAF=False)
#clinical_data = load_clinical_clean_validation(binary_BRAF=False, censor_val=6.0, drop_dcr=True)

# spns data
binary_kegg_validation = ipd.map_melanoma_kegg(type='validation')
binary_vcells_validation = ipd.map_melanoma_vcells(type='validation')
pathway_pvalues_validation = ipd.compute_pathway_pvalues(type='validation')

########## MELANOMA NETWORK (KEGG) ##########
# join clinical and network data on patient id
clinical_kegg_validation = clinical_data_validation.join(binary_kegg_validation, on='patient_ID')
clinical_kegg_validation[binary_kegg_validation.columns] = clinical_kegg_validation[binary_kegg_validation.columns].fillna(value=0)


# ########## MELANOMA NETWORK (VCELLS) ##########
# # join clinical and network data on patient id
clinical_vcells_validation = clinical_data_validation.join(binary_vcells_validation, on='patient_ID')
clinical_vcells_validation[binary_vcells_validation.columns] = clinical_vcells_validation[binary_vcells_validation.columns].fillna(value=0)


# ########## PATHWAYS P-VALUE ##########
# # join clinical and pathways pvalues on patient id
clinical_pp_validation = clinical_data_validation.join(pathway_pvalues_validation, on='patient_ID')
clinical_pp_validation[pathway_pvalues_validation.columns] = clinical_pp_validation[pathway_pvalues_validation.columns].fillna(value=0)

########## MELANOMA PATHWAY P-VALUE (KEGG) ##########
# add pvalues of melanoma pathway from kegg
clinical_pp_kegg_validation = clinical_data_validation.join(pathway_pvalues_validation['Melanoma'], on='patient_ID')
clinical_pp_kegg_validation['Melanoma'] = clinical_pp_kegg_validation['Melanoma'].fillna(value=0)


########## MELANOMA PATHWAY P-VALUE (VCELLS) ##########

# add pvalues of melanoma pathway from vcells
clinical_pp_vcells_validation = clinical_data_validation.join(pathway_pvalues_validation['Melanoma Map (Curated)'], on='patient_ID')
clinical_pp_vcells_validation['Melanoma Map (Curated)'] = clinical_pp_vcells_validation['Melanoma Map (Curated)'].fillna(value=0)


########## LOAD ALL DATASETS ##########

def load_datasets(type=None):
    """ Load all datasets. """
    if type=='validation':
        datasets = {
        'clinical_val' : clinical_data_validation,
        'clinical_kegg_val' : clinical_kegg_validation,
        'clinical_vcells_val' : clinical_vcells_validation,
        'clinical_pp_val' : clinical_pp_validation,
        'clinical_pp_kegg_val' : clinical_pp_kegg_validation,
        'clinical_pp_vcells_val' : clinical_pp_vcells_validation,
    }
    else:
        datasets = {
            'clinical' : clinical_data,
            'clinical_kegg' : clinical_kegg,
            'clinical_vcells' : clinical_vcells,
            'clinical_pp' : clinical_pp,
            'clinical_pp_kegg' : clinical_pp_kegg,
            'clinical_pp_vcells' : clinical_pp_vcells,
            'combined_data' : combined_c_k_v_p,
            'clinical_fathmm' : clinical_fathmm,
            'clinical_cscape' : clinical_cscape,
            'clinical_cscape_high' : clinical_cscape_high
        }
    return datasets