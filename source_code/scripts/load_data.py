### Module to load the clinical and mutation datasets


# Import required packages 
import sqlite3
import pandas as pd
import numpy as np
import json
import os


# Import database 

cnx = sqlite3.connect('../datafiles/melanomadb_original.db')
cnx_validation = sqlite3.connect('../datafiles/melanomadb_validation.db')

# Load raw clinical dataset
def load_clinical_raw(type = None):
    if type=='validation':
        df_patients = pd.read_sql_query("SELECT * FROM patients", cnx_validation)
    else:
        df_patients = pd.read_sql_query("SELECT * FROM patients", cnx)
    # drop irrelevant columns
    df_patients.drop(columns=['id', 'data', 'creation_datetime', 'type'], inplace=True)
    df_patients['source'] = df_patients['source'].apply(lambda x : json.loads(x)['author']) # not needed
    # # remove leukemia patient
    # df_patients.drop(index=df_patients[df_patients.drug == 'vem+cob (leukemia patient)'].index, inplace=True)
    # # remove whitespace in categories
    # df_patients['disease_control_rate'] = df_patients['disease_control_rate'].apply(lambda x: x.strip() if x is not None else x)
    # unify ambiguous missing values
    df_patients.replace(to_replace=['<NA>', 'nan', 'N.E.', None], value=np.nan, inplace=True)
    df_patients.reset_index(drop=True, inplace=True)
    return df_patients


# Split drug categories
def split_drug_categories(type = None):
    """ Split drug categories to mono and combined drug therapy categories. """
    if type == 'validation':
        drug_pats = load_clinical_raw(type='validation')[['patient_ID', 'drug']].copy()
    else:
        drug_pats = load_clinical_raw()[['patient_ID', 'drug']].copy()
    unique_drugs = ['vemurafenib', 'cobimetinib', 'dabrafenib', 'trametinib']
    drug_cats = dict()
    for _, (id, pat_drugs) in drug_pats.iterrows():
        if pat_drugs is np.nan:
            # pass
            drug_cats[id] = pd.Series([0 for _ in unique_drugs], index=unique_drugs)
        elif ' + ' in pat_drugs:
            drugs = pat_drugs.split(' + ')
            drug_cats[id] = pd.Series([1 for _ in drugs], index=drugs)
        else:
            drug_cats[id] = pd.Series([1], index=[pat_drugs])
    
    return pd.DataFrame(drug_cats).T.fillna(0)


# Additional features
def additional_features():
    """ Returns additional features of patients """
    # load and process raw data
    pat_features = pd.read_sql_query("SELECT * FROM patients", cnx)
    pat_features['source'] = pat_features['source'].apply(lambda x : json.loads(x)['author'])
    pat_features = pat_features.replace(to_replace=['<NA>', 'nan', 'N.E.', None], value=np.nan)
    pat_features = pat_features[['patient_ID', 'drug', 'source']]

    # feature for available snp data; 1 if snp data available
    pat_features['snp_data'] = pat_features.source.map({'Yibing Yan, Antoni Ribas' : 0}).fillna(1)

    # feature for type of genome sequencing; 1 if complete genome sequencing was done
    pat_features['complete_sequencing'] = pat_features.source.map({
        'Eliezer M. Van Allen, Dirk Schadendorf' : 1,
        'Federica Catalanotti, David B. Solit' : 1}).fillna(0)

    # # feature for type of drug therapy; 1 if combination therapy given
    # pat_features['comb_therapy'] = pat_features.drug.map({'vemurafenib + cobimetinib' : 1, 'dabrafenib + trametinib' : 1}).fillna(0)

    # # feature for type of drug therapy; 1 if mono therapy given
    # pat_features['mono_therapy'] = pat_features.drug.map({'vemurafenib' : 1, 'dabrafenib' : 1}).fillna(0) # high correlation

    pat_features.drop(columns=['drug', 'source'], inplace=True)
    pat_features.set_index('patient_ID', inplace=True)

    return pat_features


# Load cleaned clinical dataset
def load_clinical_clean(binary_BRAF = True, drop_BRAF_mut = True, censor_val=None, drop_dcr=False):
    df_patients = load_clinical_raw()
    # Set PFS_statut of van_allen patients to 1
    df_patients.loc[(df_patients.source == 'Eliezer M. Van Allen, Dirk Schadendorf'), 'PFS_statut']=1
    # Drop OS_month and OS_statut (not required for prediction of PFS)
    # Drop brain_metastasis, immunotherapy_treatment, M_stage ( many missing values)
    df_patients.drop(columns=['M_stage', 'brain_metastasis', 'immunotherapy_treatment', 'OS_statut', 'OS_month'], inplace=True)
    
    ## clinical data with source
    #df_patients.drop(columns=['source'], inplace=True)
    source_cat = {
        'Federica Catalanotti, David B. Solit' : 1.0,
        'Eliezer M. Van Allen, Dirk Schadendorf' : 0.0,
        }
    df_patients['source'] = df_patients.source.map(source_cat)

    # Remove missing values in PFS_statut
    df_patients.dropna(subset=['PFS_statut'], inplace=True) # (2 rows)
    # change str dtype events to numeric 
    df_patients['PFS_statut'] = pd.to_numeric(df_patients.PFS_statut)
    # # Drop rows with AJCC_stage, BRAF_mut, LDH missing
    # df_patients.dropna(subset=['AJCC_stage'], inplace=True) #(45 rows)
    # # Drop missing values in drug and disease_control_rate
    # df_patients.dropna(subset=['drug', 'disease_control_rate'], inplace=True) # (1+4 rows) # 'nan' and 'N.E.'

    # # harmonize categories in disease control rate
    # ## Pat_27 : PR ; Pat_55 : SD
    # df_patients.replace(to_replace={'PR/SD' : 'PR', 'SD/PR' : 'SD'}, inplace=True)

    # patients with complete sequencing data
    df_patients = df_patients[df_patients.sequencing_type == 'complete']
    #df_patients = df_patients[df_patients.disease_control_rate != 'PD']
    df_patients.drop(columns=['sequencing_data', 'sequencing_type'], inplace=True)
    # remove patients with only post treatment snp data (6 rows)
    df_patients = df_patients[df_patients.patient_ID.isin(['Pat_02', 'Pat_21', 'Pat_27', 'Pat_28', 'Pat_36', 'Pat_37', 'Pat_49']) == False]

    # drop feature AJCC_stage - high class imbalance
    df_patients.drop(columns=['AJCC_stage'], inplace=True)

    # drop feature disease control rate (DCR) - for comparison with DCR dataset
    #df_patients.drop(columns=['disease_control_rate'], inplace=True)

    # # patients with no mutation data removed (19 rows)
    # df_patients = df_patients[df_patients.BRAF_mut != 'MND']

    # split drug categories
    drug_table = split_drug_categories()
    df_patients = df_patients.join(drug_table, on='patient_ID')
    df_patients.drop(columns=['drug'], inplace=True)

    # # add additional features
    # add_features = additional_features()
    # df_patients = df_patients.join(add_features, on='patient_ID')

    # # remove patients with no snp data (130 Yan study patients)
    # df_patients = df_patients[df_patients.snp_data == 1]
    # df_patients.drop(columns=['snp_data'], inplace=True)

    # # remove patients with only post treatment snp data (6 rows)
    # df_patients = df_patients[df_patients.patient_ID.isin(['Pat_02', 'Pat_21', 'Pat_27', 'Pat_28', 'Pat_36', 'Pat_37', 'Pat_49']) == False]

    # change censoring value
    if censor_val:
        df_patients.loc[df_patients['PFS_month']>censor_val, 'PFS_statut']=0
        df_patients.loc[df_patients['PFS_month']>censor_val, 'PFS_month']=censor_val

    #drop_dcr
    if drop_dcr==True:
        df_patients.drop(columns=['disease_control_rate'], inplace=True)

    # clean dataset
    df_patients.reset_index(drop=True, inplace=True)
    # # BRAF gene binary column from BRAF data
    # if binary_BRAF == True:
    #     df_patients['BRAF'] = df_patients['BRAF_mut'].map({'MND' : 0}).fillna(1)
    # Drop BRAF mutation feature
    if drop_BRAF_mut == True:
        df_patients.drop(columns=['BRAF_mut'], inplace=True)
    
    return df_patients

# Load clinical clean for Blateu & Louveau patients
def load_clinical_clean_validation(binary_BRAF = True, drop_BRAF_mut = True, censor_val=None, drop_dcr=False):
    df_patients = load_clinical_raw(type='validation')
    sources_to_keep = ['Pauline Blateau, Jerome Solassol', 'Baptiste Louveau, Samia Mourah']
    df_patients = df_patients.loc[df_patients['source'].isin(sources_to_keep)]
    #########################################################
    # remove non valid PFS month patients from Blateau (after curation we choose to remove some patients with non coherent PFS month values)
    #pd_patients = df_patients[(df_patients['source']=='Pauline Blateau, Jerome Solassol') & (df_patients['disease_control_rate']=='PD')]
    #to_remove = pd_patients[pd_patients['PFS_month']>5.2]
    #pr_pat_to_remove = df_patients[df_patients['patient_ID']=='BS050']
    #df_patients = df_patients.drop(to_remove.index)
    #df_patients.loc[to_remove.index,'disease_control_rate']=np.nan
    #df_patients = df_patients.drop(pr_pat_to_remove.index)
    #invert PR & PD patients
    pr_pat = df_patients[(df_patients['source']=='Pauline Blateau, Jerome Solassol') & (df_patients['disease_control_rate']=='PR')]
    pd_pat = df_patients[(df_patients['source']=='Pauline Blateau, Jerome Solassol') & (df_patients['disease_control_rate']=='PD')]
    #df_patients.loc[pr_pat.index,'disease_control_rate']='PD'
    #df_patients.loc[pd_pat.index,'disease_control_rate']='PR'
    ##########################################################
    df_patients.drop(columns=['M_stage', 'brain_metastasis', 'immunotherapy_treatment', 'OS_statut', 'OS_month', 'source'], inplace=True)
    df_patients.dropna(subset=['PFS_statut'], inplace=True) # (2 rows)
    df_patients['PFS_statut'] = pd.to_numeric(df_patients.PFS_statut)
    # drop sequencing data columns & feature AJCC_stage - high class imbalance
    df_patients.drop(columns=['sequencing_data', 'sequencing_type','AJCC_stage'], inplace=True)
    # split drug category
    drug_table = split_drug_categories(type='validation')
    df_patients = df_patients.join(drug_table, on='patient_ID')
    df_patients.drop(columns=['drug'], inplace=True)
    # change censoring value
    if censor_val:
        df_patients.loc[df_patients['PFS_month']>censor_val, 'PFS_statut']=0
        df_patients.loc[df_patients['PFS_month']>censor_val, 'PFS_month']=censor_val
    #drop_dcr
    if drop_dcr==True:
        df_patients.drop(columns=['disease_control_rate'], inplace=True)
    # clean dataset
    df_patients.reset_index(drop=True, inplace=True)
    if drop_BRAF_mut == True:
        df_patients.drop(columns=['BRAF_mut'], inplace=True)
    return df_patients

# Load clinical dataset
def load_clinical(raw: bool = False, binary_BRAF = True, drop_BRAF_mut = True, censor_val = None, drop_dcr =False):
    if raw == True:
        return load_clinical_raw()
    else:
        return load_clinical_clean(binary_BRAF=binary_BRAF, drop_BRAF_mut=drop_BRAF_mut, censor_val = censor_val, drop_dcr=drop_dcr) 

#Load clinical_validation_dataset

def load_clinical_validation(raw: bool = False, binary_BRAF = True, drop_BRAF_mut = True, censor_val = None, drop_dcr =False):
    return load_clinical_clean_validation(binary_BRAF=binary_BRAF, drop_BRAF_mut=drop_BRAF_mut, censor_val = censor_val, drop_dcr=drop_dcr) 

# Load raw snps dataset
def load_snps(type = None):
    if type=='validation':
        df_snps = pd.read_sql_query("SELECT * FROM snps", cnx_validation)
    else:
        df_snps = pd.read_sql_query("SELECT * FROM snps", cnx)
    # get list of Catalanotti's sequenced genes
    bait = pd.ExcelFile("../datafiles/catalanotti_supplement2.xlsx").parse(0)
    list_genes = list(bait['Gene Symbol'])
    df_snps = df_snps[df_snps.HGNC.isin(list_genes)]
    # extract source
    import ast
    df_snps['source'] = df_snps['source'].apply(lambda x : ast.literal_eval(x)['author'])

    # drop irrelevant columns
    df_snps.drop(columns=['id', 'data', 'creation_datetime', 'other_prelevements'], inplace=True)
    #df_snps.drop(columns=['source'], inplace=True)
    
    # patients with pre-treatment snp data
    df_snps = df_snps[df_snps.temporality == 'pre treatment']

    # fill missing values
    snp_fills = {'consequence' : 'missing', 'variant_classification' : 'missing'}
    df_snps.fillna(value=snp_fills, inplace=True)
    df_snps.drop_duplicates(inplace=True)
    df_snps.reset_index(drop=True, inplace=True)
    return df_snps

