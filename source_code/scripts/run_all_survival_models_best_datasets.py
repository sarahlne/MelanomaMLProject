""" Main script for running the survival models on best and comparable datasets"""

######################################################
########## Import required external packages #########
######################################################

import sqlite3
import os
import sys
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from datetime import date, datetime
from timeit import default_timer
import time
from peewee import SqliteDatabase


######################################################
########## Internal imports               ############
######################################################

from datasets import load_datasets
import preprocessing_data as prep
from integrate_pathway_data import adj_melanoma_vcells
from integrate_pathway_data import adj_melanoma_kegg
from melanoma_survival_model import SurvivalModel

###############################################################################
########## #Import required packages for running survival models   ############
###############################################################################

# from sklearn
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder, StandardScaler, MinMaxScaler
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold, KFold
from sklearn.ensemble import RandomForestClassifier

# from sksurv

from sksurv.metrics import concordance_index_censored as cindex_h
from sksurv.metrics import concordance_index_ipcw as cindex_u
from sksurv.metrics import cumulative_dynamic_auc, integrated_brier_score, brier_score
from sklearn.metrics import plot_precision_recall_curve, plot_roc_curve
from sksurv.column import encode_categorical
from sksurv.metrics import concordance_index_censored
from sksurv.svm import FastSurvivalSVM
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis, ExtraSurvivalTrees
from sksurv.linear_model import CoxnetSurvivalAnalysis
#from coxnet_model import Coxnet
from math import ceil, floor
from sksurv.util import Surv
from sklearn.model_selection import GridSearchCV, KFold


#######################################################
##########  Folder names and log setting   ############
#######################################################

results_folder='../../source_code/model_results/'
sub_folder = 'main_result_pp_no_dcr_12_CAT'
results_path = results_folder+sub_folder+'/'
if not os.path.exists(results_path):
    os.makedirs(results_path)


# log initiation
log_file = results_path+"running_info.log"
if not os.path.exists(log_file):
    #os.mkdir(log_file)
    logging.basicConfig(filename=log_file, encoding='utf-8', level=logging.INFO)
logging.getLogger().setLevel(logging.INFO)
fh = logging.FileHandler(log_file, mode='w')
fh.setLevel(logging.INFO)
logger = logging.getLogger(__name__)
formatter = logging.Formatter("%(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)


##############################################################
##########  Datasets import & Models definition   ############
##############################################################


all_datasets = load_datasets()
# remove unusefull datasets
datasets_to_keep = ['clinical', 'clinical_pp', 'clinical_pp_kegg', 'clinical_pp_vcells']
all_datasets = {key: all_datasets[key] for key in datasets_to_keep}
all_models = {}
models_infos = pd.DataFrame(columns=['MODEL', 'PARAMS'])


######## Survival Random Forest
MODEL = RandomSurvivalForest(random_state=1)
MODEL_PARAMS = {
    'n_estimators' : [25, 50, 100, 150, 250],
    'max_depth' : [5, 15, 25],
}
MODEL_NAME="RSF"
rsf_model = SurvivalModel(model=MODEL, 
                          hyperparameters=MODEL_PARAMS, 
                          model_name=MODEL_NAME, 
                          datasets=all_datasets, 
                          save_models=True)
#rsf_model.results_folder= rsf_model.results_folder+sub_folder
rsf_model.results_folder=results_path
all_models['rsf'] = rsf_model
models_infos.loc[MODEL_NAME] = [MODEL, MODEL_PARAMS]

######## Gradient Boosting survival
MODEL = GradientBoostingSurvivalAnalysis(random_state=1)
MODEL_PARAMS = {
    'n_estimators' : [100, 200, 300], # 50 200 500
    #'max_depth' : [1, 2, 3],
    'learning_rate' : [0.01, 0.1, 0.5, 0.9] # 0.7, 0.9, , 0.05
}
MODEL_NAME = "GBS"
gbs_model = SurvivalModel(model=MODEL, 
                          hyperparameters=MODEL_PARAMS, 
                          model_name=MODEL_NAME, 
                          datasets=all_datasets, 
                          save_models=True)
gbs_model.results_folder= results_path
all_models['gbs'] = gbs_model
models_infos.loc[MODEL_NAME] = [MODEL, MODEL_PARAMS]

######## Extra Survival Trees
MODEL = ExtraSurvivalTrees(random_state=1)
MODEL_PARAMS = {
    'n_estimators' : [25, 50, 100, 200, 400],
    'max_depth' : [5, 15, 25],
}

MODEL_NAME = "EST"
est_model = SurvivalModel(model=MODEL, 
                          hyperparameters=MODEL_PARAMS, 
                          model_name=MODEL_NAME, 
                          datasets=all_datasets, 
                          save_models=True)
est_model.results_folder = results_path
all_models['est'] = est_model
models_infos.loc[MODEL_NAME] = [MODEL, MODEL_PARAMS]

######## COX Survival
MODEL = CoxnetSurvivalAnalysis(fit_baseline_model=True)
MODEL_PARAMS = {
    # 'n_alphas' : [25, 50, 100],
    # 'l1_ratio' : [0.001, 0.25, 0.5, 0.75, 1],
    'l1_ratio' : [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
}
MODEL_NAME = "COX"
cox_model = SurvivalModel(model=MODEL, 
                          hyperparameters=MODEL_PARAMS, 
                          model_name=MODEL_NAME, 
                          datasets=all_datasets, 
                          save_models=True)
cox_model.results_folder = results_path
all_models['cox'] = cox_model
models_infos.loc[MODEL_NAME] = [MODEL, MODEL_PARAMS]

######## Survival SVM
MODEL = FastSurvivalSVM(rank_ratio = 1.0, max_iter=1000, tol=1e-5, random_state=0)
MODEL_PARAMS={
    'alpha': 2. ** np.arange(-12,13,2)
}
MODEL_NAME="F-SVM"
svm_model = SurvivalModel(model=MODEL, 
                          hyperparameters=MODEL_PARAMS, 
                          model_name=MODEL_NAME, 
                          datasets=all_datasets, 
                          save_models=True)
svm_model.results_folder= results_path
all_models['svm'] = svm_model
models_infos.loc[MODEL_NAME] = [MODEL, MODEL_PARAMS]


# Register models parameters
models_infos.to_csv(results_folder+sub_folder+'/'+sub_folder+'.txt', sep = "\t", header=True)

##############################################################
##########  Running models   #################################
##############################################################


if __name__=='__main__':
    logger.info(f'{datetime.now().strftime("%d.%B %Y %H:%M")} Running survival models...results in {results_folder+sub_folder}:')
    start_time = time.time()
    for model in all_models:
        logger.info(f'Running {all_models[model].model_name} ...')
        print(all_models[model])
        curr_time = time.time()
        all_models[model].run_survival_model(evaluation=False)
        elapsed_time = time.time() - curr_time
        logger.info(f'...end of {all_models[model].model_name} training in {elapsed_time/60} min')
    elapsed_time = time.time() - start_time
    logger.info(f'Running complete in {elapsed_time/60} min')