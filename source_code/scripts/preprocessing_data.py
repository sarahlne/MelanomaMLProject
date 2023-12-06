### Module containing functions for preprocessing data and nested cross-validation.


# Import required packages 
import pandas as pd
import numpy as np

from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.preprocessing import OrdinalEncoder, OneHotEncoder, MinMaxScaler
from sklearn.ensemble import RandomForestClassifier

from sksurv.util import Surv
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold
from sksurv.metrics import concordance_index_ipcw as cindex_u
from sksurv.metrics import cumulative_dynamic_auc, integrated_brier_score, brier_score
from math import ceil, floor


# Variables

def get_cat(drop_dcr=False):
	if drop_dcr==True:
		CAT_COLS = ['sex', 'LDH'] #  , 'drug', 'AJCC_stage', 'disease_control_rate'
		NON_CAT_COLS = ['age']
		# All column-wise categories
		CATEGORIES = [
            ['female', 'male'],
            #['III', 'IV'],
            ['elevated', 'normal'],
            #['dabrafenib', 'dabrafenib + trametinib', 'vemurafenib', 'vemurafenib + cobimetinib'],
            #['CR', 'SD', 'PD', 'PR']
        ]
	else:
		CAT_COLS = ['sex', 'LDH','disease_control_rate'] #  ,'drug', 'AJCC_stage',
		NON_CAT_COLS = ['age']
		# All column-wise categories
		CATEGORIES = [
			['female', 'male'],
			#['III', 'IV'],
			['elevated', 'normal'],
			#['dabrafenib', 'dabrafenib + trametinib', 'vemurafenib', 'vemurafenib + cobimetinib'],
			['CR', 'SD', 'PD', 'PR']
        ]
	return CAT_COLS, NON_CAT_COLS, CATEGORIES

CAT_COLS, NON_CAT_COLS, CATEGORIES = get_cat(drop_dcr=False)
#CAT_COLS, NON_CAT_COLS, CATEGORIES = get_cat(drop_dcr=True)

#CAT_COLS = ['sex', 'LDH', 'disease_control_rate'] #  , 'drug', 'AJCC_stage', 'disease_control_rate'
#NON_CAT_COLS = ['age']
# All column-wise categories
#CATEGORIES = [
#    ['female', 'male'],
#    #['III', 'IV'],
#    ['elevated', 'normal'],
#    #['dabrafenib', 'dabrafenib + trametinib', 'vemurafenib', 'vemurafenib + cobimetinib'],
#    ['CR', 'SD', 'PD', 'PR']
#]


# # add BRAF mut col
# from load_data import BRAF_MUT_COL
# if BRAF_MUT_COL == True:
#     CAT_COLS.append('BRAF_mut')


# Seperate features and target
def seperate_x_y(data: pd.DataFrame):
    data = data.copy()
    # features
    X = data.drop(columns=['PFS_statut', 'PFS_month', 'patient_ID'])
    #X.set_index('patient_ID', inplace=True)
    # target
    # y = data[['PFS_statut', 'PFS_month']].copy(deep=True)
    # y['PFS_statut'] = y['PFS_statut'].map({'1': True, '0': False})
    # y = [tuple(i) for i in y.to_numpy()]
    # y = np.array(y, dtype=[('status', '?'), ('time_in_months', '<f8')])
    y = Surv.from_arrays(event=data['PFS_statut'], time=data['PFS_month'])

    return X, y



################### IMPUTATION ###################

# Imputation using Iterative Imputer

def impute_train_test(X_train: pd.DataFrame, X_test: pd.DataFrame):
    """ impute missing data using Iterative Imputer to fit train set and transform test set."""

    train = X_train.reset_index(drop=True)
    test = X_test.reset_index(drop=True)
    # Encode categorical values to numerical for imputer
    ord = OrdinalEncoder(categories=CATEGORIES, handle_unknown='use_encoded_value', unknown_value=np.nan)
    ord_fit = ord.fit(train[CAT_COLS]) # fit train
    # transform train
    train[CAT_COLS] = ord_fit.transform(train[CAT_COLS])
    # transform test
    test[CAT_COLS] = ord_fit.transform(test[CAT_COLS])
    # train imputer (on clinical data only)
    iimp = IterativeImputer(initial_strategy='most_frequent', estimator=RandomForestClassifier(random_state=1), max_iter=20)
    clinical_cols = CAT_COLS + NON_CAT_COLS
    train_clinical = train[clinical_cols]
    iimp_fit = iimp.fit(train_clinical) # fit train on clinical data
    # Impute missing values
    # transform train
    train_imp_cl= pd.DataFrame(iimp_fit.transform(train_clinical), columns=train_clinical.columns)
    train_imp = pd.concat([train_imp_cl, train.drop(columns=clinical_cols)], axis=1) # rest of data
    train_imp[CAT_COLS] = ord_fit.inverse_transform(train_imp[CAT_COLS])
    # transform test
    test_clinical = test[clinical_cols]
    test_imp_cl= pd.DataFrame(iimp_fit.transform(test_clinical), columns=test_clinical.columns)
    test_imp = pd.concat([test_imp_cl, test.drop(columns=clinical_cols)], axis=1) # rest of data
    test_imp[CAT_COLS] = ord_fit.inverse_transform(test_imp[CAT_COLS])

    return train_imp, test_imp



################### PREPROCESSING ###################

# Preprocessing using OneHot encoder and MinMax scaler

def preprocess_train_test(X_train: pd.DataFrame, X_test: pd.DataFrame):
    """Preprocessing using onehot encoding and scaling to fit train set and transform test set."""

    train = X_train.reset_index(drop=True)
    test = X_test.reset_index(drop=True)

    # Preprocessing step 1 : OneHot encoding
    onehot = OneHotEncoder(sparse=False, drop="first", categories = CATEGORIES)
    onehot_fit = onehot.fit(train[CAT_COLS]) # fit train categorical columns
    # transform train and concat to non categorical columns
    train_oh= pd.DataFrame(onehot_fit.transform(train[CAT_COLS]), columns=onehot_fit.get_feature_names_out())
    train_pp = pd.concat([train_oh, train.drop(columns=CAT_COLS)], axis=1)
    # transform test and concat to non categorical columns
    test_oh= pd.DataFrame(onehot_fit.transform(test[CAT_COLS]), columns=onehot_fit.get_feature_names_out())
    test_pp = pd.concat([test_oh, test.drop(columns=CAT_COLS)], axis=1)

    # Preprocessing step 2 : Scaling
    scaler = MinMaxScaler()
    scaler_fit = scaler.fit(train_pp['age'].values.reshape(-1, 1)) # fit train numeric column
    # transform train
    train_pp['age'] = scaler_fit.transform(train_pp['age'].values.reshape(-1, 1))
    # transform test
    test_pp['age'] = scaler_fit.transform(test_pp['age'].values.reshape(-1, 1))

    return train_pp, test_pp



################### NESTED CROSS-VALIDATION ###################

# Function to perform nested cross-validation given the features and 
# target variables, estimator and parameters

def perform_nested_cv(X, y, estimator, params, outer_splits=5, inner_splits=3, seed=1):
    """ Performs nested cv"""
    # store outer fold best models
    of_models = []
	# outer loop of the cross-validation procedure
    #cv_outer = KFold(n_splits=outer_splits, shuffle=True, random_state=seed)
    cv_outer = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=seed)
    for i, (train_ix, test_ix) in enumerate(cv_outer.split(X, X.source), 1):
		# split data
        X_train, X_test = X.iloc[train_ix, :], X.iloc[test_ix, :]
        y_train, y_test = y[train_ix], y[test_ix]

        # remove source feature
        feature_source = X_train.source.reset_index(drop=True)
        X_train = X_train.drop(columns='source')
        X_test = X_test.drop(columns='source')

		# impute and preprocess data
        X_train, X_test = impute_train_test(X_train, X_test)
        X_train, X_test = preprocess_train_test(X_train, X_test)

		# hyperparameter tuning inner loop
        #cv_inner = KFold(n_splits=inner_splits, shuffle=True, random_state=i)
        cv_inner = StratifiedKFold(n_splits=inner_splits, shuffle=True, random_state=i)
		# perform model selection
        search = GridSearchCV(estimator, params, cv=cv_inner.split(X_train, feature_source))
        search.fit(X_train, y_train)
        best_model = search.best_estimator_
        #print(best_model)

        # store cv best model for evaluation
        #train_test_set = (X_train, X_test, y_train, y_test) # for evaluation
        of_best_model = [i, best_model, X_train, y_train, X_test, y_test]
        of_models.append(of_best_model)
    
    return of_models



################### EVALUATION ###################

# time-intervals for time-dependent evaluation
def time_intervals(y_train, y_test):
    """ Returns time points for time-dependent evaluation"""
    # https://stackoverflow.com/questions/69439777/integrated-brier-score-for-sklearns-gridsearchcv
    y_times_tr = [i[1] for i in y_train]
    y_times_te = [i[1] for i in y_test]

    T1 = np.percentile(y_times_tr, 15, method='higher')
    T2 = np.percentile(y_times_tr, 94, method='lower')
    T3 = np.percentile(y_times_te, 15, method='higher')
    T4 = np.percentile(y_times_te, 94, method='lower')
    #print(T1, T2, T3, T4)
    return np.maximum(T1,T3), np.minimum(T2, T4)


# Evaluation of outer fold models

def evaluation(best_model, X_train, y_train, X_test, y_test):
    """ Evaluated survival models """

    ## Uno'c-index
    yhat_test = best_model.predict(X_test)
    uno_idx = cindex_u(y_train,y_test, yhat_test)[0] # Uno's c-index

    ## time-dependent scores
    t_min, t_max = time_intervals(y_train, y_test)
    event_times = np.arange(ceil(t_min), floor(t_max)+1, step=1, dtype=int)
    #print(t_min, t_max, event_times)
    #print(str(best_model.event_times_))
    # surv_func = best_model.predict_cumulative_hazard_function(X_test)
    # risk_scores = np.row_stack([func(event_times) for func in surv_func])

    ## cumulative hazard function
    chf = best_model.predict_cumulative_hazard_function(X_test)
    chf_scores = np.row_stack([func(event_times) for func in chf])

    ## survival function
    sf = best_model.predict_survival_function(X_test)
    sf_scores = np.row_stack([func(event_times) for func in sf])

    ## dynamic AUC
    td_auc = cumulative_dynamic_auc(y_train, y_test, chf_scores, event_times)[0]
    dynamic_auc = list(zip(event_times, td_auc))
    mean_auc = np.nanmean(td_auc)

    ## brier score
    _, td_bs = brier_score(y_train, y_test, sf_scores, event_times)
    brier_scores = list(zip(event_times, td_bs))

    ## integrated brier score
    ibs = integrated_brier_score(y_train, y_test, sf_scores, event_times)

    ## kaplan meier curve
    times = np.row_stack([event_times for i in range(len(sf))])
    kp_steps = list(zip(times.flatten(), sf_scores.flatten()))

    return uno_idx, ibs, mean_auc, dynamic_auc, brier_scores, kp_steps


