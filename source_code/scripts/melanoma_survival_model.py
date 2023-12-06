""" Module with class method for running the survival models """



import pandas as pd
import preprocessing_data as prep



class SurvivalModel:
    """ Method for running the melanoma survival model """

    def __init__(
        self, model, hyperparameters, datasets, model_name=None,
        save_models=False, results_folder='../../source_code/model_results/',
        n_trials=10, n_seeds=[42, 2678, 4116, 2479, 824, 3439, 7958, 4676, 761, 1827], ):

        self.model = model
        self.params = hyperparameters
        self.all_datasets = datasets
        self.n_seeds = n_seeds
        self.n_trials = n_trials
        self.model_name = str(model)[:str(model).find('(')] if model_name is None else model_name
        self.results_folder = results_folder
        self.save_models = save_models

        self.ncv_models = None
    
    def perform_nested_cv(self,):
        """performs nested cross validatation procedure"""

        ncv_models = list()

        for trail, seed in enumerate(self.n_seeds[ : self.n_trials], 1):
            print('\n','='*10, f'TRIAL {trail}','='*10)
            MODEL = self.model
            if 'random_state' in MODEL.get_params():
                    MODEL.set_params(random_state=seed) # for reproducibility
            for name, dataset in self.all_datasets.items():
                # seperate X and y
                X, y = prep.seperate_x_y(dataset)
                # perform nested cv
                of_models = prep.perform_nested_cv(X, y, MODEL, self.params, outer_splits=5, inner_splits=3, seed=seed)
                # store best models
                for model in of_models:
                    model_info = pd.Series([self.model_name, trail, name], index=['model', 'trail', 'dataset'])
                    of_model = pd.Series(model, index= ['outer_fold', 'best_model', 'X_train', 'y_train', 'X_test', 'y_test'])
                    ncv_model = pd.concat([model_info, of_model])
                    ncv_models.append(ncv_model)
                    # print('fold',model[0])
                #print(name)
        print('\n','='*10, 'END OF TRIALs','='*10)

        # all nested cv models
        self.ncv_models = pd.DataFrame(ncv_models)
        return None
    
    def save_ncv_models(self, path=None, ):
        """Stores outer CV models"""

        if self.ncv_models is None:
            raise ValueError('No model executed. Perform nested CV.')
            
        if path is None:
            path = self.results_folder+f'{self.model_name}_ncv_models.pkl'
        else:
            path = path+f'{self.model_name}_ncv_models.pkl'
        # store nested cv models
        self.ncv_models.to_pickle(path)
        print(f'{self.model_name} models saved at location: {path}')
        return None
    
    def perform_evaluation(self, path=None, ncv_models=None):
        """performs evaluation of the outer CV models"""

        if ncv_models is None and self.ncv_models is None:
            raise ValueError('No models available. Perform nested CV.')
        elif ncv_models is not None:
            NCV_MODELS = ncv_models
        else:
            NCV_MODELS = self.ncv_models

        from sksurv.metrics import concordance_index_ipcw as cindex_u
        from sksurv.metrics import cumulative_dynamic_auc, integrated_brier_score, brier_score
        from math import ceil, floor
        import numpy as np
        #import pandas as pd

        eval_models = []
        for i, row in NCV_MODELS.iterrows():
            if str(self.model)[:str(self.model).find('(')] == 'FastSurvivalSVM':
                # Specific Evaluation for Support Vector Machines
                model_info = pd.Series([row.model, row.trail, row.dataset, row.outer_fold], index=['model', 'trail', 'dataset', 'outer_fold'])

                ## Uno'c-index
                yhat_test = row.best_model.predict(row.X_test)
                uno_idx = cindex_u(row.y_train,row.y_test, yhat_test)[0] # Uno's c-index
                ## time-dependent scores
                #t_min, t_max = prep.time_intervals(row.y_train, row.y_test)
                #event_times = np.arange(ceil(t_min), floor(t_max)+1, step=1, dtype=int)
                ## dynamic AUC
                #td_auc = cumulative_dynamic_auc(row.y_train, row.y_test, yhat_test, event_times)[0]
                #dynamic_auc = list(zip(event_times, td_auc))
                #mean_auc = np.nanmean(td_auc)

                #scores = [uno_idx, np.nan, mean_auc, dynamic_auc, np.nan, np.nan]
                scores = [uno_idx, np.nan, np.nan, np.nan, np.nan, np.nan]
                results = pd.Series(scores, index=['c_index', 'ibs', 'mean_auc', 'dynamic_auc', 'brier_scores', 'kaplan_meier'])
                evaluated_model = pd.concat([model_info, results])
                eval_models.append(evaluated_model)

            else:
                try:
                    model_info = pd.Series([row.model, row.trail, row.dataset, row.outer_fold], index=['model', 'trail', 'dataset', 'outer_fold'])
                    scores = prep.evaluation(row.best_model, row.X_train, row.y_train, row.X_test, row.y_test)
                    results = pd.Series(scores, index=['c_index', 'ibs', 'mean_auc', 'dynamic_auc', 'brier_scores', 'kaplan_meier'])
                    evaluated_model = pd.concat([model_info, results])
                    eval_models.append(evaluated_model)

                except ValueError as ve:
                    print(f'Evaluation error index: {i} ', '_'.join([row.model, str(row.trail), row.dataset, str(row.outer_fold)]))
                    model_info = pd.Series([row.model, row.trail, row.dataset, row.outer_fold], index=['model', 'trail', 'dataset', 'outer_fold'])
                    event_times = np.arange(2, 20, step=1, dtype=int)
                    ## Uno'c-index
                    yhat_test = row.best_model.predict(row.X_test)
                    uno_idx = cindex_u(row.y_train,row.y_test, yhat_test)[0] # Uno's c-index
                    ## cumulative hazard function
                    #chf = row.best_model.predict_cumulative_hazard_function(row.X_test)
                    #chf_scores = np.row_stack([func(event_times) for func in chf])
                    ## dynamic AUC
                    #td_auc = cumulative_dynamic_auc(row.y_train, row.y_test, chf_scores, event_times)[0]
                    #dynamic_auc = list(zip(event_times, td_auc))
                    #mean_auc = np.nanmean(td_auc)

                    #scores = [uno_idx, np.nan, mean_auc, dynamic_auc, np.nan, np.nan]
                    scores = [uno_idx, np.nan, np.nan, np.nan, np.nan, np.nan]
                    results = pd.Series(scores, index=['c_index', 'ibs', 'mean_auc', 'dynamic_auc', 'brier_scores', 'kaplan_meier'])
                    evaluated_model = pd.concat([model_info, results])
                    eval_models.append(evaluated_model)
        
        # evaluated model results
        self.eval_scores = pd.DataFrame(eval_models)
        # save results
        if path is None:
            path = self.results_folder+f'{self.model_name}_ncv_results.pkl'
        else:
            path = path+f'{self.model_name}_ncv_results.pkl'
        self.eval_scores.to_pickle(path)
        print(f'{self.model_name} results saved at location: {path}')
        return None

    def run_survival_model(self, evaluation=True):
        """Runs the model with nested cross validatation procedure, performs evaluation 
        and stores the results"""

        self.perform_nested_cv()
        if self.save_models == True:
            self.save_ncv_models()
        if evaluation == True:
            self.perform_evaluation()

        return None

