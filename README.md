## MelanoMLProject

Source code for "Improved prediction of the response duration to MAPK inhibitors in patients with advanced melanoma using baseline genomic data and machine learning algorithms" (Accepted in NPJ Precision Oncology on 25.01.2025) (https://www.medrxiv.org/content/10.1101/2023.12.07.23299389v1.full.pdf+html)


## Abstract
Baseline genomic data have not demonstrated significant value for predicting the response duration to MAPK inhibitors (MAPKi) in patients with BRAFV600-mutated melanoma. We used machine learning algorithms and pre-processed genomic data to test whether they could contain useful information to improve the progression-free survival (PFS) prediction. This exploratory analysis compared the predictive performance of a dataset that contained clinical features alone and supplemented with baseline genomic data. Whole and partial exon sequencing data from four cohorts of patients with BRAFV600-mutated melanoma treated with MAPKi were used: two cohorts as training/evaluation set (n = 111) and two as validation set (n = 73). Genomic data were pre-processed using three strategies to generate eight different genomic datasets. Several machine learning algorithms and one statistical algorithm were employed to predict PFS. The performance of these survival models was assessed using the concordance index, time-dependent receiver operating characteristic (ROC) curve and Brier score. The cross-validated model performance improved when pre-processed genomic data, such as mutation rates, were added to the clinical features. In the validation dataset, the best model with genomic data outperformed the best model with clinical features alone. The trend towards improved prediction with baseline genomic data was maintained when data were censored according to the two clinical setting scenarios (duration of clinical benefit and progression before 12 months). Finally, our best model outperformed with baseline genomic data, increasing the number of patients with a correctly predicted relapse by between +12% and +28%. In our models, baseline genomic data improved the prediction of response duration and could be incorporated into the development of predictive models of MAPKi treatment in melanoma. 

----------------------

## Authors

- [@sarahlne](https://www.github.com/sarahlne)
- [@kritiAT](https://github.com/kritiAT)

## Source code
Contain all necessary files to run survival models from melanodb data (melanodb_original.db & melanodb_validation.db)

- **scripts**
    - load data from database
    - generate datasets
    - define survival machine learning (ML) models
    - process data for survival ML algorithms
    - run models and store results in /model_reuslts folder
    - 
- **datafiles**
Datafiles required to generated datasets

- **model_results**
Results of ML algorithms stored in pickle files

- **notebooks**
  Visualization of results and figures generation

## Usage/Examples

```python
py ./source_code/run_all_survival_models.py
```
