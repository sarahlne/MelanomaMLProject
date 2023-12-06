## MelanoMLProject

Machine Learning models using patient related clinical and genomic data for the prediction of response to MAPKi treatments. 

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