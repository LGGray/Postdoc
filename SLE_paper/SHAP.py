import pickle
import pandas as pd
import sys
import os
import numpy as np
from collections import Counter
from sklearn.ensemble import VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import matthews_corrcoef
import matplotlib.pyplot as plt
import shap

file = sys.argv[1]

cell = file.replace('pseudobulk/', '').replace('.RDS', '')

# Read in ensembe models
ensemble_models = []
for i in range(1, 11):
    model = pickle.load(open(f'pseudobulk_update/split_{i}/combined/ensemble/{cell}.sav', 'rb'))
    ensemble_models.append(model)

# Function to calculate the average or mode of hyperparameters
def average_hyperparams(models, param_path):
    values = [eval(f'model.{param_path}') for model in models]
    if isinstance(values[0], (int, float)):
        return int(np.mean(values))
    else:
        return Counter(values).most_common(1)[0][0]
    
# Averaging or taking the mode of hyperparameters
if average_hyperparams(ensemble_models, 'estimators[0][1].C') == 0:
    average_logit_C = 0.1
else:
    average_logit_C = average_hyperparams(ensemble_models, 'estimators[0][1].C')

average_logit_penalty = average_hyperparams(ensemble_models, 'estimators[0][1].penalty')
average_logit_solver = average_hyperparams(ensemble_models, 'estimators[0][1].solver')
average_logit_l1_ratio = average_hyperparams(ensemble_models, 'estimators[0][1].l1_ratio')

average_rf_max_depth = average_hyperparams(ensemble_models, 'estimators[1][1].max_depth')
average_rf_min_samples_split = average_hyperparams(ensemble_models, 'estimators[1][1].min_samples_split')
average_rf_n_estimators = average_hyperparams(ensemble_models, 'estimators[1][1].n_estimators')

if average_hyperparams(ensemble_models, 'estimators[2][1].C') == 0:
    average_svm_C = 0.1
else:
    average_svm_C = average_hyperparams(ensemble_models, 'estimators[2][1].C')

#average_gbm_max_features = average_hyperparams(ensemble_models, 'estimators[3][1].max_features')

average_mlp_alpha = average_hyperparams(ensemble_models, 'estimators[4][1].alpha')
average_mlp_max_iter = average_hyperparams(ensemble_models, 'estimators[4][1].max_iter')

# Creating the new model with the average/mode hyperparameters
ensemble = VotingClassifier(
    estimators=[
        ('logit', LogisticRegression(C=0.1, penalty=average_logit_penalty, solver=average_logit_solver, l1_ratio=average_logit_l1_ratio, class_weight='balanced', max_iter=10000, n_jobs=8, random_state=42)),
        ('RF', RandomForestClassifier(max_depth=average_rf_max_depth, min_samples_split=average_rf_min_samples_split, n_estimators=average_rf_n_estimators, class_weight='balanced', n_jobs=8)),
        ('SVM', SVC(C=0.1, class_weight='balanced', probability=True, random_state=42)),
        ('GBM', GradientBoostingClassifier(max_features='sqrt', random_state=42, subsample=1)),
        ('MLP', MLPClassifier(alpha=average_mlp_alpha, early_stopping=True, max_iter=average_mlp_max_iter, random_state=42))
    ],
    voting='soft'
)

all_features = pd.read_csv('results_update/all_features.csv')
features = all_features[all_features['celltype'] == cell]['feature'].tolist()

X_train = pd.read_csv(f'pseudobulk_update/split_1/data.splits/X_train.'+cell+'.csv', index_col=0)
y_train = pd.read_csv(f'pseudobulk_update/split_1/data.splits/y_train.'+cell+'.csv', index_col=0)
X_test = pd.read_csv(f'pseudobulk_update/split_1/data.splits/X_test.'+cell+'.csv', index_col=0)
y_test = pd.read_csv(f'pseudobulk_update/split_1/data.splits/y_test.'+cell+'.csv', index_col=0)

# Fit the ensemble model
ensemble.fit(X_train.loc[:, features], y_train['class'])

# predict the test set
y_pred = ensemble.predict(X_test.loc[:, features])

# calculate the MCC
mcc = matthews_corrcoef(y_test['class'], y_pred)

# Create a wrapper function for the predict_proba method of VotingClassifier
def voting_classifier_proba(data):
    return ensemble.predict_proba(data)

explainer = shap.KernelExplainer(voting_classifier_proba, X_train.loc[:, features])

# shap_values = explainer.shap_values(X_test.loc[:, features])

explanation = explainer(X_test.loc[:, features])

shap_values_single_class = explanation[..., 1]  # Adjust index based on the class you are interested in
shap.plots.beeswarm(shap_values_single_class, max_display=len(features))
plt.title(cell)
plt.savefig(f'results_update/SHAP/{cell}_shap.beeswarm.pdf', bbox_inches='tight')
plt.close()