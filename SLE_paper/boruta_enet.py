"""Boruta and Elastic Net feature selection pipeline.

Performs Boruta (Random Forest-based) and Elastic Net feature selection
on pseudobulk gene expression data for SLE classification.

Usage:
    python boruta_enet.py <input_RDS_file> <split_number>

Arguments:
    input_RDS_file: Path to an RDS file containing pseudobulk expression matrix.
    split_number:   Integer (1-10) indicating the cross-validation split.
"""

import sys
import os.path
import time
import pandas as pd
import numpy as np
from numpy import arange
import pyreadr
from boruta import BorutaPy
from sklearn.utils import resample
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, RepeatedKFold, GroupShuffleSplit
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import ElasticNetCV
import pickle

start_time = time.time()

# Get the file name from the command line
file = sys.argv[1]
print(os.path.basename(file))
cell = os.path.basename(file).replace('.RDS', '')

# Read in expression RDS file
df = pyreadr.read_r(file)
df = df[None]
print(df.head())

# Replace class with binary label
df['class'] = df['class'].replace({"control": 0, "disease": 1})

#### REMOVE ANCESTRY, AGE, and CELLCOUNT as features ####
df = df.drop(['ancestry', 'age', 'cellCount'], axis=1)

# Read in data splits
train = pd.read_csv(f'pseudobulk/split_{sys.argv[2]}/train_index.csv')
test = pd.read_csv(f'pseudobulk/split_{sys.argv[2]}/test_index.csv')

# Get the training and testing data
X_train = df[df['individual'].isin(train['rownames'])].drop(['class', 'individual'], axis=1)
y_train = df[df['individual'].isin(train['rownames'])]['class']
X_test = df[df['individual'].isin(test['rownames'])].drop(['class', 'individual'], axis=1)
y_test = df[df['individual'].isin(test['rownames'])]['class']

# Standard scale the data - z-scores
scaler = StandardScaler()
# Fit ONLY on training data
X_train = pd.DataFrame(
    scaler.fit_transform(X_train),
    columns=X_train.columns,
    index=X_train.index
)

with open(f'pseudobulk_update/split_{sys.argv[2]}/scaler/scaler_'+cell+'.pkl', "wb") as f:
    pickle.dump(scaler, f)

# Use the SAME scaler for test
X_test = pd.DataFrame(
    scaler.transform(X_test),
    columns=X_test.columns,
    index=X_test.index
)

# Save data splits
X_train.to_csv(f'pseudobulk_update/split_{sys.argv[2]}/data.splits/X_train.'+cell+'.csv', index=True)
y_train.to_csv(f'pseudobulk_update/split_{sys.argv[2]}/data.splits/y_train.'+cell+'.csv', index=True)
X_test.to_csv(f'pseudobulk_update/split_{sys.argv[2]}/data.splits/X_test.'+cell+'.csv', index=True)
y_test.to_csv(f'pseudobulk_update/split_{sys.argv[2]}/data.splits/y_test.'+cell+'.csv', index=True)

### Boruta feature selection ###
X = X_train.values
y = y_train.values.ravel()

# random forest classifier
param_grid = {'n_estimators': [100, 200, 300],
              'max_features': ['sqrt', 0.3],
                'max_depth': [5, 10, 15],
                'min_samples_split': [2, 5, 8]
}
clf = RandomForestClassifier(n_jobs=16, class_weight='balanced', criterion='gini')
grid_search = GridSearchCV(clf, param_grid, scoring='f1_weighted',
                           cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=42), n_jobs=16, verbose=1)

# Fit the grid search object to the training data
grid_search.fit(X, y)
# Return estimator with best parameter combination
rf = grid_search.best_estimator_
# define Boruta feature selection method
feat_selector = BorutaPy(rf, n_estimators='auto', verbose=2, random_state=42)
# find all relevant features
feat_selector.fit(X, y)

# Save the model
filename = f'pseudobulk_update/split_{sys.argv[2]}/feature.select.model/boruta_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(feat_selector, open(filename, 'wb'))

# Get feature rankings
feature_ranks = list(feat_selector.ranking_)
# Create a DataFrame with features, their importance, and ranks
feature_df = pd.DataFrame({
    'Feature': X_train.columns,
    'Rank': feature_ranks
})
# Sort the DataFrame based on feature ranks
feature_df.sort_values(by='Rank', ascending=True, inplace=True)
# Save the feature importance to file
feature_df.to_csv(f'pseudobulk_update/split_{sys.argv[2]}/features/boruta_features.'+cell+'.csv', index=False)

### Elastic net feature selection ###
ratios = [.1, .5, .7, .9, .95, .99, 1]
alphas = np.logspace(-4, 1, 50)
cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=42)
enet = ElasticNetCV(l1_ratio=ratios, alphas=alphas, max_iter=50000, cv=cv, n_jobs=16, random_state=42)
enet.fit(X_train, y_train.ravel())
print(enet)

# Create a dataframe of the features and their coefficients
enet_features = pd.DataFrame({'Feature': enet.feature_names_in_, 'coef': enet.coef_})
# Save the features to file
enet_features.to_csv(f'pseudobulk_update/split_{sys.argv[2]}/features/enet_features.'+cell+'.csv', index=False)

# Save the model
filename = f'pseudobulk_update/split_{sys.argv[2]}/feature.select.model/enet_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(enet, open(filename, 'wb'))

end_time = time.time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")