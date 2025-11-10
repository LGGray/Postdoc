# Import required libraries
import sys
import os.path
import pandas as pd
import numpy as np
import time
import pyreadr
from boruta import BorutaPy
from sklearn.model_selection import GridSearchCV, RepeatedKFold, GroupShuffleSplit
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score, cohen_kappa_score, matthews_corrcoef
from sklearn.feature_selection import RFECV
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_curve, auc, roc_auc_score
import matplotlib.pyplot as plt
from sklearn.utils import resample

start_time = time.time()

# Get the file name from the command line
file = sys.argv[1]
print(os.path.basename(file))

# Read in tune, train, test and features
X_train = pd.read_csv(f'pseudobulk/split_{sys.argv[2]}/data.splits/X_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_train = pd.read_csv(f'pseudobulk/split_{sys.argv[2]}/data.splits/y_train.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
X_test = pd.read_csv(f'pseudobulk/split_{sys.argv[2]}/data.splits/X_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
y_test = pd.read_csv(f'pseudobulk/split_{sys.argv[2]}/data.splits/y_test.'+os.path.basename(file).replace('.RDS', '')+'.csv', index_col=0)
enet_features = pd.read_csv(f'pseudobulk/split_{sys.argv[2]}/features/enet_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')
boruta_features = pd.read_csv(f'pseudobulk/split_{sys.argv[2]}/features/boruta_features.'+os.path.basename(file).replace('.RDS', '')+'.csv')

# Subset for selected and tentitive features from boruta
boruta_features = boruta_features[boruta_features['Rank'] == 1]
# Subset elastic net features to those with absolute value of coefficients in 80th percentile
threshold = np.percentile(np.abs(enet_features['coef']), 90)
enet_features = enet_features[np.abs(enet_features['coef']) >= threshold]

#### Condition for command-line argument indicating feature type ###
if sys.argv[3] == 'intersection':
    # Intersection of features selected by Boruta and Elastic Net
    features = pd.merge(enet_features, boruta_features, on='Feature', how='inner')['Feature']
    if(len(features) == 0):
        print("No common features between Boruta and Elastic Net")
        sys.exit()
elif sys.argv[3] == 'combined':
    # Features selected by Boruta and Elastic Net
    features = pd.merge(enet_features, boruta_features, on='Feature', how='outer')['Feature']
elif sys.argv[3] == 'boruta':
    features = boruta_features['Feature']
elif sys.argv[3] == 'enet':
    features = enet_features['Feature']

# Perform a grid search to find the best parameters
param_grid = {'learning_rate': [0.1, 0.05, 0.01],
              'n_estimators': [100, 200, 300],
              'max_features': ['sqrt', 0.3],
              'max_depth': [3, 5, 7]
}
clf = GradientBoostingClassifier(random_state=42, subsample=1)
grid_search = GridSearchCV(clf, param_grid, scoring='f1_weighted', 
                           cv=RepeatedKFold(n_splits=10, n_repeats=3, random_state=42), n_jobs=16, verbose=1)

# Fit the grid search object to the training data
grid_search.fit(X_train.loc[:,features], y_train['class'])

# Return estimator with best parameter combination
clf = grid_search.best_estimator_

# Predict the test set
y_pred = clf.predict(X_test.loc[:, features])
y_pred_proba = clf.predict_proba(X_test.loc[:, features])[:, 1]

# Calculate Youden's J statistic to find the optimal threshold
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
j_scores = tpr - fpr
# Find the optimal threshold
optimal_idx = np.argmax(j_scores)
optimal_threshold = thresholds[optimal_idx]
# Convert probabilities to binary predictions based on optimal threshold
y_pred = (y_pred_proba >= optimal_threshold).astype(int)

# Calculate the metrics
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
recall = recall_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred, average='weighted')
auroc = roc_auc_score(y_test, y_pred)
auprc = average_precision_score(y_test, y_pred)
kappa = cohen_kappa_score(y_test, y_pred)
mcc = matthews_corrcoef(y_test, y_pred)

# Define bootstrap parameters
n_bootstraps = 1000
confidence_level = 0.9
# Initialize an empty list to store bootstrap scores
bootstrapped_f1 = []
bootstrapped_AUC = []
bootstrapped_AUPRC = []
bootstrapped_mcc = []

# Loop over bootstrap samples
for i in range(n_bootstraps):
    # Resample with replacement
    y_test_resampled, y_pred_resampled = resample(y_test, y_pred, stratify=y_test)
    # Calculate F1 score
    f1 = f1_score(y_test_resampled, y_pred_resampled, average='weighted')
    # Calculate AUROC
    auroc = roc_auc_score(y_test_resampled, y_pred_resampled)
    # Calculate AUPRC
    auprc = average_precision_score(y_test_resampled, y_pred_resampled)
    # Calculate MCC
    mcc = matthews_corrcoef(y_test_resampled, y_pred_resampled)
    # Append score to list
    bootstrapped_f1.append(f1)
    bootstrapped_AUC.append(auroc)
    bootstrapped_AUPRC.append(auprc)
    bootstrapped_mcc.append(mcc)

# Calculate percentile for confidence intervals
lower_percentile = (1 - confidence_level) / 2 * 100
upper_percentile = (1 + confidence_level) / 2 * 100

f1_lower_bound = np.percentile(bootstrapped_f1, lower_percentile)
f1_upper_bound = np.percentile(bootstrapped_f1, upper_percentile)

auroc_lower_bound = np.percentile(bootstrapped_AUC, lower_percentile)
auroc_upper_bound = np.percentile(bootstrapped_AUC, upper_percentile)

auprc_lower_bound = np.percentile(bootstrapped_AUPRC, lower_percentile)
auprc_upper_bound = np.percentile(bootstrapped_AUPRC, upper_percentile)

mcc_lower_bound = np.percentile(bootstrapped_mcc, lower_percentile)
mcc_upper_bound = np.percentile(bootstrapped_mcc, upper_percentile)

# Create dataframe of metrics and save to file
metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1],
                        'F1_lower': [f1_lower_bound],
                        'F1_upper': [f1_upper_bound],
                        'AUC': [auroc],
                        'AUC_lower': [auroc_lower_bound],
                        'AUC_upper': [auroc_upper_bound],
                        'AUPRC': [auprc],
                        'AUPRC_lower': [auprc_lower_bound],
                        'AUPRC_upper': [auprc_upper_bound],
                        'Kappa': [kappa],
                        'MCC': [mcc],
                        'MCC_lower': [mcc_lower_bound],
                        'MCC_upper': [mcc_upper_bound],
                        'n_features': [len(features)]})
metrics.to_csv(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/metrics/GBM_metrics_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Save confusion matrix to file
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))
confusion.to_csv(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/metrics/GBM_confusion_'+os.path.basename(file).replace('.RDS', '')+'.csv', index=False)

# Define class names
classes = ['Control', 'Disease']
fig, ax = plt.subplots()
# Set the color map to 'coolwarm'
cmap = plt.cm.coolwarm
# Create the heatmap for the confusion matrix
cax = ax.matshow(confusion, cmap=cmap)
# Add color bar
plt.colorbar(cax)
# Add counts to the confusion matrix cells
confusion_values = confusion.values
for (i, j), val in np.ndenumerate(confusion_values):
    ax.text(j, i, f'{val}', ha='center', va='center', color='black')
# Set axis labels
ax.set_xlabel('Predicted labels')
ax.set_ylabel('True labels')
ax.set_xticks(range(len(classes)))
ax.set_yticks(range(len(classes)))
ax.set_xticklabels(classes)
ax.set_yticklabels(classes)
# Set the title
ax.set_title('GBM: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
# Annotate with F1 score
plt.annotate(f'MCC: {mcc:.2f}', xy=(0.5, -0.1), xycoords='axes fraction', 
             ha='center', va='center', fontsize=12, color='black')
# Adjust layout for visibility
plt.tight_layout()
# Save the figure
plt.savefig(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/confusion/GBM_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')
plt.close()

# Print the AUROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
plt.figure()
plt.plot(fpr, tpr, label='AUC-ROC (area = %0.2f)' % auroc)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('GBM: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.legend(loc="lower right")
plt.savefig(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/AUROC/GBM_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title('GBM: ' + os.path.basename(file).replace('.RDS', '').replace('.', ' '))
plt.savefig(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/PRC/GBM_'+os.path.basename(file).replace('.RDS', '')+'.pdf', bbox_inches='tight')

# Save the model
import pickle
filename = f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/GBM_model_'+os.path.basename(file).replace('.RDS', '')+'.sav'
pickle.dump(clf, open(filename, 'wb'))

end_time = time.time()
cpu_time = end_time - start_time

print(f"CPU time used: {cpu_time:.2f} seconds")