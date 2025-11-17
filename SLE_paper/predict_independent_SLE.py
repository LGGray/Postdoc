import sys
import pickle
import os.path
import pandas as pd
import numpy as np
import pyreadr
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score, cohen_kappa_score, matthews_corrcoef
import matplotlib.pyplot as plt
from sklearn.ensemble import VotingClassifier
from sklearn.inspection import permutation_importance
from sklearn.utils import resample

cell = os.path.basename(sys.argv[1]).replace('.RDS', '')

# Load in ensemble model
model_path = f'/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/SLE_paper/SLE/pseudobulk_update/split_{sys.argv[2]}/combined/ensemble/{cell}.sav'
eclf = pickle.load(open(model_path, 'rb'))
features = eclf.feature_names_in_

### Read in independent test set
test = pyreadr.read_r(sys.argv[1])
test = test[None]

# Subset for adult, child or all
if sys.argv[3] == 'adult':
    # Subset class to include aHD and aSLE
    test = test[test['class'].isin(['aHD', 'aSLE'])]
if sys.argv[3] == 'child':
    # Subset class to include cHD and cSLE
    test = test[test['class'].isin(['cHD', 'cSLE'])]
if sys.argv[3] == 'all':
    # Subset class to include all
    test = test[test['class'].isin(['cHD', 'cSLE', 'aHD', 'aSLE'])]

# Replace class labels with 0 and 1
test['class'] = test['class'].replace({"cSLE": 1, "aSLE": 1, "cHD": 0, "aHD": 0})

#test['class'] = test['class'].replace({"disease": 1, "control": 0})

# Set class as target variable
y_test = test[['class']]

# Check which features are missing and fill with 0
missing = np.setdiff1d(features, test.columns)
test[missing] = 0
X_test = test[features]

# Scale the data
scaler = StandardScaler()
X_test = pd.DataFrame(scaler.fit_transform(X_test), columns=X_test.columns, index=X_test.index)

# Predict the test set
y_pred = eclf.predict(X_test)
y_pred_proba = eclf.predict_proba(X_test)[:, 1]

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

# Create dataframe of metrics and save to file
metrics = pd.DataFrame({'Accuracy': [accuracy], 
                        'Precision': [precision], 
                        'Recall': [recall], 
                        'F1': [f1],
                        'AUC': [auroc],
                        'AUPRC': [auprc],
                        'Kappa': [kappa],
                        'MCC': [mcc],
                        'n_features': [len(features)]})
metrics.to_csv(f'/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/SLE_paper/SLE_GSE135779/results_update/metrics_{cell}_{sys.argv[3]}_split_{sys.argv[2]}.csv', index=False)

# confusion matrix
confusion = pd.DataFrame(confusion_matrix(y_test, y_pred))

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
ax.set_title(f'Split_{sys.argv[1]}: {sys.argv[3]} {cell.replace("_", " ")}')
# Annotate with F1 score
plt.annotate(f'MCC: {mcc:.2f}', xy=(0.5, -0.1), xycoords='axes fraction', 
             ha='center', va='center', fontsize=12, color='black')
# Adjust layout for visibility
plt.tight_layout()
# Save the figure
plt.savefig(f'/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/SLE_paper/SLE_GSE135779/results_update/confusion_{cell}_{sys.argv[3]}_split_{sys.argv[2]}.pdf', bbox_inches='tight')
plt.close()

# Print the PR curve
precision, recall, thresholds = precision_recall_curve(y_test, y_pred_proba)
average_precision = average_precision_score(y_test, y_pred_proba)
disp = PrecisionRecallDisplay(precision=precision, recall=recall, average_precision=average_precision)
disp.plot()
disp.ax_.set_title(f'Split_{sys.argv[2]}: {sys.argv[3]} {cell.replace("_", " ")}')
plt.savefig(f'/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/SLE_paper/SLE_GSE135779/results_update/PRcurve_{cell}_{sys.argv[3]}_split_{sys.argv[2]}.pdf', bbox_inches='tight')