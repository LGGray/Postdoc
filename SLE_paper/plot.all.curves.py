import sys
import os
import pickle
import pandas as pd
import numpy as np
import pyreadr
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc, roc_auc_score, precision_recall_curve, PrecisionRecallDisplay, average_precision_score
import matplotlib.pyplot as plt

file = sys.argv[1]
cell = file.replace('pseudobulk/', '').replace('.RDS', '')

# load the model from disk
logit = pickle.load(open(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/logit_model_'+cell+'.sav', 'rb'))
RF = pickle.load(open(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/RF_model_'+cell+'.sav', 'rb'))
SVM = pickle.load(open(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/SVM_model_'+cell+'.sav', 'rb'))
GBM = pickle.load(open(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/GBM_model_'+cell+'.sav', 'rb'))
MLP = pickle.load(open(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/MLP_model_'+cell+'.sav', 'rb'))

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

# Get the predicted probabilities
logit_pred_proba = logit.predict_proba(X_test.loc[:, features])[:, 1]
RF_pred_proba = RF.predict_proba(X_test.loc[:, features])[:, 1]
SVM_pred_proba = SVM.predict_proba(X_test.loc[:, features])[:, 1]
GBM_pred_proba = GBM.predict_proba(X_test.loc[:, features])[:, 1]
MLP_pred_proba = MLP.predict_proba(X_test.loc[:, features])[:, 1]

# Create a dictionary of model names and predicted probabilities
models = {'logit': logit_pred_proba,
          'RF': RF_pred_proba,
          'SVM': SVM_pred_proba,
          'GBM': GBM_pred_proba,
          'MLP': MLP_pred_proba}
# Loop through the dictionary and plot models
for name, proba in models.items():
    # Calculate the precision and recall for each model
    precision, recall, _ = precision_recall_curve(y_test['class'], proba)
    # Plot the curve and add the label
    plt.plot(recall, precision, label=name)
# Add labels and legend
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve: ' + cell.replace('.', ' '))
plt.legend()
plt.savefig(f'pseudobulk/split_{sys.argv[2]}/{sys.argv[3]}/ML.models/PRCurve_'+ cell +'.pdf', dpi=300)
plt.show()
plt.close()