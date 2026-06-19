import pandas as pd
import numpy as np
import imbens
import sklearn
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold


# Input / output folders
input_fold = '../result/4_DESWAN_for_gene_expression/sig_df'
output_fold = '../result/5_peak_gene_model'


# 1. First-pass feature selection: select genes with adjusted p < 1e-4 at age 6
df = pd.read_csv(f"{input_fold}/peak_sig_gene.csv", index_col=0)
df_6 = df[df['age'] == 6]
df_6[df_6['q_value'] < 0.0001].to_csv(f"{output_fold}/peak_sig_gene_X6.csv")


# 2. Load expression and clinical data for training, test, and disease cohorts

##2.1 train and healthy validation cohort
train_df = pd.read_csv(f"{input_fold}/train/train_expr_peak6_0001.csv", index_col=0)
test_df = pd.read_csv(f"{input_fold}/pj/pj_expr_peak6_0001.csv", index_col=0)
cli = pd.read_csv(f"{input_fold}/cli.csv", index_col=0)
cli_pj = pd.read_csv(f"{input_fold}/pj/pj_cli.csv", index_col=0)

## 2.2 T1D
t1d_df = pd.read_csv(f"{input_fold}/T1D/T1D_expr_peak6_0001.csv", index_col=0)
cli_t1d = pd.read_csv(f"{input_fold}/T1D/T1D_cli.csv", index_col=0)

## 2.3 SLE
sle_df = pd.read_csv(f"{input_fold}/SLE/SLE_expr_peak6_0001.csv", index_col=0)
cli_sle = pd.read_csv(f"{input_fold}/SLE/SLE_cli.csv", index_col=0)

## 2.4 CD
cd_df = pd.read_csv(f"{input_fold}/CD/CD_expr_peak6_0001.csv", index_col=0)
cli_cd = pd.read_csv(f"{input_fold}/CD/CD_cli.csv", index_col=0)

# 3. Data preprocessing

train_df.dropna(axis=1, inplace=True)
test_df.dropna(axis=1, inplace=True)

cli['age_group'] = (cli['age'] > 6).astype(int)
cli_pj = cli_pj.loc[test_df.columns]
cli_pj['age_group'] = (cli_pj['age'] > 6).astype(int)

y_train = cli['age_group']
X_train = train_df.T

y_test = cli_pj['age_group']
X_test = test_df.T


# Ensure all cohorts contain the same features (genes)
selected = list(
    set(X_train.columns)
    & set(X_test.columns)
    & set(sle_df.index)
    & set(t1d_df.index)
    & set(cd_df.index)
)
X_train = X_train[list(selected)]


# 4. Second-pass feature selection: RFECV (recursive feature elimination with CV)
strkfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
base_logistic = sklearn.linear_model.LogisticRegression(random_state=42, C=0.1, max_iter=1000)
rfecv = RFECV(estimator=base_logistic, cv=strkfold, n_jobs=-1, scoring='average_precision')
rfecv.fit(X_train, y_train)

selected_features = X_train.columns[rfecv.support_]

X_train = X_train[list(selected_features)]
X_test = X_test[list(selected_features)]


# 5. Model training using SelfPacedEnsembleClassifier from imbens
np.random.seed(42)
base_logistic = sklearn.linear_model.LogisticRegression(random_state=42, C=0.1, max_iter=1000)
logistic = imbens.ensemble.SelfPacedEnsembleClassifier(n_estimators=200, estimator=base_logistic)
logistic.fit(X_train, y_train)


import pickle
with open(f'{output_fold}/model.pkl', 'wb') as f:
    pickle.dump(logistic, f)


# 6. Score the training cohort and the validation cohort
## 6.1 Score the training cohort
score_train = logistic.predict_proba(X_train)[:, 1]
result_train = pd.DataFrame({'age': cli['age'], 'remodeling_score': score_train, 'sex': cli['sex']})
result_train.to_csv(f'{output_fold}/train_result.csv', index=False)

## 6.2 Score the validation cohorts
def IDS(df, cli_df, name):
    # df: expression matrix with genes as index and samples as columns
    # cli_df: clinical table indexed by sample IDs
    df = df.loc[selected_features, :]
    df.dropna(axis=1, inplace=True)

    cli_df.index = cli_df.index.astype(str)
    cli_sub = cli_df.loc[df.columns]
    cli_sub['age_group'] = (cli_sub['age'] > 6).astype(int)

    # Predict probabilities
    X = df.T
    y_prob = logistic.predict_proba(X)[:, 1]

    result = pd.DataFrame({
        'remodeling_score': y_prob,
        'age': cli_sub['age'],
        'sex': cli_sub['sex']
    })
    result.to_csv(f"{output_fold}/{name}_result.csv", index=False)
    return result


disease_data = {
    'test': (test_df, cli_pj),  # healthy validation cohort
    'SLE': (sle_df, cli_sle),
    'CD': (cd_df, cli_cd),
    'T1D': (t1d_df, cli_t1d)
}

results = {}
for name, (df_cohort, cli_cohort) in disease_data.items():
    results[name] = IDS(df_cohort, cli_cohort, name)


# 7. Permutation importance calculation for selected features
from sklearn.inspection import permutation_importance

perm_res = permutation_importance(
    logistic, X_test, y_test, n_repeats=100, random_state=42, scoring='average_precision'
)

imp = perm_res.importances_mean
rel_imp = np.abs(imp) / np.max(np.abs(imp))

df_imp = pd.DataFrame({
    'feature': list(selected_features),
    'importance_mean': perm_res.importances_mean,
    'relative_importance': rel_imp
}).sort_values('relative_importance', ascending=False)

df_imp[['ct', 'gene']] = df_imp['feature'].str.rsplit('_', n=1, expand=True)
df_imp.to_csv(f"{output_fold}/model_importance.csv", index=False)