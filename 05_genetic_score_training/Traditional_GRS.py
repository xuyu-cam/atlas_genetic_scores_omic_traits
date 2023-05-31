
import numpy as np
from sklearn.metrics import r2_score,explained_variance_score
from scipy.stats import pearsonr
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr


def get_beta_vec_by_vars_ids(beta_file,vars):
    df = pd.read_csv(beta_file,sep='\t',index_col=0)
    betas = df.loc[vars]['effect']
    return np.array(betas)


def traditional_GRS_selected_vars(beta_file,X,y,vars):
    beta_vec = get_beta_vec_by_vars_ids(beta_file,vars)
    y_pred = X.dot(beta_vec)
    y_pred = stats.zscore(y_pred)
    return pearsonr(y, y_pred)[0],r2_score(y, y_pred),explained_variance_score(y,y_pred),spearmanr(y, y_pred)[0]

