from sklearn.model_selection import train_test_split
import random
from sklearn.linear_model import SGDRegressor
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score,explained_variance_score
from scipy.stats import pearsonr
from sklearn.linear_model import BayesianRidge
from sklearn.model_selection import KFold
from scipy.stats import spearmanr




def get_BayesianRidge_prediction(x_train, y_train, x_val, y_val, alpha_1, alpha_2, lambda_1, lambda_2):
    model = BayesianRidge(alpha_1=alpha_1, alpha_2=alpha_2, lambda_1=lambda_1, lambda_2=lambda_2)
    model.fit(x_train, y_train)
    y_pred = model.predict(x_val)
    return model,pearsonr(y_val, y_pred)[0]



def full_fit_BayesianRidge(x_train, x_test, y_train, y_test, alpha_1, alpha_2, lambda_1, lambda_2):
    # Bayesian Ridge Regression Method
    #
    # Note prior Gamma distribution are set as ( alpha_1, alpha_2) and (lambda_1, lambda_2) which were selected via cross-validation step using training data (see codes below: full_fit_BayesianRidge_para_turing;
    # all traits shared the same)
    # X_train: training genotype data  (Numpy matrix - samples X viriants)
    # X_test: testing genotype data
    # y_train: training trait value data (Numpy vector - 1 X N)
    # y_test: testing trait value data
    # return the learned model, r, explained variance score and spearmanr performance
    model, r = get_BayesianRidge_prediction(x_train, y_train, x_test, y_test, alpha_1, alpha_2, lambda_1, lambda_2)
    y_pred = model.predict(x_test)
    return model,pearsonr(y_test, y_pred)[0],r2_score(y_test, y_pred),explained_variance_score(y_test,y_pred),spearmanr(y_test,y_pred)[0]



#hyper-parameter Tuning - finding the best prior Gamma distributions on the training set of a trait
#Grid search on (1e10, 1e5, 1e1, 0, -1e1, -1e5, -1e10)
#return the best 'alpha_1', 'alpha_2' 'lambda_1' 'lambda_2'
def full_fit_BayesianRidge_para_turing(x_train, x_val, y_train, y_val,para_file):
    nums = (1e10, 1e5, 1e3, 1e1, 0, -1e1, -1e3,-1e5,-1e10)
    f=open(para_file,'w')
    alpha_1 = nums
    alpha_2 = nums
    lambda_1 = nums
    lambda_2 = nums
    best_model = None
    best_r = 0
    for a1 in alpha_1:
        for a2 in alpha_2:
            for l1 in lambda_1:
                for l2 in lambda_2:
                    model,r = get_BayesianRidge_prediction(x_train,y_train,x_val,y_val,a1,a2,l1,l2)
                    text = "Training BayesianRidge with alpha_1: {}, alphs_2: {}, lambda_1: {}, lambda_2:{} - r score {}\n".format(a1,a2,l1,l2,r)
                    print(text)
                    f.write(text)
                    if best_model == None or r > best_r:
                        best_model = model
                        best_r = r
                        best_params =  {'alpha_1':a1, 'alpha_2': a2, 'lambda_1': l1, 'lambda_2': l2}
    print("Best Para: {}".format(best_params))
    return best_model




