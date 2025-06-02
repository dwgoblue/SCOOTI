"""
regressorCollection.py
====================================================================
A collection of regression models for meta-learner regression models
"""

# Essential packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
#sys.path.append('../GeneralMethods/')
#from AnalysisKits import *
from SCOOTI.MatplotProp import CanvasStyle, PltProps
from SCOOTI.regressionAnalyzer import *
import networkx as nx
import pingouin as pg
PltProps()
import warnings; warnings.simplefilter('ignore')
from statsmodels.stats.multitest import fdrcorrection
import scipy.stats as ss
import os
from tqdm import tqdm


# regression functions for metabolic modeling results
from sklearn.linear_model import Lasso, LassoLars, ElasticNet, LarsCV, LassoCV, LassoLarsIC, ElasticNetCV, QuantileRegressor, LinearRegression
from sklearn.feature_selection import RFECV
from sklearn.svm import SVR


# packages for super learner models for regression
from math import sqrt
from numpy import hstack
from numpy import vstack
from numpy import asarray
from sklearn.datasets import make_regression
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import ElasticNet


# packages referred to "regressors" that is used to calculate pvalues for sklearn models
#import SCOOBI.regression_stats as regss


class regression_methods:
    
    def __init__(self, df_variables, df_response):
        self.df_var = df_variables
        self.df_res = df_response
    
    @staticmethod
    def selector_models(X, y, features, model_select='Coefficient'):
        # define a model
        # RFElinearRegression
        estimator = LinearRegression(positive=True, fit_intercept=False)
        selector = RFECV(estimator, cv=5, step=1)
        selector = selector.fit(X, y)
        if model_select=='Ranking':      
            return dict(zip(X.columns, selector.ranking_))
        else: # coefficient
            model = LinearRegression(positive=True, fit_intercept=False).fit(X[X.columns[selector.support_]], y)
            reg_coef = np.zeros(len(X.columns))
            reg_coef[selector.support_] = model.coef_
            return dict(zip(X.columns, reg_coef)), model, selector.support_
    
    @staticmethod
    def sklearn_models(X, y, model_select='LinearRegression'):
        # define a model
        if model_select=='LassoCV':
            #model = Lasso(
            #    positive=True, fit_intercept=False, random_state=0
            #).fit(X, y)
            model = LassoCV(
                cv=5, positive=True, fit_intercept=False, random_state=0
            ).fit(X, y)
        elif model_select=='ElasticNetCV':
            #model = ElasticNet(positive=True, fit_intercept=False, random_state=0).fit(X, y)
            model = ElasticNetCV(cv=5, positive=True, fit_intercept=False, random_state=0).fit(X, y)
        elif model_select=='LassoLarsIC':
            #model = LassoLars(
            #    positive=True, normalize=False, fit_intercept=False,
            #).fit(X, y)
            model = LassoLarsIC(
                criterion='aic', positive=True, normalize=False, fit_intercept=False,
            ).fit(X, y)
        elif model_select=='QuantileRegressor':
            model = QuantileRegressor(quantile=0.5, fit_intercept=False).fit(X, y)
        else:
            model = LinearRegression(positive=True, fit_intercept=False).fit(X, y)
        
        # get pvalue of the model
        #pv = regss.coef_pval(model, X, y)
        # temporily cancel pvalues cuz bugs
        return dict(zip(X.columns, model.coef_)), 1, model#dict(zip(X.columns, pv)), model
        

    def get_models(self, feature_sel, method_select, get_pvalues=False):

        model_coef = {}
        #for i in tqdm(range(len(self.df_res.columns))):
        for i in range(len(self.df_res.columns)):
            # data
            X, y = self.df_var, self.df_res.iloc[:, i]
            # model settings
            if feature_sel==True:
                model_coef[i] = self.selector_models(
                        X, y, len(X.columns)//2, model_select=method_select
                    )[0]
                
            else:
                if get_pvalues==True:
                    model_coef[i] = self.sklearn_models(
                            X, y, model_select=method_select
                            )[1]
                else:
                    model_coef[i] = self.sklearn_models(
                            X, y, model_select=method_select
                            )[0]
        
        return pd.DataFrame(model_coef)

    def coef_plot(self, df, binarize=False):

        fig, ax = plt.subplots(1,1,figsize=(20,5))
        if binarize==True:
            sns.heatmap(df>0, ax=ax, cmap='mako')
        else:
            sns.heatmap(df, ax=ax, cmap='mako')



    # collect out of fold predictions form k-fold cross validation
    def get_out_of_fold_predictions(self, X, y):
        meta_X, meta_y = list(), list()
        # define split of data
        kfold = KFold(n_splits=5, shuffle=True, random_state=0)
        # enumerate splits
        #i = 1
        for train_ix, test_ix in kfold.split(X):
            #if i==1:
            #    test_ix, train_ix = test_ix[:-1], np.append(test_ix[-1], train_ix)
            #    i = 0
            print(len(test_ix), len(train_ix))
            fold_yhats = list()
            # get data
            train_X, test_X = X.iloc[train_ix,:], X.iloc[test_ix,:]
            train_y, test_y = y.iloc[train_ix], y.iloc[test_ix]
            meta_y.extend(test_y)
            # model 1: linear regression with RFE
            # fit and make predictions with each sub-model
            _, model, features = self.selector_models(
                    train_X, train_y, len(X.columns)//2, model_select='Coefficient' # was using X, y
                    )
            yhat = model.predict(test_X[test_X.columns[features]])
            # store columns
            fold_yhats.append(yhat.reshape(len(yhat),1))

            # model 2: LASSO regression
            # fit and make predictions with each sub-model
            for model_sel in ['LassoCV', 'ElasticNetCV', 'LassoLarsIC']:
                model = self.sklearn_models(
                        train_X, train_y, model_select=model_sel    # was using X, y
                        )[2]
                yhat = model.predict(test_X)
                # store columns
                fold_yhats.append(yhat.reshape(len(yhat),1))

            # store fold yhats as columns
            meta_X.append(hstack(fold_yhats))

        return vstack(meta_X), asarray(meta_y)

    def SuperLearner(self):
        
        # check if the model is overfitting or not
        def meta_model_CV(model):
            from sklearn.model_selection import cross_val_score
            scores = cross_val_score(model, X, y, cv=5)
            return scores

        # fit a meta model
        def fit_meta_model(X, y):
            model = LinearRegression(positive=True, fit_intercept=False)
            model.fit(X, y)
            scores = meta_model_CV(model)
            # print(model.coef_)
            return model, scores
        
        # get models
        weights_lm = self.get_models(True, 'Coefficient', get_pvalues=False)
        model_res_list = [weights_lm]
        # get sklearn models
        for model_sel in ['LassoCV', 'ElasticNetCV', 'LassoLarsIC']:
            weights_lasso = self.get_models(False, model_sel, get_pvalues=False)
            model_res_list.append(weights_lasso)
        
        # fit meta-learner models
        meta_coef = {}
        scores_df = {}
        #for i in tqdm(range(len(self.df_res.columns))):
        for i in range(len(self.df_res.columns)):
            # data
            X, y = self.df_var, self.df_res.iloc[:, i]
            # get out of fold predictions
            meta_X, meta_y = self.get_out_of_fold_predictions(X, y)
            print('Meta ', meta_X.shape, meta_y.shape)
            # fit the meta model
            meta_model, scores = fit_meta_model(meta_X, meta_y)
            meta_coef[i] = meta_model.coef_
            scores_df[i] = scores

        
        # calculate the final coefficients by the new weight
        scores_df = pd.DataFrame(scores_df)
        print(scores_df.shape)
        print(weights_lm.shape)
        


        scores_df.columns = weights_lm.columns
        integrate_res = weights_lm.copy()
        for col in weights_lm.columns:
            integrate_res[col] = np.sum([model_w[col].mul(meta_coef[col][i]).to_numpy() for i, model_w in enumerate(model_res_list)], axis=0)
            #weights_lm[col].mul(meta_coef[col][0])+weights_lasso[col].mul(meta_coef[col][1])
        #RFElm_res = weights_lm
        #lasso_res = weights_lasso
        model_res_list.append(integrate_res)
        model_res_list.append(scores_df)



        return model_res_list#meta_coef, weights_lm, weights_lasso
    

    ## evaluate a list of models on a dataset
    #def evaluate_models(X, y, models):
    #    for model in models:
    #        yhat = model.predict(X)
    #        mse = mean_squared_error(y, yhat)
    #        print('%s: RMSE %.3f' % (model.__class__.__name__, sqrt(mse)))

    ## make predictions with stacked model
    #def super_learner_predictions(X, models, meta_model):
    #    meta_X = list()
    #    for model in models:
    #    yhat = model.predict(X)
    #    meta_X.append(yhat.reshape(len(yhat),1))
    #    meta_X = hstack(meta_X)
    #    # predict
    #    return meta_model.predict(meta_X)
