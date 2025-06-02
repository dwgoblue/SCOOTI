"""
ADALINEMetaLearner.py
=================================================================================
Version 3.0: A collection of regression models for meta-learner regression models

-----
Major changes from V2.0:
    - Optimize the code with ChatGPT
    - Move the meta-learner model from linear regression to ADALINE
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from sklearn.linear_model import (
    LassoCV, ElasticNetCV, LassoLarsIC, LinearRegression
)
from sklearn.feature_selection import RFECV
from sklearn.model_selection import KFold, cross_val_score
from sklearn.metrics import mean_squared_error
from sklearn.exceptions import ConvergenceWarning
from tqdm import tqdm

warnings.simplefilter('ignore', ConvergenceWarning)
sns.set_style("whitegrid")


class Adaline:
    def __init__(self, learning_rate=0.001, epochs=10000):
        self.learning_rate = learning_rate
        self.epochs = epochs

    def fit(self, X, y):
        # inittiate the weights with random array
        self.weights = np.random.random(X.shape[1]+1)
        #self.weights = np.zeros(X.shape[1] + 1)

        stop = 0.001
        Error=[stop +1]
        step = 0
        # check the stop condition for the network
        while (Error[-1] > stop or Error[-1]-Error[-2] > 0.0001) and step<self.epochs:
            step += 1
        #for _ in range(self.epochs):
            output = self.net_input(X)
            errors = y - output
            self.weights[1:] += self.learning_rate * X.T.dot(errors)
            self.weights[0] += self.learning_rate * errors.sum()
            # Store sum of square errors
            Error.append((errors**2).sum())    
        print('Error :',Error[-1], Error[1])

    #def fit(self, X, y):
    #    self.weight = np.random.random(X.shape[1]+1)
    #    stop = 0.001
    #    Error=[stop +1]
    #    # check the stop condition for the network
    #    while Error[-1] > stop or Error[-1]-Error[-2] > 0.0001:
    #        error = []
    #        for i in range(X.shape[0]):
    #            Y_input = sum(weight*Input[i]) + bias
    #             
    #            # Update the weight
    #            for j in range(Input.shape[1]):
    #                weight[j]=weight[j] + lr*(Target[i]-Y_input)*Input[i][j]
 
    #            # Update the bias
    #            bias=bias + lr*(Target[i]-Y_input)
    #             
    #            # Store squared error value
    #            error.append((Target[i]-Y_input)**2)
    #        # Store sum of square errors
    #        Error.append(sum(error))
    #        print('Error :',Error[-1])
    #    return weight, bias


    def net_input(self, X):
        return np.dot(X, self.weights[1:]) + self.weights[0]

    def predict(self, X):
        return self.net_input(X)


class regression_methods:
    def __init__(self, df_variables, df_response, cluster_label_path=''):
        self.df_var = df_variables
        self.df_res = df_response

        if cluster_label_path:
            cluster_labels = pd.read_csv(cluster_label_path)
            self.df_var = self.df_var.loc[:, self.df_var.columns.isin(cluster_labels.iloc[:, 0])]
            self.cluster_labels = cluster_labels['cluster'].values
            self.cluster_label_names = cluster_labels.iloc[:, 0].values
        else:
            self.cluster_labels = []
            self.cluster_label_names = []

        self.models = none
        self.meta_model = none

    @staticmethod
    def selector_models(X, y):
        """Integrate feature selections into regression models."""
        # define a linear regression model
        estimator = LinearRegression(positive=True, fit_intercept=False)
        # RFE feature selection
        selector = RFECV(estimator, cv=5, step=1)
        # fit a model
        selector = selector.fit(X, y)
        # coefficient
        Xtmp = X.copy()
        # make unselected features 0
        Xtmp[Xtmp.columns[selector.support_]] = 0
        # fit the model again with modified variables
        model = LinearRegression(
                positive=True, # non-negative coefficients
                fit_intercept=False
                ).fit(Xtmp, y)
        return model

    @staticmethod
    def sklearn_models(X, y, model_select='LinearRegression'):
        """Build regression models with different loss functions."""
        # define models
        if model_select=='LassoCV':
            # LASSO
            model = LassoCV(
                cv=5, positive=True, fit_intercept=False, random_state=0
            ).fit(X, y)
        elif model_select=='ElasticNetCV':
            # ElasticNet
            model = ElasticNetCV(
                    cv=5, positive=True, fit_intercept=False, random_state=0
                    ).fit(X, y)
        elif model_select=='LassoLarsIC':
            # LASSO-Lars
            model = LassoLarsIC(
                criterion='aic', positive=True, normalize=False, fit_intercept=False,
            ).fit(X, y)
        else:
            # typical linear regression without feature selections
            model = LinearRegression(positive=True, fit_intercept=False).fit(X, y)

        return model

    def get_models(self, X, y, return_dict=False):
        models = {}
        models['rfelm'] = self.selector_models(X, y)
        for model_sel, key in zip(['LassoCV', 'ElasticNetCV', 'LassoLarsIC'], ['lasso', 'EN', 'LL']):
            models[key] = self.sklearn_models(X, y, model_select=model_sel)
        return models if return_dict else list(models.values())

    def get_out_of_fold_predictions(self, X, y):
        meta_X, meta_y = [], []
        kfold = KFold(n_splits=5, shuffle=True, random_state=0)
        for train_ix, test_ix in kfold.split(X):
            fold_yhats = []
            train_X, test_X = X.iloc[train_ix, :], X.iloc[test_ix, :]
            train_y, test_y = y.iloc[train_ix], y.iloc[test_ix]
            meta_y.extend(test_y)
            models = self.get_models(train_X, train_y)
            for model in models:
                yhat = model.predict(test_X)
                fold_yhats.append(yhat.reshape(-1, 1))
            meta_X.append(np.hstack(fold_yhats))

        return np.vstack(meta_X), np.asarray(meta_y)

    def fit_meta_model(self, X, y):
        print('Fitting Adaline meta model...')
        if self.learner=='L':
            print('fitting meta model...')
            model = LinearRegression(positive=True, fit_intercept=False)
            model.fit(X, y)

        else: # ADALINE
            meta_model = Adaline()
            meta_model.fit(X, y)
            self.meta_model = meta_model

        # 5 fold CV
        kfold = KFold(n_splits=5, shuffle=True, random_state=1000)
        scores = []
        for train_ix, test_ix in kfold.split(X):
            train_X, test_X = X[train_ix, :], X[test_ix, :]
            train_y, test_y = y[train_ix], y[test_ix]
            y_pred = meta_model.predict(test_X)
            mse = mean_squared_error(test_y, y_pred)
            scores.append(mse)  # Just to keep consistent with the structure
        #print(f'Mean Squared Error: {mse}')
        return meta_model, scores

    def fit_base_models(self):
        res_collect = []
        df_res_tmp = pd.DataFrame(self.df_res)
        for i in range(len(df_res_tmp.columns)):
            X, y = self.df_var, df_res_tmp.iloc[:, i]
            model_dict = self.get_models(X, y)
            res_collect.append(model_dict)

        model_coef_collect = []
        for j in range(len(model_dict)):
            model_coef = {}
            for i in range(len(res_collect)):
                model_coef[i] = res_collect[i][j].coef_
            model_coef_collect.append(pd.DataFrame(model_coef))

        self.models = res_collect
        return model_coef_collect

    def evaluate_models(self, models):
        for model in models:
            yhat = model.predict(self.df_var)
            mse = mean_squared_error(self.df_res, yhat)
            print(f'{model.__class__.__name__}: RMSE {np.sqrt(mse):.3f}')

    def super_learner_predictions(self, X, models, meta_model):
        meta_X = [model.predict(X).reshape(-1, 1) for model in models]
        meta_X = np.hstack(meta_X)
        return meta_model.predict(meta_X)

    def clusterSuperLearner(self):
        # save the original data
        ori_var = self.df_var.copy()
        ori_res = self.df_res.copy()
        cluster_model_res_list = []
        # list to collect cluster superlearners
        for cluster in np.unique(self.cluster_labels):
            # make predictions based on training data
            self.df_var = ori_var.loc[:, self.cluster_labels == cluster]
            # fit superlearner model
            model_res_list = self.super_learner()
            # collect coefficients
            cluster_model_res_list.append(model_res_list[-1])
            # replace the variables with the original table
            self.df_var, self.df_res = ori_var.copy(), ori_res.copy()

        meta_coef = {}
        scores_collect = {}

        for i in tqdm(range(len(self.df_res.columns))):
            meta_X, meta_y = [], []
            kfold = KFold(n_splits=5, shuffle=True, random_state=0)
            X, y = self.df_var.copy(), self.df_res.iloc[:, i].copy()

            for train_ix, test_ix in kfold.split(X):
                fold_yhats = []
                train_X, test_X = X.iloc[train_ix, :], X.iloc[test_ix, :]
                train_y, test_y = y.iloc[train_ix], y.iloc[test_ix]
                meta_y.extend(test_y)

                for cluster in np.unique(self.cluster_labels):
                    self.df_var = train_X.loc[:, self.cluster_labels == cluster]
                    self.df_res = train_y
                    model_res_list = self.fit_base_models()
                    model_res_list = self.super_learner()
                    cluster_test_X = test_X.loc[:, self.cluster_labels == cluster]
                    yhat = self.super_learner_predictions(cluster_test_X, self.models[0], self.meta_model)
                    fold_yhats.append(yhat.reshape(-1, 1))

                meta_X.append(np.hstack(fold_yhats))

            meta_X, meta_y = np.vstack(meta_X), np.asarray(meta_y)
            self.df_var, self.df_res = ori_var.copy(), ori_res.copy()
            meta_model, scores = self.fit_meta_model(meta_X, meta_y)
            meta_coef[i] = meta_model.weights[1:]  # Adaline does not have a coef_ attribute
            scores_collect[i] = scores

        meta_df = pd.DataFrame(meta_coef)
        integrate_res = []
        for i in range(len(cluster_model_res_list)):
            integrate_res.append(cluster_model_res_list[i].mul(meta_df.iloc[i, :]))
        integrate_res = pd.concat(integrate_res, axis=0)

        cluster_model_res_list.append(integrate_res)
        cluster_model_res_list.append(meta_df)

        return cluster_model_res_list

    def SuperLearner(self):
        # get base models
        model_res_list = self.fit_base_models()
        # containers to save results
        meta_coef = {}
        scores_collect = {}
        meta_model_collect = []
        df_res_tmp = pd.DataFrame(self.df_res)
        # fit models
        for i in range(len(df_res_tmp.columns)):
            X, y = self.df_var, df_res_tmp.iloc[:, i]
            meta_X, meta_y = self.get_out_of_fold_predictions(X, y)
            meta_model, scores = self.fit_meta_model(meta_X, meta_y)
            meta_model_collect.append(meta_model)
            meta_coef[i] = meta_model.weights[1:]  # Adaline does not have a coef_ attribute
            scores_collect[i] = scores

        self.meta_model = meta_model_collect[0]

        scores_df = pd.DataFrame(scores_collect)
        integrate_res = model_res_list[0].copy()
        for col in model_res_list[0].columns:
            integrate_res[col] = np.sum(
                [model_w[col].mul(meta_coef[col][i]).to_numpy() for i, model_w in enumerate(model_res_list)], axis=0)
        model_res_list.append(integrate_res)

        for model_res in model_res_list:
            model_res.index = self.df_var.columns

        model_res_list.append(scores_df)
        return model_res_list

    #def SuperLearner(self):
    #    # get models
    #    model_res_list = self.fit_base_models()
    #    
    #    # fit meta-learner models
    #    meta_coef = {}
    #    scores_collect = {}
    #    meta_model_collect = []
    #    df_res_tmp = pd.DataFrame(self.df_res)
    #    #for i in tqdm(range(len(df_res_tmp.columns))):
    #    for i in range(len(df_res_tmp.columns)):
    #        # data
    #        X, y = self.df_var, df_res_tmp.iloc[:, i]
    #        # get out of fold predictions
    #        meta_X, meta_y = self.get_out_of_fold_predictions(X, y)
    #        #print('Meta ', meta_X.shape, meta_y.shape)
    #        # fit the meta model
    #        meta_model, scores = self.fit_meta_model(meta_X, meta_y)
    #        meta_model_collect.append(meta_model)
    #        meta_coef[i] = meta_model.coef_
    #        scores_collect[i] = scores

    #    
    #    # only works for one-column response
    #    self.meta_model = meta_model
    #    
    #    # calculate the final coefficients by the new weight
    #    scores_df = pd.DataFrame(scores_collect)
    #    print(scores_df.shape)
    #    print(model_res_list[0].shape)
    #    scores_df.columns = model_res_list[0].columns
    #    integrate_res = model_res_list[0].copy()
    #    for col in model_res_list[0].columns:
    #        integrate_res[col] = np.sum([model_w[col].mul(meta_coef[col][i]).to_numpy() for i, model_w in enumerate(model_res_list)], axis=0)
    #    model_res_list.append(integrate_res)
    #    # replace index with the metabolite name
    #    for i in range(len(model_res_list)):
    #       model_res_list[i].index = self.df_var.columns

    #    # append scores
    #    model_res_list.append(scores_df)

    #    return model_res_list

# Usage example:
# X, y = np.random.rand(100, 10), np.random.rand(100)
# reg_methods = RegressionMethods(df_variables=pd.DataFrame(X), df_response=pd.DataFrame(y))
# cluster_learners = reg_methods.cluster_super_learner()
