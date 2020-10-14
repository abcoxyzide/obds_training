#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 23:40:51 2018

@author: edward
"""

import numpy as np
from sklearn.datasets import load_iris
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

## Load Iris data
iris = load_iris()
print(iris.data)
print(iris.target)

### Create a classifer and predict===================
## Create classifier object
clf = DecisionTreeClassifier()
## Fit the data with the classifier. We are using the default parameters
clf = clf.fit(iris.data, iris.target)

## Predict a few value from iris data
clf.predict(iris.data[:2, :])
## Check this is correct
print(iris.target[:2]) 

## Create my own data
fake_data = np.array([[20, 31, 44, 11],
          [100, 200, 3 ,1]])
## Predict with the tree
clf.predict(fake_data)

## Is it over fitting? ===========================
data_train, data_test, labels_train, labels_test = train_test_split(iris.data, iris.target, test_size=0.25)

# First score
clf.score(data_test, labels_test)

clf = clf.fit(data_train, labels_train)
## Predict a few value from iris data
clf.score(data_test, labels_test)
## ===============================================

## Random forest classifier (overfits less)
clf_forest = RandomForestClassifier() # Default parameters, these often need tuning
clf_forest = clf_forest.fit(data_train, labels_train)
clf_forest.predict(iris.data[:2, :])
clf_forest.score(data_test, labels_test)







