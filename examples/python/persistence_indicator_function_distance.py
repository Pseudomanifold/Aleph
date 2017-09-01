#!/usr/bin/env python3

import aleph
import numpy
import sys

from sklearn import svm
from sklearn import model_selection

from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score

def kernel_1(I1, I2):
  return -abs(I1-I2).integral

def kernel_2(I1, I2):
  return -abs(I1-I2).pow(2).integral

def run_classification(K,C):
  clf              = svm.SVC(kernel='precomputed', probability=True)
  ss               = model_selection.ShuffleSplit(n_splits=10, test_size=0.5)
  accuracy_scores  = list()
  precision_scores = list()
  recall_scores    = list()

  for train, test in ss.split(K):
    kTrain = K[train][:,train]
    cTrain = C[train]

    # For the test kernel matrix, the values between *all* training
    # vectors and all test vectors must be provided.
    #kTest  = K[train, test]
    kTest = K[test]
    kTest = kTest[:,train]

    cTest = C[test]

    clf.fit(kTrain, cTrain)
    prediction  = clf.predict(kTest)

    accuracy_scores.append( accuracy_score(cTest, prediction) )
    precision_scores.append( precision_score(cTest, prediction) )
    recall_scores.append( recall_score(cTest, prediction) )

  print("Average accuracy: ", sum(accuracy_scores)/len(accuracy_scores), file=sys.stderr)
  print("Average precision:", sum(precision_scores)/len(precision_scores), file=sys.stderr)
  print("Average recall:   ", sum(recall_scores)/len(recall_scores), file=sys.stderr)

if __name__ == "__main__":
  persistenceDiagrams = []

  for filename in sys.argv[1:]:
    diagram = aleph.load_persistence_diagram(filename)

    diagram.removeDiagonal()
    diagram.removeUnpaired()

    persistenceDiagrams.append( diagram )

  persistenceIndicatorFunctions = []

  for diagram in persistenceDiagrams:
    persistenceIndicatorFunctions.append( aleph.make_persistence_indicator_function( diagram ) )

  n  = len(persistenceIndicatorFunctions)
  K1 = numpy.zeros((n,n), dtype=float)
  K2 = numpy.zeros((n,n), dtype=float)

  for i in range(0,n):
    f = persistenceIndicatorFunctions[i]
    for j in range(i+1,n):
      g       = persistenceIndicatorFunctions[j]
      K1[i,j] = kernel_1(f,g)
      K1[j,i] = K1[i,j]
      K2[i,j] = kernel_2(f,g)
      K2[j,i] = K2[i,j]

  # TODO: make configurable
  labels = [1]*50 + [-1]*50
  labels = numpy.array(labels)

  run_classification(K1, labels)
  run_classification(K2, labels)
