#!/usr/bin/env python3

import aleph
import numpy
import sys

from sklearn import svm
from sklearn import model_selection

from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_recall_curve

def kernel_1(I1, I2):
  return -abs(I1-I2).integral

def kernel_2(I1, I2):
  return -abs(I1-I2).pow(2).integral

def run_classification(K,C):
  clf    = svm.SVC(kernel='precomputed', probability=True)
  ss     = model_selection.ShuffleSplit(n_splits=10, test_size=0.5)
  scores = list()

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
    score       = accuracy_score(cTest, prediction)
    scores.append(score)

    p,r,t = precision_recall_curve(cTest, clf.decision_function(kTest))

    for x,y in zip(p,r):
      print("%f\t%f" % (x,y))

    print("\n")

  print("Average accuracy:", sum(scores)/len(scores), file=sys.stderr)

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
