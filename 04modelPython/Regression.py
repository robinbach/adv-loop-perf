from sklearn import svm
from sklearn import linear_model
from sklearn.kernel_ridge import KernelRidge
import numpy as np
import sys
import random

import matplotlib.pyplot as plt

numTrain = 11

def readFile(fPath):
	data = np.genfromtxt(fPath, delimiter=',')
	random.shuffle(data)
	performance = data.T[-2]
	distortion = data.T[-1]
	numX = len(data.T) - 2
	A = data.T[0:numX]
	for i in range(len(A)):
		A[i] = A[i] / max(max(A[i]), 1.0)
	A = A.T
	ATrain = A[0:numTrain]
	ATest = A[numTrain + 1:]
	performanceTrain = performance[0:numTrain]
	performanceTest = performance[numTrain + 1:]
	distortionTrain = distortion[0:numTrain]
	distortionTest = distortion[numTrain + 1:]
	return ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest

def linearRegression(ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest):
	lr = linear_model.LinearRegression()
	lr.fit(ATrain, performanceTrain)
	performancePred = lr.predict(ATest)
	performanceErr = sum(abs(performancePred - performanceTest)) / len(performanceTest)
	print 'linear regression performance error: ', performanceErr
	
	lr.fit(ATrain, distortionTrain)
	distortionPred = lr.predict(ATest)
	distortionErr = sum(abs(distortionPred - distortionTest)) / len(distortionTest)
	print 'linear regression distortion error: ', distortionErr
	histoPlot(performancePred, performanceTest)
	histoPlot(distortionPred, distortionTest)

def SVR(ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest):
	clf = svm.SVR(C=100, epsilon=0.001)
	clf.fit(ATrain, performanceTrain)
	performancePred = clf.predict(ATest)

	performanceErr = sum(abs(performancePred - performanceTest)) / len(performanceTest)
	print 'SVR performance error: ', performanceErr
	clf.fit(ATrain, distortionTrain)
	distortionPred = clf.predict(ATest)
	distortionErr = sum(abs(distortionPred - distortionTest)) / len(distortionTest)
	print 'SVR distortion error: ', distortionErr

	histoPlot(performancePred, performanceTest)
	histoPlot(distortionPred, distortionTest)

def ridgeRegression(ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest):
	model = KernelRidge(alpha=0.01, kernel='sigmoid')
	model.fit(ATrain, performanceTrain)
	performancePred = model.predict(ATest)
	performanceErr = sum(abs(performancePred - performanceTest)) / len(performanceTest)
	print 'Kernel ridge performance error: ', performanceErr
	model.fit(ATrain, distortionTrain)
	distortionPred = model.predict(ATest)
	distortionErr = sum(abs(distortionPred - distortionTest)) / len(distortionTest)
	print 'Kernel ridge distortion error: ', distortionErr
	histoPlot(performancePred, performanceTest)
	histoPlot(distortionPred, distortionTest)

def robustRegression(ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest):
	model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
	model_ransac.fit(ATrain, performanceTrain)
	model_ransac.predict(ATest)
	temp = model_ransac.predict(ATest)
	performancePred = []
	for data in temp:
		performancePred.append(data[0])
	model_ransac.fit(ATrain, distortionTrain)
	model_ransac.predict(ATest)
	temp = model_ransac.predict(ATest)
	distortionPred = []
	for data in temp:
		distortionPred.append(data[0])	
	histoPlot(performancePred, performanceTest)
	histoPlot(distortionPred, distortionTest)

def histoPlot(pred, actual):
	x = np.arange(len(actual))
	plt.hold(True)
	rects1 = plt.bar(x, pred, 0.2, color='r')
	x = x + 0.2
	rects2 = plt.bar(x, actual, 0.2)
	plt.legend((rects1[0], rects2[0]), ('Prediction', 'Actual'), fontsize=20)
	plt.xlabel('Data Point', fontsize=30)
	plt.ylabel('Value', fontsize=30)
	performanceErr = sum(abs(pred - actual)) / len(actual)
	print 'Error: ', performanceErr	
	
	plt.title('Mean error: ' + ('%.3f' % performanceErr), fontsize=30)
	plt.hold(False)
	plt.show()

def main():
	dataPath = sys.argv[1]
	ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest = readFile(dataPath)
	linearRegression(ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest)
	SVR(ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest)
	ridgeRegression(ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest)
	robustRegression(ATrain, performanceTrain, distortionTrain, ATest, performanceTest, distortionTest)
if __name__ == '__main__':
	main()
	