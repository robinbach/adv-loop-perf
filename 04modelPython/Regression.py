import numpy as np
import os, sys

def readFile(fPath):
	data = np.genfromtxt(fPath, delimiter=',')
	y = data.T[-1]
	numX = len(data.T) - 1
	A = data.T[0:numX].T
	return A, y

def main():
	dataPath = sys.argv[1]
	A, y = readFile(dataPath)
	print np.linalg.lstsq(A, y)[0]

if __name__ == '__main__':
	main()