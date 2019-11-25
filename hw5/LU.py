import numpy as np
 
def LU(A, n):
    L = [[0 for x in range(n)] for y in range(n)]
    U = [[0 for x in range(n)] for y in range(n)]  
    for i in range(n, 0, -1):
	    for tmp1 in range(n - i, n):
		    sum = 0
		    for tmp2 in range(n - i):
			    sum += L[tmp1][tmp2] * U[tmp2][n - i]
		    L[tmp1][n - i] = A[tmp1][n - i] - sum
	    for tmp1 in range(n - i + 1, n):
	    	sum = 0
	    	for tmp2 in range(n - i):
	    		sum += L[n - i][tmp2] * U[tmp2][tmp1]
	    	U[n - i][tmp1] = (A[n - i][tmp1] - sum) / float(L[n - i][n - i])
    for i in range(n):
    	U[i][i] = 1
    #print("L:")
    print(L, 3)	
    #print("U")
    print(U, 3)
    return L, U
    
def LUSolve(A, b, n):
	L, U = LU(A, n)
	Y = [0 for x in range(n)]
	X = [0 for x in range(n)]
	Y[0] = b[0][0] / L[0][0]
	for i in range(1, n):
		sum = 0
		for tmp in range(0, i):
			sum += L[i][tmp] * Y[tmp]
		Y[i] = (b[0][i] - sum) / float(L[i][i])
	X[n - 1] = Y[n - 1]
	for i in range(n - 2, -1, -1):
		sum = 0
		for tmp in range(i + 1, n):
			sum += U[i][tmp] * X[tmp]
		X[i] = Y[i] - sum
	return X
 
if __name__ == "__main__":
    M = np.array([[17, 24, 1, 8, 15], [23, 5, 7, 14, 16], [4, 6, 13, 20, 22], [10, 12, 19, 21, 3], [11, 18, 25, 2, 9]])
    #b = np.array([1, 2, 3, 4, 8])
    #b = np.array([2, 3, 4, 3, 6])
    b = np.array([3, 4, 5, 6, 8])
    #M = np.array([[3, 2], [1, 4]])
    #M = np.array([[1.000, 1.000],[1.001, 0.999]])
    #M = M / 2

    #b = np.array([1,1])
    #b = np.array([0.99, 1.01])
    #b = b / 2

    print("M:")
    print(M)
    A = M
    #print("A:")
    #print(A)
    n = M.shape[1]
    print(b)
    b = np.reshape(b, (1, n))
    result = LUSolve(A, b, n)
    print(result)

