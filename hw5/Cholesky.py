#!/usr/bin/env python3
import numpy as np

def Cholesky(A):
    if len(np.shape(A)) != 2 or np.shape(A)[0]!=np.shape(A)[1]:
        print("error shape")
        return
    
    for i in range(np.shape(A)[0]):
        for j in range(np.shape(A)[1]):
            if A[i][j] != A[j][i]:
                print("The input matrix should be symmetric")
                return 
       
    n = A.shape[0]
    l = np.eye(n)
    d = np.zeros((n,n))
    for k in range(n):
        if k == 0:
            d[k][k] = A[k][k]
            if d[k][k] == 0:
                print('error matrix type with 0 sequential principal minor determinant')
                return
            for m in range(n):
                l[m][k] = A[m][k]/d[k][k]
        else:
            temp_sum1 = 0
            for m in range(k):
                temp_sum1 += A[k][m]*l[k][m]
            d[k][k] = A[k][k] - temp_sum1
            for j in range(k+1,n):
                temp_sum2 = 0
                for m in range(k):
                    temp_sum2 += A[j][m]*l[k][m]
                A[j][k] = A[j][k] - temp_sum2
                l[j][k] = A[j][k]/d[k][k]
    return l,d


if __name__ == "__main__":
    A = np.array([[34, 47, 5, 18, 26], [47, 10, 13, 26, 34], [5, 13, 26, 39, 47],[18, 26, 39, 42, 5],[26, 34, 47, 5, 18]])
    L, D = Cholesky(A)
    print(L)
    print(D)
    print(L.T)