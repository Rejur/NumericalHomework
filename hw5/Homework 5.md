## Homework 5

###### 1、

_solution_:

```python
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
```

when $b=b_1$, the answer is $[0.02301282051282051, 0.01724358974358974, 0.2887179487179487, -0.09814102564102564, 0.046089743589743604]$

、$b=b_2$, the answer is $[0.025897435897435886, 0.0028205128205128164, 0.18871794871794878, -0.05871794871794883, 0.11820512820512828]$

and $b=b_3$, the answer is $[0.04698717948717948, 0.04506410256410258, 0.24666666666666665, 0.006602564102564129, 0.054679487179487174]$

###### 2、

_solution_:

```python
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
```

The answer is ![image-20191122170649308](/Users/hulin/Library/Application Support/typora-user-images/image-20191122170649308.png).

###### 3、

_solution_:

When $A=A_1$:

(1) $b = b$: the answer is $x = [0.19999999999999998, 0.2]$

(2)$b=\bar{b}$: the answer is $x=[0.19400000000000006, 0.20399999999999996]$

When $A=A_2$:

(1)$b=b$: the answer is $x = [0.5000000000000278, 0.49999999999997224]$

(2)$b=\bar{b}$: the answer is $x = [10.49500000000056, -9.50500000000056]$

When A isn't ill-condition. $A_1^{-1}=[[0.8, -0.4], [-0.2, 0.6]]$, so $||A_1||_1=3,||A_1^{-1}||_1=0.6,cond(A)=1.8$, so $\delta x$ is little.

And when A is ill-condition. $A_2^{-1}=[[-999., 1000.], [1001., -1000.]]$, so $||A_2||_1=2.001,||A_2^{-1}||_1=2000,cond(A)=4000>>1$. Because $\frac{||\delta x||}{||x||}\approx cond(A)(\frac{||\delta b||}{||b||}+\frac{||\delta A||}{||A||})$, so the $x$ change a lot.

