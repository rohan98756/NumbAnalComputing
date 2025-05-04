import numpy as np

L = np.zeros((5, 5))
U = np.zeros((5, 5))


A = np.array([[14.0, 14.0, -9.0, 3.0, -5.0],
              [14.0, 52.0, -15.0, 2.0, -32.0],
              [-9.0, -15.0, 36.0, -5.0, 16.0],
              [3.0, 2.0, -5.0, 47.0, 49.0],
              [-5.0, -32.0, 16.0, 49.0, 79.0]])

cond_A = np.linalg.cond(A)
print("Condition Number of A:", cond_A)

b = np.array([[-15.0], [-100.0], [106.0], [329.0], [463.0]])

n = 5

x = np.zeros(n)
y = np.zeros(n)

def backSubForUpper(n, u, c):
    for j in range(n-1, -1, -1):
        if u[j, j] == 0:
            return None
        x[j] = c[j] / u[j,j]
        for i in range(0, j):
            c[i] -= u[i, j] * x[j]
    return x

def forwardSubsForLower(n, l, c):
    for j in range(0,n):
        if l[j, j] == 0:
            return None
        y[j] = c[j] / l[j,j]
        for i in range(j+1, n):
            c[i] -= l[i, j] * y[j]
    return y

def computeLUInPlace(A):
    for k in range(0,n-1): 
        if A[k, k] == 0:
            return None
        for i in range(k+1, n):
            A[i, k] = A[i, k] / A[k, k]
        for j in range(k+1, n):
            for i in range(k+1, n):
                A[i, j] -= A[i, k] * A[k, j]
    return A



A = computeLUInPlace(A)

# Extract L and U from A
for i in range(n):
    L[i, i] = 1  
    L[i+1:, i] = A[i+1:, i]  
    U[i, i:] = A[i, i:]  


print("Here is L: ")
print(L)
print("\n")

print("Here is U: ")
print(U)
print("\n")

y = forwardSubsForLower(n, L, b)
print("\n")
print(y)

x = backSubForUpper(n, U, y)
print("\n")
print(x)




