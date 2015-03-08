from __future__ import division


__author__ = 'nathanael'


# coding: utf-8

# """
# gauss_elim_1.py
#
# Includes code for functions that do basic vector and
# matrix arithmetic.  Most of these functions support
# the function ge_1(aug) which takes an n by n+1
# augmented matrix and returns a row-reduced matrix and
# an approximate solution vector of the corresponding linear
# system.  It uses gaussian elimination with a naive pivot
# strategy.  That is, at each stage in the row reduction it
# chooses, as the pivot, the first nonzero entry that lies
# on or below the diagonal in the current column.
#
# revision 2.00 03/12/12  changed filename, commented out tests [ds]
# revision 2.00 02/17/15  fixed test for zero entry in findPivotrow1 [ds]
#
# matrix examples: [[28, 12, 20, 28], [4, 32, 28, 16]] , 2 by 4
# [[28, 12, 20, 28]] , 1 by 4
#                  [[[28], [12], [20], [28]] , 4 by 1
#
#
# vector example [28, 12, 20, 28]
#
# """

# # Supporting functions for Naive Gaussian Elimination function ge_1

# In[2]:



# In[3]:

def rows(mat):
    #   "return number of rows"
    return (len(mat))


def cols(mat):
    #    "return number of cols"
    return (len(mat[0]))


def zero(m, n):
    #    "Create zero matrix"
    new_mat = [[0 for col in range(n)] for row in range(m)]
    return new_mat


def transpose(mat):
    #    "return transpose of mat"
    new_mat = zero(cols(mat), rows(mat))
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            new_mat[col][row] = mat[row][col]
    return (new_mat)


def dot(A, B):
    #    "vector dot product"
    if len(A) != len(B):
        print("dot: list lengths do not match")
        return ()
    dot = 0
    for i in range(len(A)):
        dot = dot + A[i] * B[i]
    return (dot)


def getCol(mat, col):
    #    "return column col from matrix mat"
    return ([r[col] for r in mat])


def getRow(mat, row):
    #    "return row row from matrix mat"
    return (mat[row])


def matMult(mat1, mat2):
    #    "multiply two matrices"
    if cols(mat1) != rows(mat2):
        print("matMult: mismatched matrices")
        return ()
    prod = zero(rows(mat1), cols(mat2))
    for row in range(rows(mat1)):
        for col in range(cols(mat2)):
            prod[row][col] = dot(mat1[row], getCol(mat2, col))
    return (prod)


def vectorQ(V):
    #    "mild test to see if V is a vector"
    if type(V) != type([1]):
        return (False)
    if type(V[0]) == type([1]):
        return (False)
    return (True)


def scalarMult(a, mat):
    #    "multiply a scalar times a matrix"
    if vectorQ(mat):
        return ([a * m for m in mat])
    for row in range(rows(mat)):
        for col in range(cols(mat)):
            mat[row][col] = a * mat[row][col]
    return (mat)


def addVectors(A, B):
    #   "add two vectors"
    if len(A) != len(B):
        print("addVectors: different lengths")
        return ()
    return ([A[i] + B[i] for i in range(len(A))])


def swaprows(M, i, j):
    #    "swap rows i and j in matrix M"
    N = copyMatrix(M)
    T = N[i]
    N[i] = N[j]
    N[j] = T
    return N


def copyMatrix(M):
    return ([[M[row][col] for col in range(cols(M))] for row in
             range(rows(M))])


def addrows(M, f, t, scale=1):
    #    "add scale times row f to row t"
    N = copyMatrix(M)
    T = addVectors(scalarMult(scale, N[f]), N[t])
    N[t] = T
    return (N)


def show(mat):
    #    "Print out matrix"
    for row in mat:
        print(row)


# # vectors vs rowVectors and colVectors
# # the latter are matrices

# In[4]:

def vec2rowVec(vec):
    #    "[a,b,c] -> [[a,b,c]]"
    return ([vec])


def vec2colVec(vec):
    return (transpose(vec2rowVec(vec)))


def colVec2vec(mat):
    rowVec = transpose(mat)
    return (rowVec[0])


def augment(mat, vec):
    #    "given nxn mat and n length vector return augmented matrix"
    amat = []
    for row in range(rows(mat)):
        amat.append(mat[row] + [vec[row]])
    return (amat)

