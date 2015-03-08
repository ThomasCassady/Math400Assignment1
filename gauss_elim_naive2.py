import linear_algebra as la

__author__ = 'nathanael'


debugging = True
printing = True


def findPivotRow1(mat, col):
    #    Finds index of the first row with nonzero entry on or
    #    below diagonal.  If there isn't one return(-1).

    epsilon = 10**(-17)
    for row in range(col, la.rows(mat)):
        #        if mat[row][col] != 0:
        if abs(mat[row][col]) > epsilon:
            return(row)
    return(-1)


def rowReduce(M, type, row_maxes):
    #    return row reduced version of M

    N = la.copyMatrix(M)
    cs = la.cols(M)-2   # no need to consider last two cols
    rs = la.rows(M)
    for col in range(cs+1):
        if type == 1:
            j = findPivotRow1(N,col)
        elif type == 2:
            j = findPartialPivotRow(N,col)
            if debugging:
                print("pivot val = ", N[col][j])
        elif type == 3:
            j = findScaledPivotRow(N, col, row_maxes)
        if j < 0:
            print("\nrowReduce: No pivot found for column index %d "%(col))
            return(N)
        else:
            #---- printing matrix
            if printing:
                print("Matrix before reduction")
                la.show(N)
            # ----End Dprinting

            if j != col:
                N = la.swaprows(N,col,j)
            scale = -1.0 / N[col][col]
            for row in range(col+1,rs):
                N = la.addrows(N, col, row, scale * N[row][col])

            # --- check calculations
            if printing:
                print("Matrix after reduction")
                la.show(N)
                # --- End check

    return(N)


def backSub(M):

    #   given a row reduced augmented matrix with nonzero
    #   diagonal entries, returns a solution vector


    cs = la.cols(M)-1 # cols not counting augmented col
    sol = [0 for i in range(cs)] # place for solution vector
    for i in range(1,cs+1):
        row = cs-i # work backwards
        sol[row] = ((M[row][cs] - sum([M[row][j]*sol[j] for
                                       j in range(row+1,cs)])) / M[row][row])
    return(sol)


def diag_test(mat):

    #   Returns True if no diagonal element is zero, False
    #   otherwise.


    for row in range(la.rows(mat)):
        if mat[row][row] == 0:
            return(False)
    else:
        return(True)

''' modified: added type '''

def ge_1(aug, pivot_type):

    #   Given an augmented matrix it returns a list.  The [0]
    #   element is the row-reduced augmented matrix, and
    #   ge_1(aug)[1] is a solution vector.  The solution
    #   vector is empty if there is no unique solution.
    scaled_pivot = 3
    row_maxes = []
    '''
    If scaled pivoting get row scalings
    '''
    if pivot_type == scaled_pivot:
        row_maxes = findRowScalings(aug)
    aug_n = rowReduce(aug, pivot_type, row_maxes)
    if diag_test(aug_n):
        sol = backSub(aug_n)
    else:
        print("\nge_1(): There is no unique solution")
        sol = []
    results = [aug_n, sol]
    return(results)


# # The next two functions support checking a solution.

# In[14]:

def getAandb(aug):
    #   Returns the coef. matrix A and the vector b of Ax=b
    m = la.rows(aug)
    n = la.cols(aug)
    A = la.zero(m,n-1)
    b = la.zero(m,1)
    for i in range(m):
        for j in range(n-1):
            A[i][j] = aug[i][j]

    for i in range(m):
        b[i] = aug[i][n-1]
    Aandb = [A,b]
    return(Aandb)

def checkSol_1(aug,x):
    #   For aug=[A|b], returns Ax, b, and b-Ax as vectors
    A  = getAandb(aug)[0]
    b  = getAandb(aug)[1]
    x_col_vec = la.vec2colVec(x)
    Ax = la.matMult(A,x_col_vec)
    r  = la.addVectors(b, la.scalarMult(-1.0, la.colVec2vec(Ax)))
    L  = [Ax,b,r]
    return(L)




def LUfactorization2(mat, type):
    #TODO: currently doesn't handle
    # row permutations.
    ''' Function returns the the matrices
    [PL, U] where PL combines the permutation
    matrix and the lower triangular factor matrix
    and U is the upper triangular factor matrix.

    '''

    len_mat = len(mat)
    #scaled_pivot = 3
    #row_maxes = []
    Ab = getAandb(mat)
    A = Ab[0]
    if debugging:
        print("first A")
        la.show(A)
    # permute vector determines permute vector
    permute_vector= []
    for i in range(len_mat):
        permute_vector.append(i)

    L = getZeroMatrix(len_mat)
    U = getZeroMatrix(len_mat)

    for i in range(len_mat):
        L[i][i] = 1
    # ----
    for j in range(len_mat):
        U[0][j] = A[0][j]

    for col in range(0,len_mat -1):
        j = findPivotRow1(A,col)
        if (j != col):
            swapItems(j, col, permute_vector)
            A = la.swaprows(A,col,j)

        for row in range(col+1, len_mat):
            scale = -1.0 / A[col][col]
            #print("here", row, col)
            L[row][col] = (-scale) * A[row][col]
            print("scale", scale)
            A = la.addrows(A, col, row, scale * A[row][col])
            print("A is")
            la.show(A)

    return [L,A]



#--- help function for swapping
def swapItems(i,j, vec):
    vec[i] = j
    vec[j] = i
    return vec


def getZeroMatrix(N):
    row = [0]*N
    mat = []
    for i in row:
        mat.append(list(row))
    return mat




#----- Additional functions for partial and scaled Pivoting

def findRowScalings(mat):
    ''' Function returns vector of values representing
    the max vals in each row with index values corresponding
    to row numbers'''
    max_vals = []
    max_val = 0
    for row in mat:
        for i in row:
            curr_val = abs(i)
            #  print("curr_val ", curr_val)
            if curr_val > max_val:
                max_val = curr_val
        #print("max_val", max_val)
        max_vals.append(max_val)
        max_val = 0
    return (max_vals)


def findScaledPivotRow(mat,col, row_maxes):
    #    Finds index of the first row with nonzero entry on or
    #    below diagonal.  If there isn't one return(-1).

    epsilon = 10**(-17)
    "initialzie max_val to epsilon"
    max_val = epsilon
    max_row = col

    for row in range(col, la.rows(mat)):
        "additions for partial pivot"
        temp_val = abs(mat[row][col])/row_maxes[col]
        if temp_val > max_val:
            max_val = temp_val
            max_row = col
    '''
    max_row is row with max_val -- return if
    max_row if max_val > epsilon
    '''
    if debugging:
        print ("max_val", max_val)
    if max_val > epsilon:
        return(max_row)
    return(-1)




def findPartialPivotRow(mat,col):
    #    Finds index of the first row with nonzero entry on or
    #    below diagonal.  If there isn't one return(-1).

    epsilon = 10**(-17)
    "initialzie max_val to epsilon"
    max_val = epsilon
    max_row = col

    for row in range(col, la.rows(mat)):
        "additions for partial pivot"
        temp_val = abs(mat[row][col])
        if temp_val > max_val:
            max_val = temp_val
            max_row = col
    '''
    max_row is row with max_val -- return if
    max_row if max_val > epsilon
    '''
    if debugging:
        print ("max_val", max_val)
    if max_val > epsilon:
        return(max_row)
    return(-1)

