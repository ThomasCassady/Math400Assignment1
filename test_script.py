import copy
import linear_algebra as la
# import hilbert_1 as hil
import gauss_elim_naive2 as ge


__author__ = 'nathanael'

'''
Basic test scripts
Pivot Type Key

Basic pivot: basic_pivot = = 1
Partial Pivot : partial_pivot = 2
Partial Scaled Pivot : part_scale == 3

'''

'''
Test functions - could be moved to separate file later
'''

def basic_test(aug_mat, pivot_type):
    '''
        Function takes an augmented matrix - aug_mat,
        runs Gaussian elimination with
        pivot_type, prints and prints the given
        error vector, error norm, and reduced
        form of input matrix- aug_matrix
    '''
    res_list = ge.ge_1(aug_mat, pivot_type)
    sol_check = ge.checkSol_1(aug_mat, res_list[1])

    print "error vector is ", sol_check[2]
    print "error vector Lp-max norm is:", max(sol_check[2])
    print "Reduce form of input matrix is", res_list[0]



'''
Testing script begin
'''
basic = 1
partial = 2
scaled = 3

basic_pivot = 1
partial_pivot = 2
scaled_pivot = 3


B = [[1, 3, -4, 6], [1, -2, 4, 3], [4, -5, 3, 2]]
aug_mat = B

pivot_type = partial_pivot

basic_test(aug_mat, pivot_type)

''' Testing for problem 2 '''

A = [[10, 10, 10, 10 ** 17, 10 ** 17],
     [1, 10 ** (-3), 10 ** (-3), 10 ** (-3), 1],
     [1, 1, 10 ** (-3), 10 ** (-3), 2],
     [1, 1, 1, 10 ** (-3), 3]]

AandB = ge.ge_1(A, scaled_pivot)

#----- checkSol returns List = Ax, b, r
ge.checkSol_1(A, AandB[1])


''' LU Factorization Test
'''


res = ge.LUfactorization2(B, pivot_type)
#res
L = res[0]
U = res[1]
#show(L)
#show(U)

C = la.matMult(L, U)
print "compare L * U to original Matrix "
la.show(C)
la.show(B)

#x = getZeroMatrix(4)
