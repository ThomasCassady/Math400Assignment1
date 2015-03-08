"""
hilbert_1.py

Includes code for functions that support basic vector and
matrix arithmetic.  It also includes functions that generate
Hilbert Matrices, and the inhomogeneous b-vector needed in
Assignment 1.

    Revision 1.00 03/12/13 ds

"""
from math import *


### Hilbert matrix and related functions

def hilbert_b(n):
    # creates "Hilbert" b vector
    b_vec = [0 for k in range(n)]
    b_vec[0] = 2/pi
    b_vec[1] = 1/pi
    for k in range(2,n):
        b_vec[k] = 1/pi -(k*(k-1)/pi**2)*b_vec[k-2]
    return(b_vec)

def hilbert(n):
    # creates a Hilbert matrix
    h_n = zero(n,n)
    for row in range(n):
        for col in range(n):
            h_n[col][row] = 1.0/(1.0+row+col)
    return(h_n)

def n_C_k(n,k):
    # computes the binomial coefficien nCk
    if k > n:
        n_ch_k = 0
    prod = 1.0
    for j in range(k):
        prod = prod*(1.0*(n-j)/(k-j))
    return(prod)

def hilbert_inv(n):
    # creates the inverse of a Hilbert matrix
    h_inv = zero(n,n)
    for k in range(n):
        for m in range(n):
            h_inv[k][m] = ((-1)**(k+m))*((k+m+1)*
                                          (n_C_k(n+k,n-m-1))*
                                          (n_C_k(n+m,n-k-1))*
s                                          (n_C_k(k+m,k))**2 )
    return(h_inv)




