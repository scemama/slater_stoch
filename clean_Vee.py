#!/usr/bin/env python3

import os
import numpy as np

def read_one(filename):
    with open(filename,'r') as f:
        l = []
        for line in f:
           try:
             i, j, x = line.replace('D','E').split()
           except:
             i, j, x, _ = line.replace('D','E').split()
           l.append( (int(i)-1, int(j)-1, float(x)) )

    ao_num = l[-1][0]+1

    V = np.zeros( (ao_num, ao_num) )
    for i, j, x in l:
        V[i,j] = x

    return V


def read_two(filename, ao_num):
    with open(filename, 'r') as f:
        V = np.zeros( (ao_num, ao_num, ao_num, ao_num) )
        for line in f:
            i, j, k, l, x, _ = line.replace('D','E').split()
            i, j, k, l = int(i)-1, int(j)-1, int(k)-1, int(l)-1
            if j<l or i<k:
               continue
            elif i == j and k<l:
               continue
            else:
               x = float(x)
#              if (i,j) != (k,l):
#                 if x > 0.: x = max(x-err, 0.)
#                 else: x = min(x+err, 0.)
#              if (i,j) == (k,l):
#                 if x > 0.: x = x+err
#                 else: x = x+err

               V[i, j, k, l] = x
               V[k, j, i, l] = x
               V[k, l, i, j] = x
               V[i, l, k, j] = x
               V[j, i, l, k] = x
               V[j, k, l, i] = x
               V[l, k, j, i] = x
               V[l, i, j, k] = x

    return(V)

def read_two_err(filename, ao_num):
    with open(filename, 'r') as f:
        V = np.zeros( (ao_num, ao_num, ao_num, ao_num) )
        for line in f:
            i, j, k, l, _, x = line.replace('D','E').split()
            i, j, k, l = int(i)-1, int(j)-1, int(k)-1, int(l)-1
            if j<l or i<k:
               continue
            elif i == j and k<l:
               continue
            else:
               x = float(x)
               V[i, j, k, l] = x
               V[k, j, i, l] = x
               V[k, l, i, j] = x
               V[i, l, k, j] = x
               V[j, i, l, k] = x
               V[j, k, l, i] = x
               V[l, k, j, i] = x
               V[l, i, j, k] = x

    return(V)

def write_two(filename, V):
    ao_num = V.shape[3]
    with open(filename, 'w') as f:
        for l in range(ao_num):
            for k in range(ao_num):
                for j in range(l,ao_num):
                    for i in range(k,ao_num):
                        if i == j and k<l:
                            continue
                        elif i >= j:
                            x = V[i,j,k,l]
                            if x != 0.:
                                f.write("%5d  %5d  %5d  %5d  %20.15E\n"% (i+1, j+1, k+1, l+1, x))



def rSVD(X,r,q,p):
    # Sample column space of X with P matrix
    ny = X.shape[1]
    P = np.random.randn(ny,r+p)
    Z = X @ P
    for k in range(q):
       Z = X @ (X.T @ Z)

    Q, _ = np.linalg.qr(Z, mode='reduced')

    # Compute SVD on projected Y
    Y = Q.T @ X
    UY, S, VT = np.linalg.svd(Y,full_matrices=0)
    U = Q @ UY

    return (U, S, VT)

def ldl(A):
    n = A.shape[0]
    L = np.eye  ( n )
    d = np.zeros( n )
    c = np.zeros( n )
    for j in range(n):
        print(j,'/',n)
        Ld = L[j,:j] * d[:j]
        s = np.sum(L[j,:j] * Ld)
        d[j] = A[j,j] - s
        c = 1./d[j]
        for i in range(j,n):
            s = np.sum(L[i,:j] * Ld)
            L[i,j] = (A[i,j] - s) * c
    return L, d


def main() :
    S   = read_one('S.qp')
    ao_num = S.shape[0]
    
    print('Read')
    Vee = read_two('bielec_ao', ao_num)
    Vee.shape = (ao_num*ao_num, ao_num*ao_num)

    algo = 5
    if algo == 1:
        print('Eigh')
        d, U = np.linalg.eigh(Vee)
        print(d)
        before = sum(d)
        for k in range(d.shape[0]):
            if d[k] < 0.:
               d[k] = 0.        # -149.6135404692439
#              d[k] = -d[k]     # -149.6302730994687
        after = sum(d)
        print("Removed: %f"%(before-after))
        d = d * before/after
        print(d)
        Vee = U @ np.diag(d) @ U.T
        print('Done')
    
    elif algo == 2:
        print('SVD')
        U, d, Vt = np.linalg.svd(Vee)
        Vee = U @ np.diag(d) @ U.T
#       Vee = np.einsum('ik,k,jk->ij', U , d , U)
        # -149.4806751082805
#       Vee = Vt.T @ np.diag(d) @ Vt
        # -149.4806751076974
        print('Done')
    
    elif algo == 3:
        print('rSVD')
        U, d, Vt = rSVD(Vee,1,30,ao_num*ao_num)
        Vee = U @ np.diag(d) @ U.T
#       Vee = np.einsum('ik,k,jk->ij', U , d , U)
        print('Done')
    
    elif algo == 4:
        print('L.D.Lt')
        L, d = ldl(Vee)
        print(d)
        before = sum(d)
        for k in range(d.shape[0]):
            if d[k] < 0.:
               d[k] = 0.
        after = sum(d)
        print("Removed: %f"%(before-after))
        d = d * before/after
        print(d)
        Vee = L @ np.diag(d) @ L.T
#       Vee = np.einsum('ik,k,jk->ij', L , d , L)
        #
        print('Done')
      
    elif algo == 5:

        print('V @ V.T')
        X = Vee @ Vee.T
        print('Eigh(V @ V.T)')
        d, U = np.linalg.eigh(X)
        print(d)
        for k in range(d.shape[0]):
            if d[k] < 0.:
               d[k] = 0.
        d = np.sqrt(d)
        print(d)
        Vee = U @ np.diag(d) @ U.T
#       Vee = np.einsum('ik,k,jk->ij', U , d , U)
        # -149.4806886397564
        print('Done')

    Vee.shape = (ao_num, ao_num, ao_num, ao_num)
    write_two('W.qp', Vee)



main()
