# definition d'une fonction donnant les indices des ancetres dans la redistribution multinomiale
import numpy as np
import numpy.random as rnd

def resampling_multi(w,N):
    u_tild = np.zeros((N))
    expo = np.zeros((N))
    alpha = np.zeros((N))
    u_ord = np.zeros((N))
    uu = np.zeros((N+1))
    s = np.zeros((N))
#
    w = w/w.sum()
    s = np.cumsum(w)
    u_tild = rnd.uniform(0,1,N)
#
    for i in range(N):
        alpha[i] = u_tild[i]**(1/float(i+1))
    alpha = np.cumprod(alpha)
    u_ord = alpha[N-1]/alpha
    u = np.append(u_ord,float("inf"))
#
    ancestor = np.zeros(N,dtype=int)
    offsprings = np.zeros(N,dtype=int)
    i = 0
    for j in range(N):
        o = 0
        while u[i]<=s[j]:
            ancestor[i] = j
            i = i+1
            o = o+1
        offsprings[j] = o
    return ancestor
