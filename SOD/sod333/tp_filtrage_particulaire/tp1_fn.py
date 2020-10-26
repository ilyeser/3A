import matplotlib.pyplot as plt
import matplotlib.animation as anim
import scipy.io as sio
import numpy as np
import numpy.random as rnd
import math
from resampling_multi import resampling_multi

r0 = (-6000 , 2000)
v0 = (120 , 0)

X1MIN = -10000
X1MAX = 10000
X2MIN = -10000
X2MAX = 10000

sigma_r0 = 100
sigma_v0 = 10
sigma_INS = 7
sigma_ALT = 10
sigma_BAR = 20
sigma_V = np.sqrt(sigma_BAR**2 + sigma_ALT**2)
T = 100

delta = 1

M = np.eye(4)
M[0,2] = delta
M[1,3] = delta

map = sio.loadmat("mnt.mat")["map"]
N1 = map.shape[1]
N2 = map.shape[0]

traj = sio.loadmat("traj.mat")
rtrue = traj["rtrue"]
vtrue = traj["vtrue"]

#------------QUESTION 2
a_INS = sio.loadmat("ins.mat")["a_INS"]
nmax = 101

r_INS = np.zeros(rtrue.shape)
v_INS = np.zeros(vtrue.shape)

r_INS[: , 0] = r0
v_INS[: , 0] = v0

for k in range(1,T+1):
    r_INS[: , k] = r_INS[: , k-1] + delta*v_INS[: , k-1]
    v_INS[: , k] = v_INS[: , k-1] + delta*a_INS[: , k-1]

#-------------

#-------------QUESTIONS 4 et 5
def h(x):
    i = min(N2-1,max(0,math.floor(N2*(X2MAX - x[1])/(X2MAX - X2MIN))))
    j = min(N1-1,max(0,math.floor(N1*(x[0] - X1MIN)/(X1MAX - X1MIN))))
    return(map[i,j])


h_ALT = sio.loadmat('alt.mat')['h_ALT'][0]

H = np.zeros([T+1])
for k in range(T+1):
    H[k] = h(np.array([rtrue[0,k] , rtrue[1,k]]))


#--------------

#--------------QUESTION 7
def f(X,w):
    return (np.dot(M,X) - delta * np.array([0,0,w[0],w[1]]))

def hp(k,X):
    x1 = r_INS[0,k]+X[0]
    x2 = r_INS[1,k]+X[1]
    return (h(np.array([x1,x2])))

def q_k(v):
    return (np.exp(-0.5*(v/sigma_V)**2))
#--------SIR
def filtrage_SIR(N):
    ksi = np.zeros([4,N])
    ksi_tilde = np.zeros([4,N])
    particules = np.zeros([4,N,T+1])
    moyenne = np.zeros([2,T+1])
    
    w = np.zeros([N])
    poids = np.zeros([N,T+1])
    
    #initialisation
    for i in range(N):
        ksi[0,i] = rnd.normal(0,sigma_r0)
        ksi[1,i] = rnd.normal(0,sigma_r0)
        ksi[2,i] = rnd.normal(0,sigma_v0)
        ksi[3,i] = rnd.normal(0,sigma_v0)
        
        particules[:,:,0] = np.copy(ksi)
    
    for i in range(N):
        X = np.array([ksi[0,i] , ksi[1,i]])
        w[i] = q_k(hp(0,X) - h_ALT[0])
    w = w/w.sum()
    
    poids[:,0] = np.copy(w)
    
    for k in range(1,T+1):
        resample = resampling_multi(w,N)  #rééchantillonage
        
        for i in range(N):
            ksi_tilde[:,i] = np.copy(ksi[:,resample[i]])   #resélection
            
        for i in range(N):
            
            #calcul des nouveaux ksi
            ksi[: , i] = f(ksi_tilde[: , i] , rnd.normal(0 , sigma_INS , 2))
        for i in range(N):
            #calcul des nouveaux poids
            X = np.array([ksi[0,i] , ksi[1,i]])
            v = hp(k , X) - h_ALT[k]
            w[i] = q_k(v)
        w = w/w.sum()
        
        particules[:,:,k] = np.copy(ksi)
        poids[:,k] = np.copy(w)
        
    for k in range(np.shape(moyenne)[1]):
        moyenne[:,k] = np.dot(particules[0:2,:,k] , poids[:,k])
        
    return(particules, moyenne, poids)


#---------SIS
def filtrage_SIS(N):
    ksi = np.zeros([4,N])
    ksi_tilde = np.zeros([4,N])
    particules = np.zeros([4,N,T+1])
    moyenne = np.zeros([2,T+1])
    
    w = np.zeros([N])
    poids = np.zeros([N,T+1])
    
    #initialisation
    for i in range(N):
        ksi[0,i] = rnd.normal(0,sigma_r0)
        ksi[1,i] = rnd.normal(0,sigma_r0)
        ksi[2,i] = rnd.normal(0,sigma_v0)
        ksi[3,i] = rnd.normal(0,sigma_v0)
        
        particules[:,:,0] = np.copy(ksi)
    
    for i in range(N):
        X = np.array([ksi[0,i] , ksi[1,i]])
        w[i] = q_k(hp(0,X) - h_ALT[0])
    w = w/w.sum()
    
    poids[:,0] = np.copy(w)
    
    for k in range(1,T+1):
            
        for i in range(N):
            #calcul des nouveaux ksi
            ksi[: , i] = f(np.copy(ksi[: , i]) , rnd.normal(0 , sigma_INS , 2))
        for i in range(N):
            #calcul des nouveaux poids
            X = np.array([ksi[0,i] , ksi[1,i]])
            v = hp(k , X) - h_ALT[k]
            w[i] = q_k(v)
        w = w/w.sum()
        particules[:,:,k] = np.copy(ksi)
        poids[:,k] = np.copy(w)
        
    for k in range(np.shape(moyenne)[1]):
        moyenne[:,k] = np.dot(particules[0:2,:,k] , poids[:,k])
        
    return(particules, moyenne, poids)

#---------adaptatif
def filtrage_adaptatif(N,c):
    ksi = np.zeros([4,N])
    ksi_tilde = np.zeros([4,N])
    particules = np.zeros([4,N,T+1])
    moyenne = np.zeros([2,T+1])
    
    w = np.zeros([N])
    poids = np.zeros([N,T+1])
    
    #initialisation
    for i in range(N):
        ksi[0,i] = rnd.normal(0,sigma_r0)
        ksi[1,i] = rnd.normal(0,sigma_r0)
        ksi[2,i] = rnd.normal(0,sigma_v0)
        ksi[3,i] = rnd.normal(0,sigma_v0)
        
        particules[:,:,0] = np.copy(ksi)
    
    for i in range(N):
        X = np.array([ksi[0,i] , ksi[1,i]])
        w[i] = q_k(hp(0,X) - h_ALT[0])
    w = w/w.sum()
    
    poids[:,0] = np.copy(w)
    
    for k in range(1,T+1):
        
        N_eff = 0   #test adaptatif
        for i in range(N):
            N_eff = N_eff + w[i]**2
        N_eff = 1./N_eff
        
        if N_eff < c*N:    #cas de rééchantillonage (SIR)
            resample = resampling_multi(w,N)  #rééchantillonage
            for i in range(N):
                ksi_tilde[:,i] = np.copy(ksi[:,resample[i]])   #resélection
        else:              #cas de non-rééchantillonage (SIS)
            ksi_tilde[:,:] = np.copy(ksi)
            
        for i in range(N):
            #calcul des nouveaux ksi
            ksi[: , i] = f(ksi_tilde[: , i] , rnd.normal(0 , sigma_INS , 2))
        for i in range(N):
            #calcul des nouveaux poids
            X = np.array([ksi[0,i] , ksi[1,i]])
            v = hp(k , X) - h_ALT[k]
            w[i] = q_k(v)
        w = w/w.sum()
        particules[:,:,k] = np.copy(ksi)
        poids[:,k] = np.copy(w)
        
    for k in range(np.shape(moyenne)[1]):
        moyenne[:,k] = np.dot(particules[0:2,:,k] , poids[:,k])
        
    return(particules, moyenne, poids)

def filtrage(N,c,option=2):
    if option==0:
        print("Filtrage SIR")
        (particules, moyenne, poids) = filtrage_SIR(N)
    elif option==1:
        print("Filtrage SIS")
        (particules, moyenne, poids) = filtrage_SIS(N)
    elif option==2:
        print("Filtrage adaptatif")
        (particules, moyenne, poids) = filtrage_adaptatif(N,c)
    else:
        print("Erreur : Option incorrecte")
    fig = plt.figure()
    k = 1
    while k < T+1:
        plt.cla()
        plt.xlim(-10000, 10000)
        plt.ylim(-10000, 10000)
        plt.imshow(map, cmap="jet", extent=[X1MIN, X1MAX, X2MIN, X2MAX])
        plt.plot(particules[0,:,k] + r_INS[0,k] , particules[1,:,k] + r_INS[1,k] ,"w.", label="particules") #trajectoires particulaires
        plt.plot(rtrue[0,:],rtrue[1,:],"r-")  #trajectoire réelle
        plt.plot(r_INS[0,:],r_INS[1,:],"g-")  #mesure intertielle de la trajectoire
        plt.plot(moyenne[0, :] + r_INS[0,:] , moyenne[1, :] + r_INS[1,:] , "k-", markersize=8)
        plt.plot(r_INS[0, k], r_INS[1, k], "g*", label=" r_INS", markersize=8) #trajectoire inertielle
        plt.plot(rtrue[0, k], rtrue[1, k], "r*", label="rtrue", markersize=8)  #vraie trajectoire
        plt.plot(moyenne[0, k] +r_INS[0,k] , moyenne[1, k] + r_INS[1,k] , "k*", markersize=8, label="r moyen") #moyenne de la trajectoire particulaires
    
        plt.legend()
        plt.pause(0.0001)
        k = k+4
    plt.pause(5)
