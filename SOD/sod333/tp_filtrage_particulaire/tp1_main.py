#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 15:28:02 2020

@author: ilyes ER-RAMMACH / farès EL BOUZIANI
"""


from tp1_fn import *

#--------------affichage carte
#plt.imshow(map , cmap='jet', extent=[X1MIN,X1MAX,X2MIN,X2MAX])

#------------QUESTION 1
#plt.plot(rtrue[0,:],rtrue[1,:],"r-")

#------------QUESTION 2
#plt.plot(r_INS[0,:],r_INS[1,:],"m-")

#------------QUESTIONS 4 et 5
#plt.plot(np.linspace(0,100,101) , H)   #relief réel
#plt.plot(np.linspace(0,100,101) , h_ALT , color = 'red' , marker = 'x' , linestyle='none')    #relief mesuré

#------------QUESTIONS 7 et 8
#   mode d'emploi : 
# - choisir N, le nombre de particules
# - choisir l'option (0 pour SIR, 1 pour SIS, et 2 pour adaptatif)
# - choisir c, le taux de rééchantillonage (NECESSAIRE POUR FILTRAGE ADAPTATIF)
N = 500
c = 0.75
option = 2
filtrage(N,c,option)