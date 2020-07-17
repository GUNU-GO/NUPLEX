#!/usr/bin/env python
# coding: utf-8

import numpy as np

kB = 1.38e-23 #[J/K] [m2 kg K-1 s-2] Boltzmann constant
e = 1.602e-19 #[C] electronic charge
M = 1.67e-27 #[kg] mass of H atom
m = 9.1e-31 #[kg] mass of electorn
ro = 55e-3 #[m] radius of chamber
l = 240e-3 #[m] chamber length
Tg = 300 #[K] room temperature
sigma_i = 5e-19 #[m2]
rec = 0.1 #Recombination Factor
V = np.pi*ro**2*l #[m^3] discharge volume
A = 2*np.pi*ro*l+2*np.pi*ro**2 #[m^2] loss area
v0 = (8*Tg*kB/(2*M*np.pi))**0.5 #[m/s] mean velocity of H atom
LAMBDAeff = ((2.405/ro)**2+(np.pi/l)**2)**-0.5 #[m]
D_Kn = v0 * LAMBDAeff/3 #[m2/s]
Deff = D_Kn
T1 = LAMBDAeff**2/Deff #[s]

reaction13_coefficient = np.array([-48.92557576,   8.32912175,  -2.18562087])
reaction6_coefficient = np.array([-47.74494377,   4.28536395,  -0.52848558,
                                  -48.49135121,   4.40141876,  -0.54330457, 
                                  -48.58376177,   4.39139036,  -0.54179572, 
                                  -48.58376177,   4.39139036,  -0.54179572, 
                                  -48.67270095,   4.41017493,  -0.54441162,
                                  -48.67270095,   4.41017493,  -0.54441162,
                                  -48.78511339,   4.39766071,  -0.54209366,
                                  -48.94528926,   4.42846597,  -0.54438239, 
                                  -49.18794196,   4.48166038,  -0.55165091, 
                                  -49.13626643,   4.37163242,  -0.53909847])

reaction3_coefficient = np.array([-39.47244639,   2.31509725,  -1.03031797, -44.49719343,
         5.07364316,  -1.8382544 , -46.42448557,   4.95301281,
        -1.8105022 , -52.08688305,   9.59660697,  -3.01319477,
       -55.72982684,  12.00995523,  -3.53911475, -58.08563116,
        13.02636735,  -3.65982417])



reaction13_coefficient = np.reshape(reaction13_coefficient,(int(len(reaction13_coefficient)/3),3))
reaction6_coefficient = np.reshape(reaction6_coefficient,(int(len(reaction6_coefficient)/3),3))
reaction3_coefficient = np.reshape(reaction3_coefficient,(int(len(reaction3_coefficient)/3),3))