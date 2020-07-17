#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import quad
from constants import *

Te = np.logspace(-2,3,200)
M = 1.673e-27 #[kg] hydrogen mass
m = 9.109e-31 #[kg] electron mass
k = 1.38064852e-23 #[m2 kg s-2 K-1] Boltzmann constant
e = 1.6e-19 #[C] electron charge
data = pd.read_excel('XS_test.xlsx',header=1)

def maxwellian_with_point_xs(v, Te, reaction_name): 
    return (m/(2*np.pi*e*Te))**1.5*point_xs(v, reaction_name)*v*np.exp(-(m*(v)**2)/(2*e*Te))*4*np.pi*v**2
    #return (m/(2*np.pi*e*Te))**1.5*np.exp(-(m*v**2)/(2*e*Te))*4*np.pi*v**2

def maxwellian_with_analytic_xs(v,Te, reaction_name):
    return (m/(2*np.pi*e*Te))**1.5*analytic_xs(v, reaction_name)*v*np.exp(-(m*(v)**2)/(2*e*Te))*4*np.pi*v**2

def point_xs(v,reaction_name):
    Te = 1/2*m*v**2/e
    Te_data = data[str(reaction_name)+'[eV]']
    xs_data = data[str(reaction_name)+'[cm2]']
    f = interp1d(Te_data, xs_data, fill_value='extrapolate')
    if f(Te) > 0:
        return f(Te)*1e-4
    else:
        return 0
    
def analytic_xs(v,reaction_name):
    Te = 1/2*m*v**2/e
    if 'reaction1_' in reaction_name: # Janev 43p
        dE = 15.42
        if Te < dE:
            return 0
        C_0 = 2.05*dE
        x = Te/dE
        F_0v = np.array([0.092, 0.162, 0.176, 0.155, 0.121])
        xs_ndiss_ion = 1.828/x*(1-1/x**0.92)**2.19*np.log(C_0*x)*1e-16
        i = int(reaction_name[-1])
        return F_0v[i]*xs_ndiss_ion*1e-4
    
    if 'reaction3_' in reaction_name:
        idx = int(reaction_name[-1])-1
        a_list = reaction3_coefficient[idx]
        xs = 0
        for i,a in enumerate(a_list):
            xs += a*np.log(Te)**i
        return np.exp(xs)*1e-4
        
    if 'reaction4_' in reaction_name: #Janev 34p and Shakhatov(2011)
        a_list = np.array([-4.063959689566e+01, 1.636189705461e+01, -3.342841685940e+01, 3.479549344686e+01, -2.082704506646e+01, 7.301916128338e+00, -1.477679988432e+00, 1.596127782326e-01, -7.118499383243e-03])
        xsv0_1 = 0
        for idx,a in enumerate(a_list):
            xsv0_1 += a*np.log(Te)**idx
        v = int(reaction_name[-1])
        return np.exp(xsv0_1)*np.exp(-3*1)/(1+0.05*(v+1))*1e-4
            
    
    if 'reaction5_' in reaction_name: #Janev 47p
        E_th = np.array([3.72,3.21,2.72,2.26,1.83,1.43,1.36,0.713,0.397,0.113,-0.139,-0.354,-0.529,-0.659,-0.736])
        i = int(reaction_name[-1])
        if Te < E_th[i]:
            return 0
        Sigma_v0 = np.array([3.22e-5, 5.18e-4, 4.16e-3 ,2.20e-2 ,1.22e-1 , 4.53e-1, 1.51, 4.48, 10.1, 13.9, 11.8, 8.87, 7.11, 5.00, 3.35])*1e-16
        E_0 = 0.45
        return Sigma_v0[i]*np.exp(-(Te-abs(E_th[i]))/E_0)*1e-4

    if 'reaction6_' in reaction_name:
        idx = int(reaction_name[-1])
        a_list = reaction6_coefficient[idx]
        xs = 0
        for i,a in enumerate(a_list):
            xs += a*np.log(Te)**i
        return np.exp(xs)*1e-4
    
    if 'reaction13' in reaction_name:
        a_list = reaction13_coefficient[0]
        xs = 0
        for i,a in enumerate(a_list):
            xs += a*np.log(Te)**i
        return np.exp(xs)*1e-4

    else:
        return 0
    
    
def rate_constant_with_point_xs(Te, reaction_name):
    return quad(maxwellian_with_point_xs,1,3e8, args=(Te, reaction_name),epsabs=1e-26,limit=100)[0]

def rate_constant_with_analytic_xs(Te, reaction_name):
    return quad(maxwellian_with_analytic_xs,1,3e8, args=(Te, reaction_name),epsabs=1e-36,limit=300)[0]

def quick_rate_constant_with_analytic_xs(Te,reaction_name):
    Te_list = np.array(list(rate_const_dict[reaction_name].keys()))
    rate_const = np.array(list(rate_const_dict[reaction_name].values()))
    if Te < Te_list[0]:
        return 1e-30
    
    elif Te > Te_list[-1]:
        return rate_const[-1]
    
    else:
        f = interp1d(Te_list,rate_const)
        return float(f(Te))

#### reaction_list 생성
reaction_idx = np.array(['reaction1', 'reaction3', 'reaction4', 'reaction5', 'reaction6'])
reaction_list = np.array([])
for reaction in reaction_idx:
    if reaction[-1] =='1':
        for num in range(5):
            reaction_list = np.append(reaction_list,reaction+'_'+str(num))
    if reaction[-1] =='3':
        for num in range(1,7):
            reaction_list = np.append(reaction_list,reaction+'_'+str(num))
    if reaction[-1] =='4':
        for num in range(9):
            reaction_list = np.append(reaction_list,reaction+'_'+str(num))
    if reaction[-1] =='5':
        for num in range(10):
            reaction_list = np.append(reaction_list,reaction+'_'+str(num))
    if reaction[-1] =='6':
        for num in range(10):
            reaction_list = np.append(reaction_list,reaction+'_'+str(num))
reaction_list = np.append(reaction_list,'reaction13')
#### rate_const_list 생성
rate_const_list = list()      
for reaction_name in reaction_list:
    temp_dict = dict()
    rate_const = np.array(list(map(lambda Te:rate_constant_with_analytic_xs(Te,reaction_name),Te))) 
    for key, value in dict(zip(Te,rate_const)).items():
        if value > 1e-28:
            temp_dict[key] = value
    rate_const_list.append(temp_dict)
    
#### rate_const_dict 생성    
rate_const_dict = dict(zip(reaction_list,rate_const_list))