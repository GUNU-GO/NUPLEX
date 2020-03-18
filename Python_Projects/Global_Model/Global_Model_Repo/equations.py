import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import isclose
from xs import *
from setup import *

#Initialization
def exp_condition(p, input_power, duty, period, time_resolution = 1e-8):
    ng = (p/7.5)/(Tg*kB)*1e-6 #[cm^-3]
    lambda_i = 1/(ng*sigma_i) #[cm] ion-neutral mean free path
    hl = 0.86*(3+l/2/lambda_i)**-0.5
    hR = 0.8*(4+ro/lambda_i)**-0.5

    Aeff = 2*np.pi*ro*(l*hR+ro*hl)  #[cm^2] effective area
    deff = V/Aeff #[cm]
    print('Condition : {}mTorr, {}W, {}ms, {}'.format(p, input_power/6.241509e18, period*1000, duty))

    t_long = None
    iteration_number = None

#Balance equation calculate function
def electron_balance_eqn(density, t, power):
    Te = density[0]
    nH = density[1]
    nHp = density[2]
    nH2p = density[3]
    nH3p = density[4]
    ne = nHp + nH2p + nH3p
    uB = np.sqrt(kB/kB1*Te/M)*100 #[cm/s] #kB/kB1 = e
    uB2 = np.sqrt(kB/kB1*Te/2/M)*100
    uB3 = np.sqrt(kB/kB1*Te/3/M)*100
    #Vs = -Te*np.log(4/ne/np.sqrt(8*e*Te/np.pi/m)*(nHp*uB+nH2p*uB2+nH3p*uB3))
    Vs = Te*np.log(np.sqrt(M/(2*np.pi*m)))
    t0 = V/Aeff*np.sqrt(M/(kB/kB1*Te))/100 #[s] Characteristic transit time of H+ ion 
    #k8,k9,k11의 Te가 매우 작을때의 Cross section값을 구해야한다. (k2는 괜찮음)
    
    ##### Rate coefficient calculation #####
    # ref. R. K. Janev, et al., Elementary Processes in Hydrogen-Helium Plasmas, Springer (1987)
    # H + e -> H+ + 2e, Reaction 2.1.5 E = 13.6
    k1 = np.exp(-3.271396786375e+01+1.353655609057e+01*np.log(Te)-5.739328757388e+00*(np.log(Te))**2+1.563154982022e+00*(np.log(Te))**3-2.877056004391e-01*(np.log(Te))**4+3.482559773737e-02*(np.log(Te))**5-2.631976175590e-03*(np.log(Te))**6+1.119543953861e-04*(np.log(Te))**7-2.039149852002e-06*(np.log(Te))**8)
    # H+ + e -> H + hv, Reaction 2.1.8 E = Te
    k2 = 3.92e-14*(13.6/Te)**1.5/(13.6/Te+0.35) # n = 1s
    # H2 + e -> 2H + e, Reaction 2.2.5 E = 10
    k3 = np.exp(-2.858072836568e+01+1.038543976082e+01*np.log(Te)-5.383825026583e+00*(np.log(Te))**2+1.950636494405e+00*(np.log(Te))**3-5.393666392407e-01*(np.log(Te))**4+1.006916814453e-01*(np.log(Te))**5-1.160758573972e-02*(np.log(Te))**6+7.411623859122e-04*(np.log(Te))**7-2.001369618807e-05*(np.log(Te))**8)
    # H2 + e -> H2+ + 2e, Reaction 2.2.9 E = 15.4
    k4 = np.exp(-3.568640293666e+01+1.733468989961e+01*np.log(Te)-7.767469363538e+00*(np.log(Te))**2+2.211579405415e+00*(np.log(Te))**3-4.169840174384e-01*(np.log(Te))**4+5.088289820867e-02*(np.log(Te))**5-3.832737518325e-03*(np.log(Te))**6+1.612863120371e-04*(np.log(Te))**7-2.893391904431e-06*(np.log(Te))**8)
    # H2 + e -> H+ + H + 2e, Reaction 2.2.10 E = 18
    k5 = np.exp(-3.834597006782e+01+1.426322356722e+01*np.log(Te)-5.826468569506e+00*(np.log(Te))**2+1.727940947913e+00*(np.log(Te))**3-3.598120866343e-01*(np.log(Te))**4+4.822199350494e-02*(np.log(Te))**5-3.909402993006e-03*(np.log(Te))**6+1.738776657690e-04*(np.log(Te))**7-3.252844486351e-06*(np.log(Te))**8)
    # H2+ + e -> 2H+ + 2e, Reaction 2.2.11 E = 15.5
    k6 = np.exp(-3.746192301092e+01+1.559355031108e+01*np.log(Te)-6.693238367093e+00*(np.log(Te))**2+1.981700292134e+00*(np.log(Te))**3-4.044820889297e-01*(np.log(Te))**4+5.352391623039e-02*(np.log(Te))**5-4.317451841436e-03*(np.log(Te))**6+1.918499873454e-04*(np.log(Te))**7-3.591779705419e-06*(np.log(Te))**8)
    # H2+ + e -> H+ + H + e, Reaction 2.2.12 E = 10.5
    k7 = np.exp(-1.781416067709e+01+2.277799785711e+00*np.log(Te)-1.266868411626e+00*(np.log(Te))**2+4.296170447419e-01*(np.log(Te))**3-9.609908013189e-02*(np.log(Te))**4+1.387958040699e-02*(np.log(Te))**5-1.231349039470e-03*(np.log(Te))**6+6.042383126281e-05*(np.log(Te))**7-1.247521040900e-06*(np.log(Te))**8)
    # H2+ + e -> 2H, Reaction 2.2.14 E = Te
    k8 = np.exp(-1.670435653561e+01-6.035644995682e-01*np.log(Te)-1.942745783445e-08*(np.log(Te))**2-2.005952284492e-07*(np.log(Te))**3+2.962996104431e-08*(np.log(Te))**4+2.134293274971e-08*(np.log(Te))**5-6.353973401838e-09*(np.log(Te))**6+6.152557460831e-10*(np.log(Te))**7-2.025361858319e-11*(np.log(Te))**8)
    # H3+ + e -> H2 + H, Reaction 2.2.15 E = Te
    k9 = np.exp(-1.700270758355e+01-4.050073042947e-01*np.log(Te)+1.018733477232e-08*(np.log(Te))**2-1.695586285687e-08*(np.log(Te))**3+1.564311217508e-10*(np.log(Te))**4+1.979725412288e-09*(np.log(Te))**5-4.395545994733e-10*(np.log(Te))**6+3.584926377078e-11*(np.log(Te))**7-1.024189019465e-12*(np.log(Te))**8)
    # H3+ + e -> H+ + 2H + e, Reaction 2.2.16 E = 14
    k10 = np.exp(-3.078408636631e+01+1.509421488513e+01*np.log(Te)-7.349167207324e+00*(np.log(Te))**2+2.320966107642e+00*(np.log(Te))**3-4.818077551719e-01*(np.log(Te))**4+6.389229162737e-02*(np.log(Te))**5-5.161880953089e-03*(np.log(Te))**6+2.303985092606e-04*(np.log(Te))**7-4.344846146197e-06*(np.log(Te))**8)
    # H2+ + H2 -> H3+ + H, Reaction 4.3.3 E = 0
    k11 = 2.1e-9

    if Te < 0.025:
        # H+ + e -> H + hv, Reaction 2.1.8 E = Te
        k2 = 9.137053951846942e-13
        # H2+ + e -> 2H Janev++ 73p
        k8 = 5.156170153467892e-07
        # H3+ + e -> H2 + H Janev++ 98p
        k9 = 1.8393447268390669e-07

    ##### Energy Loss per Reaction #####
    E1 = 13.6
    E2 = Te
    E3 = 10
    E4 = 15.4
    E5 = 18
    E6 = 15.5
    E7 = 10.5
    E8 = Te
    E9 = Te
    E10 = 14
    E11 = 0
    E12 = 12.1
    #Quasi-Neutrality eqn
    ne = nHp + nH2p + nH3p

    #Hydrogen atom conservation eqn
    nH2 = ng - (0.5*(nH+nHp)+nH2p+1.5*nH3p)

    #Particle balance eqn for electron
    dne_dt = (k1*ne*nH)+(k4*ne*nH2)+(k5*ne*nH2)+(k6*ne*nH2p)-(k2*ne*nHp)-(k8*ne*nH2p)-(k9*ne*nH3p)-ne*uB*Aeff/V

    #Power balance eqn for electron
    dTe_dt = 2/(3*ne)*(power(t)/V - (Vs+2.5*Te)*ne*uB*Aeff/V - 3/2*Te*dne_dt\
    - (k1*nH*E1*ne + k2*nHp*E2*ne + k3*nH2*E3*ne + k4*nH2*E4*ne + k5*nH2*E5*ne + k6*nH2p*E6*ne + k7*nH2p*E7*ne\
    + k8*nH2p*E8*ne + k9*nH3p*E9*ne + k10*nH3p*E10*ne + k11*nH2p*E11*nH2))

    #Particle balance eqn for other species except electron
    dnH_dt = (k2*nHp*ne)+2*(k3*nH2*ne)+(k5*nH2*ne)+(k7*nH2p*ne)+2*(k8*nH2p*ne)\
    +(k9*nH3p*ne)+2*(k10*nH3p*ne)+(k11*nH2p*nH2)+(nHp/(t0))+(nH3p/(np.sqrt(3)*t0))-(k1*ne*nH)-(nH*g/T1)

    dnHp_dt = (k1*ne*nH)+(k5*nH2*ne)+2*(k6*nH2p*ne)+(k10*nH3p*ne)-(k2*nHp*ne)-(nHp/(t0))

    dnH2p_dt = (k4*nH2*ne)-(k6*nH2p*ne)-(k7*nH2p*ne)-(k8*nH2p*ne)-(k11*nH2p*nH2)-(nH2p/(np.sqrt(2)*t0))

    dnH3p_dt = (k11*nH2p*nH2)-(k9*nH3p*ne)-(k10*nH3p*ne)-(nH3p/(np.sqrt(3)*t0))

    return [dTe_dt, dnH_dt, dnHp_dt, dnH2p_dt, dnH3p_dt]

#Pulsed power generate function
def power(t):
    if t <= duty*period:
        return input_power
    else:
        return 0

#Temperature & Density Calculation
def calculation(period, time_resolution):
    density = [1.5,1e11,1e11,1e11,1e11] #Te, H, H+, H2+, H3+
    t = np.linspace(0, period, period/time_resolution)
    args = (power,)
    ans1 = odeint(electron_balance_eqn, density, t, args, rtol=10**-3, mxstep=10**6)
    T = ans1[:,0]
    H = ans1[:,1]
    Hp = ans1[:,2]
    H2p = ans1[:,3]
    H3p = ans1[:,4]
    ne = Hp + H2p + H3p
    H2 = ng-(0.5*(H+Hp)+H2p+1.5*H3p)

#Iteration
def iteration():
    iteration_number = 0
    H3p_compare = 1
    Hp_compare = 1
    H2p_compare = 1
    while(not isclose(H2p[-1], H2p_compare, rel_tol=1e-2) or not isclose(H3p[-1], H3p_compare, rel_tol=1e-2) or not isclose(Hp[-1], Hp_compare, rel_tol =1e-2)):
        if iteration_number > 150:
            print('did not converge')
            break
        x0 = [T[-1], H[-1], Hp[-1], H2p[-1], H3p[-1]] #Te, nH, nHp, nH2p, nH3p
        args = (power,)
        H2p_compare = H2p[-1]
        H3p_compare = H3p[-1]
        Hp_compare = Hp[-1]
        ans2 = odeint(electron_balance_eqn, x0, t, args, rtol=10**-2, mxstep=10**3)
        T = np.append(T, ans2[:,0])
        H = np.append(H, ans2[:,1])
        Hp = np.append(Hp, ans2[:,2])
        H2p = np.append(H2p, ans2[:,3])
        H3p = np.append(H3p, ans2[:,4])
        ne = Hp + H2p + H3p
        H2 = ng - (0.5*(H+Hp)+H2p+1.5*H3p)
        iteration_number += 1
    print('iteration count :' + str(iteration_number))
    print('---------------------------------------')
    t_long = np.linspace(0, (iteration_number+1)*period-time_resolution, (iteration_number+1)*int(period/time_resolution))    
    data = np.vstack([H,Hp,H2p,H3p,ne,H2,T])