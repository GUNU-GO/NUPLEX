#!/usr/bin/env python
#coding=utf8

from numpy import sin, pi, sqrt, arccos, log
from pandas import read_excel

e = 1.602e-19 # [C] electron charge
r_p = 0.15e-3 # [m] probe radius
l_p = 1e-3 # [m] probe length
h = 0.5e-3 # [m] Hole radius
s = 0.7e-3 # [m] Rotation center to Hole edge
R = 0.6e-3 # [m] Rotation center to Wire center
m_i = (19*2+10) * 1.67e-27 #[kg] mass of BF2+
k = 1.38e-23  #[m2kg/s2K] Boltzmann const
alpha = pi/2 # [rad] angle between B-field and Rotation center to Wire center
gamma = (1+0.5)/(2+0.5)

class Machprobe():
    def __init__(self, ne, Te, Ti, m_i, I):
        self.ne = ne
        self.I = I
        self.Cs =sqrt(e*(Te+Ti)/(m_i))
        d_alpha = arccos((s**2 + R**2 - h**2)/(2*s*R))
        self.A_eff = l_p*(R*sin(alpha)+r_p-max(R*sin(alpha)-r_p, s*sin(alpha-d_alpha)))
        #print('Te : {} eV'.format(Te))
        #print('Ti : {} eV'.format(Ti))
        #print('Effective area : {} m2'.format(self.A_eff))
        #print('Ion sound speed : {} m/s'.format(self.Cs))
        
    def perp_current(self):
        self.I_D = (r_p/l_p)*(1-gamma)*self.A_eff     # diffusion current calculation
        self.I_sat = gamma*e*self.A_eff*self.Cs*self.ne     # saturation current calculation
        self.I_perp = self.I - self.I_D - self.I_sat             # perpendicular current calculation
        
        #print('diffusion current : ',self.I_D)
        #print('saturation current : ',self.I_sat)
        #print('perp current : ',self.I_perp)

file_path = 'Mach_Probe_datasheet.xlsx'
data = read_excel(file_path, encoding='utf8')

for i in range(len(data)):
    I_upstream = data.loc[i,'I_Upstream [A]']
    I_downstream = data.loc[i,'I_Downstream [A]']
    ne = data.loc[i,'Electron density [m-3]']
    Te = data.loc[i,'Electron temperature [eV]']
    Ti = data.loc[i,'Ion temperature [eV]']

    upstream = Machprobe(ne, Te, Ti, m_i, I_upstream)
    downstream = Machprobe(ne, Te, Ti, m_i, I_downstream)
    
    upstream.perp_current()
    downstream.perp_current()
    
    data.loc[i,'Mach number'] = 0.73*log(upstream.I_perp/downstream.I_perp)
    data.loc[i,'I_sat (Up) [A]'] = upstream.I_sat
    data.loc[i,'I_D (Up) [A]'] = upstream.I_D
    data.loc[i,'I_Perp (Up) [A]'] = upstream.I_perp
    data.loc[i,'I_sat (Down) [A]'] = downstream.I_sat
    data.loc[i,'I_D (Down) [A]'] = downstream.I_D
    data.loc[i,'I_Perp (Down) [A]'] = downstream.I_perp
    
data.to_excel('Result_' + file_path, encoding='utf8')





