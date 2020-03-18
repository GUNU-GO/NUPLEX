import pandas as pd
import numpy as np

def data_process():
    data = np.vstack([H,Hp,H2,H2p,H3p,ne,T]).T
    data = pd.DataFrame(data, index = t_long*1e6)
    data.columns= ['H','H+','H2','H2+','H3+','ne','Te']
    data.index.name='Time[us]'
    #data.to_csv('Result.csv')
    idx = []
    for i in range(iteration_number+1):
        idx.append(int(period/time_resolution*(i+1)-1))
    Hp_frac = Hp[-1]/ne[-1]
    H2p_frac = H2p[-1]/ne[-1]
    H3p_frac = H3p[-1]/ne[-1]
    frac_data = [Hp_frac, H2p_frac, H3p_frac,T[-1],ne[-1]]
    #print('{},{},{},{},{}'.format(Hp_frac,H2p_frac,H3p_frac,T[-1],ne[-1]*1e-12))

    return data.iloc[idx[-1]], data.iloc[idx], frac_data

