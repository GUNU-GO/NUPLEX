#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# data는 [[Te],[nH],[nHm], .... ] 형식으로 들어옴
# [Te, ne, nH, nH_2s, nH2_v0, nH2_v1, nH2_v2, nH2_v3, nH2_v4, nH2_v5, nH2_v6, nH2_v7, nH2_v8, nH2_v9, nHp, nH2p, nH3p, nHm]
# exp_condition은 [p, input_power, duty, period, time_resolution] 과 같은 list 형식으로 들어올거임
path = 'model_result/'

def plot(data, input_species, add_path, file_save = False, xlim = (), ylim = ()): #species 는 'Te,ne,nH' 이런식으로 들어올거임
    file_name = ''
    input_species = input_species.replace(' ','').split(',')
    p, input_power, duty, period, time_resolution = data[1]
    time_interval = np.linspace(0, period, int(period/time_resolution))
    plt.figure(figsize=(14,7))
    
    
    if 'alpha' in input_species:
        ylabel = 'Electronegativity'
        plt.plot(time_interval*1e6,data[0]['nHm']/data[0]['ne'])
        plt.xlabel('Time(μs)')
        plt.ylabel(ylabel)
        if xlim == () and ylim == ():
            pass
        else:
            plt.xlim(xlim)
            plt.ylim(ylim)
        plt.grid(True)
        file_name += 'alpha' + ','
        
    for species in input_species:
        if species == 'alpha':
            break
        
        if species == 'Te':
            plt.plot(time_interval*1e6,data[0][species])
            plt.xlabel('Time(μs)')
            plt.ylabel('Temperature$(eV)$')
            if xlim == () and ylim == ():
                pass
            else:
                plt.xlim(xlim)
                plt.ylim(ylim)
            plt.grid(True)
            file_name += species + ',' 
            break
            
        else:
            plt.plot(time_interval*1e6,data[0][species])
            plt.yscale('log')
            plt.xlabel('Time(μs)')
            plt.ylabel('Density$(m^{-3})$')
            if xlim == () and ylim == ():
                pass
            else:
                plt.xlim(xlim)
                plt.ylim(ylim)
            plt.grid(True)
            file_name += species + ',' 
            
    plt.legend(input_species)
    if file_save:
        plt.savefig(path + add_path + file_name+'(p-{}mTorr,power-{}W,duty-{},period-{}s)'.format(p, input_power/6.241509e18, duty, period) + '.png')
    plt.show()          

def result_to_csv(data,add_path):
    df = pd.DataFrame.from_dict(data[0])
    p, input_power, duty, period, time_resolution = data[1]
    df.to_csv(path + add_path +'(p-{}mTorr,power-{}W,duty-{},period-{}s)'.format(p, input_power/6.241509e18, duty, period) + '.csv')
    
def result_to_csv_select(data,select_list,add_path):
    df = pd.DataFrame.from_dict(data[0])
    file_name = select_list
    p, input_power, duty, period, time_resolution = data[1]
    select_list = select_list.replace(' ','').split(',')
    df[select_list].to_csv(path + add_path + str(file_name) +'(p-{}mTorr,power-{}W,duty-{},period-{}s)'.format(p, input_power/6.241509e18, duty, period) + '.csv')
    
def result_to_csv_select_quick(data,select_list,add_path):
    df = pd.DataFrame.from_dict(data[0])
    file_name = select_list
    p, input_power, duty, period, time_resolution = data[1]
    select_list = select_list.replace(' ','').split(',')
    df[select_list].loc[0:1].to_csv(path + add_path + 'quick' +str(file_name) +'(p-{}mTorr,power-{}W,duty-{},period-{}s)'.format(p, input_power/6.241509e18, duty, period) + '.csv')
    
def collector(data_list,select_list,add_path):
    file_name = select_list
    for data in data_list:
        df = pd.DataFrame.from_dict(data[0])
        p, input_power, duty, period, time_resolution = data[1]
        
        select_list = select_list.replace(' ','').split(',')
        df[select_list].to_csv(path + add_path + 'collect'+str(file_name) +'(p-{}mTorr,power-{}W,duty-{},period-{}s)'.format(p, input_power/6.241509e18, duty, period) + '.csv')
    