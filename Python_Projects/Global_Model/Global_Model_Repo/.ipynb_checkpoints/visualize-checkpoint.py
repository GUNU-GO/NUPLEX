import matplotlib.pyplot as plt
from equations import *
from setup import *

def visualization1():
    plt.figure(figsize=(16,16))
    Power_list = []
    for i in t:
        Power_list.append(power(i)/6.241509e18)

    plt.subplot(421)
    plt.plot(t*1e6,T)
    plt.xlabel('Time (us)')
    plt.ylabel('Temperature (eV)')
    #plt.xlim(440,560)
    #plt.ylim(0,10)
    plt.title('Electron Temperature')
    plt.grid(True)

    plt.subplot(422)
    plt.plot(t*1e6,ne,'brown')
    plt.yscale('log')
    plt.xlabel('Time (us)')
    plt.ylabel('Density (cm-3)')
    plt.title('Electron Density')
    plt.grid(True)
    plt.subplots_adjust(hspace = 0.5)
    plt.savefig(path + 'vis1.png')
    plt.show()

def visualization2(): 
    plt.figure(figsize=(16,10))
    plt.plot(t_long*1e6,H,t_long*1e6,Hp,t_long*1e6,H2,t_long*1e6,H2p,t_long*1e6,H3p,t_long*1e6,ne)
    plt.yscale('log')
    plt.xlabel('Time (us)')
    plt.ylabel('Density (cm-3)')
    #plt.ylim(1e8,1e16)
    plt.legend(['H','H+','H2','H2+','H3+','e'],loc = 'upper right')
    plt.title('Density of All Species')
    plt.grid(True)
    plt.savefig(path + 'vis2.png')
    plt.show()

def visualization3():
    legend_list=['H','H+','H2+','H3+','Electron']
    color_list = ['tab:blue','tab:orange','tab:red','tab:purple','tab:brown']
    start = int(iteration_number*period/time_resolution)
    end = int((iteration_number+1)*period/time_resolution)-1

    plt.figure(figsize=(16,10))
    for i in range(5):
        plt.subplot(4,2,i+1)
        plt.plot(t_long[start:end]*1e6,data[i][start:end],color_list[i])
        plt.yscale('log')
        plt.xlabel('Time (us)')
        plt.ylabel('Density (cm-3)')
        plt.legend([legend_list[i]], loc = 'upper right')
        #plt.ylim(1e8,1e16)
        plt.grid(True)
    plt.savefig(path + 'vis3-1.png')

    plt.figure(figsize=(8,5))
    plt.plot(t_long[start:end]*1e6,H[start:end],'tab:blue',t_long[start:end]*1e6,Hp[start:end],'tab:orange'\
             ,t_long[start:end]*1e6,H2p[start:end],'tab:red',t_long[start:end]*1e6,H3p[start:end],'tab:purple'\
             ,t_long[start:end]*1e6,ne[start:end],'tab:brown',t_long[start:end]*1e6,H2[start:end],'tab:green')
    plt.yscale('log')
    plt.xlabel('Time (us)')
    plt.ylabel('Density (cm-3)')
    plt.title('Density of All Species')
    plt.legend(legend_list, loc = 'upper right')
    plt.grid(True)
    plt.savefig(path + 'vis3-2.png')
    plt.show()
