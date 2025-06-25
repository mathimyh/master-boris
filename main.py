import plotting
import transport
import dispersion_relations
import params
import numpy as np

import sys
import os
sys.path.insert(0, 'C:/Users/matmy/OneDrive/Documents/Boris Data/BorisPythonScripts')

from NetSocks import NSMultiClient   # type: ignore
from NetSocks import NSClient # type:ignore

def main():
    
    # Dimensions (nm)
    Lx = 4000
    Ly = 50
    Lz = 5
    cellsize = 5
    meshdims = (Lx, Ly, Lz)

    # Parameters
    t = 100 # ps
    V = -0.01 # mV
    damping = 4e-4 
    MEC = 0
    ani = 'IP'
    T = 0.3 # K
    type = 'AFM'
    hard_axis = 1
    critical_H = 4.711e6
    Hfield = 1*hard_axis*critical_H
    x_vals = [2000, 2100, 2200, 3000, 3500, 4000]

    common_params = {'meshdims' : meshdims, 'cellsize' : cellsize, 't': t, 'V' : V, 'damping' : damping,
              'MEC' : MEC, 'ani' : ani, 'T' : T, 'type' : type, 'hard_axis' : hard_axis, 'Hfield' : Hfield}


    steadystate = params.Steadystate(**common_params, x_vals=x_vals)
    
    timeAvgSA = params.TimeAvgSA(**common_params, x_start=Lx/2 + 20, x_stop=Lx, direction='x')
    
    magnonDispersion = params.MagnonDispersion(**common_params, component='y', axis='x', 
                                               steadystate=0, triple=0)

    criticalT = params.CriticalT(**common_params, max_T = 100)
        
        
    # try:
    #     ns = NSClient(); ns.configure(True, False)
    # except:
    #     try:
    #         nsm = NSMultiClient(scriptserverports = range(1000,1002), cudaDevices = range(0,2))
    #     except:
    #         print('No Boris open!') 
        
        
    plotting.plot_magnon_dispersion(magnonDispersion, zoom=1, analytical=1, clim_max=5000)   
        
    # for i in range(1, 9):
    #     temp2 = params.TimeAvgSA([6000, 5, 5*i], cellsize, 100, V*i, damping,
    #                              MEC, ani, T, type, hard_axis, Hfield, x_start=3020, x_stop=6000)
    #     temp3 = params.Steadystate([6000, 5, 5*i], cellsize, 150, V*i, damping, MEC,
    #                                ani, T, type, hard_axis, Hfield)
        # transport.save_steadystate(ns, temp3)
        # transport.time_avg_SA(ns, temp2)
        

if __name__ == '__main__':
    main()