import plotting
import transport
import dispersion_relations
import params
import numpy as np

import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

from NetSocks import NSMultiClient   # type: ignore
from NetSocks import NSClient # type:ignore

def main():
    
    # Dimensions (nm)
    Lx = 4000
    Ly = 50
    Lz = 10
    cellsize = 5
    meshdims = (Lx, Ly, Lz)

    # Parameters
    t = 1000 # ps
    V = -0.425 # mV
    damping = 4e-4 
    MEC = 0
    ani = 'IP'
    T = 0.3 # K
    type = 'AFM'
    hard_axis = 0
    critical_H = 4.711e6
    Hfield = 0*critical_H
    x_vals = [2000, 2100, 2200, 3000, 3500, 4000]
    # x_vals = [520, 600, 700, 800, 900, 1000]
    # x_vals = [1000, 1200, 1300, 1500, 1700, 2000]

    common_params = {'meshdims' : meshdims, 'cellsize' : cellsize, 't': t, 'V' : V, 'damping' : damping,
              'MEC' : MEC, 'ani' : ani, 'T' : T, 'type' : type, 'hard_axis' : hard_axis, 'Hfield' : Hfield}


    steadystate = params.Steadystate(**common_params, x_vals=x_vals)
    
    timeAvgSA = params.TimeAvgSA(**common_params, x_start=2020, x_stop=4000)
    
    magnonDispersion = params.MagnonDispersion(**common_params, component='y', axis='x', steadystate=True, triple=True)
    
    magnonDispersionSinc = params.MagnonDispersionSinc(**common_params, component='y', axis='x')
    
    criticalT = params.CriticalT(**common_params, max_T = 100)

    # ns = NSClient(); ns.configure(True, False)

    # transport.Init_AFM(ns, timeAvgSA)
    # transport.time_avg_SA(meshdims, cellsize,t,V, damping, MEC, ani, T, type, hard_axis, 520, 1000)
    # plotting.plot_tAvg_SA(timeAvgSA)
    # plotting.fft_transport_underneath(*params)
    # transport.save_steadystate(ns, steadystate)
    # transport.time_avg_SA(ns, timeAvgSA)
    # plotting.plot_plateau(steadystate)
    # dispersion_relations.critical_T(ns, criticalT)
    # dispersion_relations.magnon_dispersion(ns, magnonDispersion)
    # dispersion_relations.critical_T(ns, criticalT, 100)
    # plotting.plot_magnon_dispersion(magnonDispersion, clim_max=1000)
    # plotting.plot_magnon_dispersion_triple_with_zoom(magnonDispersion, 1000)
    # plotting.plot_magnon_dispersion_overlay(magnonDispersion, clim_max=10000)
    # transport.current_density(*params)
    # plotting.plot_critical_T(criticalT)
    

    ### Parallel computing needs this enabled ###
    nsm = NSMultiClient(scriptserverports = range(1000,1002), cudaDevices = range(0,2))
    nsm.configure(True, False)


    #### MULTI DISPERSIONS ####

    # Jc_dict = {10: -0.12, 15 : -0.18, 20 : -0.25, 25 : -0.32,
    #              30 : -0.38, 35 : -0.45, 40 : -0.52, 45 : -0.58}

    # dispersions = []
    # Vs = [-0.400, -0.425]
    # for v in Vs:
    #     temp = params.MagnonDispersion(meshdims, cellsize, t, v, damping, 
    #                                     MEC, ani, T, type, hard_axis, Hfield, 'y',
    #                                     'x', steadystate=True, triple=True)
    #     dispersions.append(temp)

    # nsm.Run(dispersion_relations.magnon_dispersion, dispersions)

    #### MULTI TIME AVG SA ####

    # meshes = []
    # for i in range(10, 50, 5):
    #     meshes.append([4000, 50, i])

    Vs = []
    now = 0.025
    while now < 0.7:
        Vs.append(round(now, 3))
        now += 0.025

    factor = 0.5

    timeAvgSAs = []
    steadyStates = []
    for v in Vs:
        temp1 = params.Steadystate(meshdims, cellsize, 300, -v*factor, damping,
                                  MEC, ani, T, type, hard_axis, Hfield)
        temp2 = params.TimeAvgSA(meshdims, cellsize, 100, -v*factor, damping,
                                 MEC, ani, T, type, hard_axis, Hfield, x_start=2020, x_stop=4000)
        steadyStates.append(temp1)
        timeAvgSAs.append(temp2)

    nsm.Run(transport.save_steadystate, steadyStates)
    nsm.Run(transport.time_avg_SA, timeAvgSAs)

    hard_axis = 1
    Hfield = 1*critical_H

    timeAvgSAs = []
    steadyStates = []
    for v in Vs:
        temp1 = params.Steadystate(meshdims, cellsize, 300, -v*factor, damping,
                                  MEC, ani, T, type, hard_axis, Hfield)
        temp2 = params.TimeAvgSA(meshdims, cellsize, 100, -v*factor, damping,
                                 MEC, ani, T, type, hard_axis, Hfield, x_start=2020, x_stop=4000)
        steadyStates.append(temp1)
        timeAvgSAs.append(temp2)

    nsm.Run(transport.save_steadystate, steadyStates)
    nsm.Run(transport.time_avg_SA, timeAvgSAs)

if __name__ == '__main__':
    main()