import plotting
import transport
import dispersion_relations
import params
import numpy as np

import sys
import os
sys.path.insert(0, 'C:/Users/matmy/OneDrive/Documents/Boris Data/BorisPythonScripts')

# from NetSocks import NSMultiClient   # type: ignore
# from NetSocks import NSClient # type:ignore

def main():
    
    # Dimensions (nm)
    Lx = 4000
    Ly = 50
    Lz = 75
    cellsize = 5
    meshdims = (Lx, Ly, Lz)

    # Parameters
    t = 1000 # ps
    V = -2.175 # mV
    damping = 4e-4 
    MEC = 0
    ani = 'IP'
    T = 0.3 # K
    type = 'AFM'
    hard_axis = 0
    critical_H = 4.711e6
    Hfield = hard_axis*critical_H
    x_vals = [2000, 2100, 2200, 3000, 3500, 4000]
    # x_vals = [520, 600, 700, 800, 900, 1000]
    # x_vals = [1000, 1200, 1300, 1500, 1700, 2000]

    common_params = {'meshdims' : meshdims, 'cellsize' : cellsize, 't': t, 'V' : V, 'damping' : damping,
              'MEC' : MEC, 'ani' : ani, 'T' : T, 'type' : type, 'hard_axis' : hard_axis, 'Hfield' : Hfield}


    steadystate = params.Steadystate(**common_params, x_vals=x_vals)
    
    timeAvgSA = params.TimeAvgSA(**common_params, x_start=Lx/2 + 20, x_stop=Lx, direction='z')
    
    magnonDispersion = params.MagnonDispersion(**common_params, component='y', axis='x', 
                                               steadystate=0, triple=1)
    
    magnonDispersionSinc = params.MagnonDispersionSinc(**common_params, component='y', axis='x')
    
    criticalT = params.CriticalT(**common_params, max_T = 100)

    # ns = NSClient(); ns.configure(True, False)
    
    # timeavgs = []
    # V0 = 0.06
    # for i in range(1, 10):
    #     temp2 = params.TimeAvgSA([6000, 50, 5*i], cellsize, 100, -V0*i, damping,
    #                              MEC, ani, T, type, hard_axis, Hfield, x_start=3020, x_stop=6000)
    #     timeavgs.append(temp2)
    
    # timeavgs = []
    # factor = Lz/20
    # now = 0.025
    # while now < 1.4:
    #     temp2 = params.TimeAvgSA(meshdims, cellsize, 100, -now*factor, damping,
    #                              MEC, ani, T, type, hard_axis, Hfield, x_start=3020, x_stop=6000)
    #     timeavgs.append(temp2)
    #     now += 0.025

    # # transport.Init_AFM(ns, timeAvgSA)
    # # transport.time_avg_SA(meshdims, cellsize,t,V, damping, MEC, ani, T, type, hard_axis, 520, 1000)
    # for timeavg in timeavgs:
        # plotting.plot_tAvg_SA(timeavg)
    plotting.plot_tAvg_SA_both_systems(timeAvgSA)
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
    # plotting.plot_magnon_dispersion_triple(magnonDispersion, zoom=1, clim_max=1500)
    # plotting.plot_magnon_dispersion_separate(magnonDispersion, 4000)
    # plotting.plot_tAvg_SA(timeAvgSA)
    # plotting.plot_magnon_dispersion(magnonDispersion, zoom=0)

    ### Parallel computing needs this enabled ###
    # nsm = NSMultiClient(scriptserverports = range(1000,1002), cudaDevices = range(0,2))
    # nsm.configure(True, False)

    #### MULTI DISPERSIONS ####

    # dispersions = []
    # for i in [(0, 0), (1, critical_H)]:
    #     temp = params.Steadystate(meshdims, cellsize, t, V, damping, 
    #                                     MEC, ani, T, type, i[0], i[1], False)
    #     dispersions.append(temp)

    
    # nsm.Run(transport.save_steadystate, dispersions)

    # nsm.Run(dispersion_relations.magnon_dispersion, dispersions)

    #### MULTI TIME AVG SA ####

    # meshes = []
    # for i in range(10, 50, 5):
    #     meshes.append([4000, 50, i])

    # Vs = []
    # now = 0.025
    # while now < 0.7:
    #     Vs.append(round(now, 3))
    #     now += 0.025

    # factor = 0.5

    # timeAvgSAs = []
    # steadyStates = []
    # for v in Vs:
    #     temp1 = params.Steadystate(meshdims, cellsize, 300, -v*factor, damping,
    #                               MEC, ani, T, type, hard_axis, Hfield)
    #     temp2 = params.TimeAvgSA(meshdims, cellsize, 100, -v*factor, damping,
    #                              MEC, ani, T, type, hard_axis, Hfield, x_start=2020, x_stop=4000)
    #     steadyStates.append(temp1)
    #     timeAvgSAs.append(temp2)

    # nsm.Run(transport.save_steadystate, steadyStates)
    # nsm.Run(transport.time_avg_SA, timeAvgSAs)

    # hard_axis = 1
    # Hfield = 1*critical_H

    # timeAvgSAs = []
    # steadyStates = []
    # for v in Vs:
    #     temp1 = params.Steadystate(meshdims, cellsize, 300, -v*factor, damping,
    #                               MEC, ani, T, type, hard_axis, Hfield)
    #     temp2 = params.TimeAvgSA(meshdims, cellsize, 100, -v*factor, damping,
    #                              MEC, ani, T, type, hard_axis, Hfield, x_start=2020, x_stop=4000)
    #     steadyStates.append(temp1)
    #     timeAvgSAs.append(temp2)

    # nsm.Run(transport.save_steadystate, steadyStates)
    # nsm.Run(transport.time_avg_SA, timeAvgSAs)

    # factor = 20/20
    # Vs1 = np.linspace(0.200, 0.575, 6)
    # Vs2 = np.linspace(0.175,0.550, 6)
    # Vs = np.concatenate((Vs1, Vs2))

    # for i in range(9):
    #     timeAvgSAs = []
    #     steadyStates = []
    #     for v in Vs:
    #         temp1 = params.Steadystate([6000,50,20], cellsize, 200, -v*factor, damping,
    #                                 MEC, ani, T, type, hard_axis, Hfield)
    #         temp2 = params.TimeAvgSA([6000,50,20], cellsize, 100, -v*factor, damping,
    #                                 MEC, ani, T, type, hard_axis, Hfield, x_start=3020, x_stop=6000)
    #         steadyStates.append(temp1)
    #         timeAvgSAs.append(temp2)
            
    #     for elem in timeAvgSAs:
    #         plotting.plot_tAvg_SA(elem, sim_num=i+1)

if __name__ == '__main__':
    main()