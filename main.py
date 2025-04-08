import plotting
import transport
import dispersion_relations
import params

import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

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
    t = 200 #ps
    V = -0.0115 # mV
    damping = 4e-4 
    MEC = 0
    ani = 'IP'
    T = 0.3
    type = 'AFM'
    hard_axis = 0
    critical_H = 4.711e6
    Hfield = 1*critical_H
    x_vals = [2000, 2100, 2200, 3000, 3500, 4000]
    # x_vals = [520, 600, 700, 800, 900, 1000]
    # x_vals = [1000, 1200, 1300, 1500, 1700, 2000]

    common_params = {'meshdims' : meshdims, 'cellsize' : cellsize, 't': t, 'V' : V, 'damping' : damping,
              'MEC' : MEC, 'ani' : ani, 'T' : T, 'type' : type, 'hard_axis' : hard_axis, 'Hfield' : Hfield}


    steadystate = params.Steadystate(**common_params, x_vals=x_vals)
    
    timeAvgSA = params.TimeAvgSA(**common_params, x_start=2020, x_stop=4000)
    
    magnonDispersion = params.MagnonDispersion(**common_params, component='y', axis='x', steadystate=True, triple=False)
    
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
    # plotting.plot_magnon_dispersion_with_zoom(magnonDispersion, 1000)
    # plotting.plot_magnon_dispersion_overlay(magnonDispersion, clim_max=10000)
    # transport.current_density(*params)
    # plotting.plot_critical_T(criticalT)
    
    nsm = NSMultiClient(scriptserverports = range(1000,1002), cudaDevices = range(0,2))
    nsm.configure(True, False)


    #### MULTI DISPERSIONS ####

    dispersions = []
    for i in [(0, 0), (1, critical_H)]:
        temp = params.Steadystate(meshdims, cellsize, t, V, damping, 
                                        MEC, ani, T, type, i[0], i[1], False)
        dispersions.append(temp)

    
    nsm.Run(transport.save_steadystate, dispersions)

    # nsm.Run(dispersion_relations.magnon_dispersion, dispersions)


    #### MULTI TIME AVG SA ####

    # meshes = []
    # for i in range(10, 50, 5):
    #     meshes.append([4000, 50, i])

    # Vs = [-0.12, -0.18, -0.25, -0.32, -0.38, -0.45, -0.52, -0.58]

    # timeAvgSAs = []
    # steadyStates = []
    # for i, V in enumerate(Vs):
    #     temp1 = params.Steadystate(meshes[i], cellsize, 200, V, damping,
    #                               MEC, ani, T, type, hard_axis, Hfield)
    #     temp2 = params.TimeAvgSA(meshes[i], cellsize, 100, V, damping,
    #                              MEC, ani, T, type, hard_axis, Hfield, x_start=2020, x_stop=4000)
    #     steadyStates.append(temp1)
    #     timeAvgSAs.append(temp2)

    # nsm.Run(transport.save_steadystate, steadyStates)
    # nsm.Run(transport.time_avg_SA, timeAvgSAs)

if __name__ == '__main__':
    main()