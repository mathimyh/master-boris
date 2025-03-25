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
    t = 100 #ps
    V = -2.9 # mV
    damping = 4e-4 
    MEC = 0
    hard_axis = 0
    ani = 'IP'
    type = 'AFM'
    T = 0.6
    x_vals = [2000, 2100, 2200, 3000, 3500, 4000]
    # x_vals = [520, 600, 700, 800, 900, 1000]
    # x_vals = [1000, 1200, 1300, 1500, 1700, 2000]

    steadystate = params.Steadystate(meshdims, cellsize, t, V, damping, MEC, ani, 
                                                T, type, hard_axis, x_vals)
    
    timeAvgSA = params.TimeAvgSA(meshdims, cellsize, t, V, damping, MEC, ani, T, 
                                                type, hard_axis, 2020, 4000)
    
    magnonDispersion = params.MagnonDispersion(meshdims, cellsize, t, V, damping, 
                                               MEC, ani, T, type, hard_axis, 'y', 'x', False)
    
    criticalT = params.CriticalT(meshdims, cellsize, t, V, damping, 
                                               MEC, ani, T, type, hard_axis,100)

    # ns = NSClient(); ns.configure(True, False)

    # transport.Init_AFM(timeAvgSA)
    # transport.time_avg_SA(meshdims, cellsize,t,V, damping, MEC, ani, T, type, hard_axis, 520, 1000)
    # plotting.plot_tAvg_SA(timeAvgSA)
    # plotting.fft_transport_underneath(*params)
    # transport.save_steadystate2(ns, steadystate)
    # transport.time_avg_SA_underneath(ns, timeAvgSA)
    # plotting.plot_plateau(steadystate)
    # dispersion_relations.critical_T(ns, meshdims, cellsize, t, damping, MEC, ani, type, max_T=100)
    # dispersion_relations.magnon_dispersion(ns, magnonDispersion)
    # dispersion_relations.critical_T(ns, criticalT, 100)
    # plotting.plot_magnon_dispersion_with_zoom(magnonDispersion, 1000)
    # transport.current_density(*params)
    # plotting.plot_critical_T(meshdims, damping, MEC, ani, type)
    
    
    nsm = NSMultiClient(scriptserverports = range(1000,1002), cudaDevices = range(0,2))
    nsm.configure(True, False)
    # meshes = []
    # for i in range(5, 50, 5):
    #     meshes.append([4000, 50, i])

    # Vs = [-0.06, -0.12, -0.18, -0.25, -0.32, -0.38, -0.45, -0.52, -0.58]
    # tot = len(meshes)
    dispersions = []
    for i in [0.6,0.8,1.,3.]:
        temp = params.MagnonDispersion(meshdims, cellsize, t, V, damping, 
                                               MEC, ani, i, type, hard_axis, 'y', 'x', False)
        dispersions.append(temp)

    # nsm.Run(transport.complete_simulation_AFM, meshes, [cellsize]*tot, [500]*tot, [100]*tot, Vs, [damping]*tot, [MEC]*tot, [ani]*tot, [T]*tot, [hard_axis]*tot)
    nsm.Run(dispersion_relations.magnon_dispersion, dispersions)
    # nsm.Run(transport.time_avg_SA_z, meshes, [cellsize]*tot, [t]*tot, Vs, [damping]*tot, [MEC]*tot, [ani]*tot, [T]*tot, [type]*tot)
    # nsm.Run(transport.save_steadystate(ns, [meshdims]*2, [cellsize]*2, [t]*2, [0.002, 0.003], [damping]*2, [MEC]*2, [ani]*2, [T]*2, [type]*2, [x_vals]*2))
    # params = [meshdims, cellsize, t, V, damping, MEC, ani, T, type]
    # nsm.Run(dispersion_relations.magnon_dispersion, meshes, [cellsize]*tot, [t]*tot, Vs, [damping]*tot, [MEC]*tot, [ani]*tot, [T]*tot, [type]*tot, [hard_axis]*tot, ['y']*tot, ['x']*tot, [False]*tot)

    # meshdims, cellsize, t, V, damping, MEC, ani, T, type

if __name__ == '__main__':
    main()