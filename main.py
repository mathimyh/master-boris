import plotting
import transport
import dispersion_relations

import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

from NetSocks import NSMultiClient   # type: ignore
from NetSocks import NSClient # type:ignore

def main():
    
    # Dimensions
    Lx = 4000
    Ly = 50
    Lz = 5
    cellsize = 5
    meshdims = (Lx, Ly, Lz)

    # Parameters
    t = 100
    V = -0.06
    damping = 4e-4
    MEC = 0
    hard_axis = 1
    ani = 'IP'
    type = 'AFM'
    T = 0.3
    x_vals = [2000, 2100, 2200, 3000, 3500, 4000]
    # x_vals = [520, 600, 700, 800, 900, 1000]
    # x_vals = [1000, 1200, 1300, 1500, 1700, 2000]

    params = [meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis]
    ns = NSClient(); ns.configure(True, False)

    # transport.time_avg_SA(*params, 2020, 4000)
    plotting.plot_tAvg_SA(*params, 2020, 4000)
    # plotting.fft_transport_underneath(*params)
    # transport.save_steadystate(ns, *params, x_vals)
    # transport.time_avg_SA_underneath(ns, *params, 2020,4000)
    # plotting.plot_tAvg_SA(*params, 2020,4000)
    # transport.Init_FM(meshdims, cellsize, damping, MEC, ani, T)
    # transport.time_avg_SA(*params, 2020, 4000)
    # plotting.plot_plateau(meshdims, V, damping, x_vals, MEC, ani, T, type)
    # dispersion_relations.critical_T(ns, meshdims, cellsize, t, damping, MEC, ani, type, max_T=100)
    # dispersion_relations.magnon_dispersion(ns, *params, 'y', 'x', False)
    # plotting.plot_magnon_dispersion(*params, 'y','x', False, 1000)
    # transport.current_density(*params)
    # plotting.plot_critical_T(meshdims, damping, MEC, ani, type)
    
    # nsm = NSMultiClient(scriptserverports = range(1000,1002), cudaDevices = range(0,2))
    # nsm.configure(True, False)
    # meshes = []
    # for i in range(5, 50, 5):
    #     meshes.append([4000, 50, i])

    # Vs = [-0.06, -0.12, -0.18, -0.25, -0.32, -0.38, -0.45, -0.52, -0.58]
    # tot = len(meshes)

    # nsm.Run(transport.complete_simulation_AFM_z, meshes, [cellsize]*tot, [1e6]*tot, [0]*tot, Vs, [damping]*tot, [MEC]*tot, [ani]*tot, [T]*tot)
    # nsm.Run(dispersion_relations.magnon_dispersion, meshes, [cellsize]*tot, [t]*tot, Vs, [damping]*tot, [MEC]*tot, [ani]*tot, [T]*tot, [type]*tot, ['y']*tot, ['x']*tot, [True]*tot)
    # nsm.Run(transport.time_avg_SA_z, meshes, [cellsize]*tot, [t]*tot, Vs, [damping]*tot, [MEC]*tot, [ani]*tot, [T]*tot, [type]*tot)
    # nsm.Run(transport.save_steadystate(ns, [meshdims]*2, [cellsize]*2, [t]*2, [0.002, 0.003], [damping]*2, [MEC]*2, [ani]*2, [T]*2, [type]*2, [x_vals]*2))
    # params = [meshdims, cellsi4ze, t, V, damping, MEC, ani, T, type]

    meshdims, cellsize, t, V, damping, MEC, ani, T, type

if __name__ == '__main__':
    main()