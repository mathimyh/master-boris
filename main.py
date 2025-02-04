import plotting
import transport
import dispersion_relations

import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

from NetSocks import NSMultiClient   # type: ignore

def simulate_2_Ts(ns, steady_params, saving_params):
    transport.save_steadystate(*steady_params)
    transport.time_avg_SA(*saving_params)

def main():
    
    # Dimensions
    Lx = 4000
    Ly = 50
    Lz = 5
    cellsize = 5
    meshdims = (Lx, Ly, Lz)

    # Parameters
    t = 100
    V = -0.03
    data = '<mxdmdt>'
    damping = 4e-4
    MEC = 0
    ani = 'IP'
    type = 'AFM'
    T = 20
    x_vals = [2020, 2300, 2600, 3000, 3500, 4000]
    # x_vals = [520, 600, 700, 800, 900, 1000]

    params = [meshdims, cellsize, t, V, damping, MEC, ani, T, type]

    transport.save_steadystate(*params, x_vals)
    # transport.time_avg_SA(*params, 1520, 3000)
    # # plotting.plot_plateau(meshdims, V, damping, x_vals, MEC, ani, type)
    # dispersion_relations.critical_T(meshdims, cellsize, t, damping, MEC, ani, type, 100)
    # dispersion_relations.critical_T(meshdims, cellsize, t, damping, MEC, ani, type, max_T=500)
    # dispersion_relations.magnon_dispersion(*params, 'y', 'x')
    # plotting.plot_magnon_dispersion(*params, 'y','x', True, 1000)

    # nsm = NSMultiClient(scriptserverports = range(1000,1002), cudaDevices = range(0,2))
    # nsm.configure(True, False)
    # nsm.Run(transport.complete_simulation_AFM, [meshdims]*2, [cellsize]*2, [400]*2, [100]*2, [V]*2, [damping]*2, [MEC]*2, [ani]*2, [10,20])
    # plotting.plot_tAvg_SA(*params, 1520, 3000)

if __name__ == '__main__':
    main()