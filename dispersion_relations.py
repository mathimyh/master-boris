import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import plotting
import transport
import params

def magnon_dispersion(ns, magnonDispersion):

    int_dir = 0

    if magnonDispersion.component == 'x':
        int_dir = 1
    elif magnonDispersion.component == 'y':
        int_dir = 2
    elif magnonDispersion.component == 'z':
        int_dir = 3
    else:
        print('Choose direction')
        exit()


    if magnonDispersion.type == 'AFM':
        time_step = 0.1e-12
        total_time = magnonDispersion.t*1e-12
    elif magnonDispersion.type == 'FM':
        time_step = 1e-12
        total_time = magnonDispersion.t*1e-12

    Ms = 2.1e3

    if magnonDispersion.triple:
        # If triple then measure at the edges and in the middle
        y_vals = [magnonDispersion.cellsize,magnonDispersion.meshdims[1]/2,magnonDispersion.meshdims[1]-magnonDispersion.cellsize]
    else:
        # If not triple just find in the middle of the system
        y_vals = [magnonDispersion.meshdims[1]/2]
    
    # Either load a simulation in steadystate or initialize one in ground state
    if magnonDispersion.steadystate:
        sim_name = magnonDispersion.simname()
        ns.loadsim(sim_name)
    
    else:
        if magnonDispersion.type == 'AFM':
            M = transport.Init_AFM(ns, magnonDispersion)
        elif magnonDispersion.type == 'FM':
            M = transport.Init_FM(magnonDispersion)
        else:
            print('Choose type!')
            exit()

    output_files = magnonDispersion.cachename()
    params.make_folder(output_files[0])

    ns.reset()

    time = 0.0
    try:
        ns.cuda(1)
    except:
        print('No cuda!')

    for output_filex in output_files:
        ns.dp_newfile(output_filex)

    while time < total_time:
        
        if magnonDispersion.steadystate:
            ns.V([0.001*magnonDispersion.V, 'time', time + time_step])
        else:
            ns.Relax(['time', time + time_step])
        
        for i, output_filex in enumerate(output_files):
            ns.dp_getexactprofile((np.array([magnonDispersion.cellsize/2, y_vals[i], magnonDispersion.meshdims[2]-magnonDispersion.cellsize])*1e-9), (np.array([magnonDispersion.meshdims[0] - magnonDispersion.cellsize/2, y_vals[i], magnonDispersion.meshdims[2]])*1e-9), magnonDispersion.cellsize*1e-9, 0)
            
            # If we have the biaxial system, then every other measurement is for z-component
            if magnonDispersion.hard_axis and i % 2 != 0:
                ns.dp_div(3, Ms)
                ns.dp_saveappendasrow(output_filex, 3)
            else:
                ns.dp_div(int_dir, Ms)
                ns.dp_saveappendasrow(output_filex, int_dir)
        time += time_step

def critical_T(ns, criticalT):
    
    if criticalT.type == 'AFM':
        M = transport.Init_AFM(ns, criticalT)
        measuring_t = 10e-12
        step_t = 40e-12
    elif criticalT.type == 'FM':
        M = transport.Init_FM(criticalT)
        measuring_t = 100e-12
        step_t = 400e-12
    else:
        print('Choose material type!')
        exit()

    ns.reset()
    ns.iterupdate(200)

    # # Okay so both a temperature increase and a time average over each temperature
    # ns.setdata('<T>', type, np.array([0,0,0,meshdims[0],meshdims[1],meshdims[2]])*1e-9)
    # ns.adddata('<M2>', type, np.array([0,0,0,meshdims[0],meshdims[1],meshdims[2]])*1e-9) # average magnetization of whole mesh (M2 because I want positive values)
    # ns.editdatasave(0, 'time', (t*1e-12)/500) # 10 times for each step, then average over these

    output_file = criticalT.cachename() 
    params.make_folder(output_file)
    data = []

    avs = 30

    # Increase temperature over time 
    for i in range(criticalT.max_T):
        ns.temperature(i)
        ns.Relax(['time', step_t])
        m = [0,0,0]
        for j in range(avs):
            temp = ns.showdata('<M>', type, np.array([0,0,0,criticalT.meshdims[0],criticalT.meshdims[1],criticalT.meshdims[2]])*1e-9)
            m[0] += temp[0]
            m[1] += temp[1]
            m[2] += temp[2]
            ns.Relax(['time', measuring_t])
            ns.reset()
        m[0] /= avs
        m[1] /= avs
        m[2] /= avs
        listy = [i]
        listy.extend(m)
        data.append(listy)

    with open(output_file, 'w') as f:
        for d in data:
            f.write(f"{d}\n")

    plotting.plot_critical_T(criticalT)