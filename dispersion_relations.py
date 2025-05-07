import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

from NetSocks import NSClient   # type: ignore
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import plotting
import transport
import params

path = 'C:/Users/mathimyh/master/master-boris/'

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

    y_vals = [5, 5, 5]
    
    
    if magnonDispersion.steadystate:
        sim_name = magnonDispersion.simname()
        ns.loadsim(sim_name)
        y_vals = [5,25,45]
    
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
    ns.cuda(1)

    for output_filex in output_files:
        ns.dp_newfile(output_filex)

    while time < total_time:
        
        if magnonDispersion.steadystate:
            ns.V([0.001*magnonDispersion.V, 'time', time + time_step])
        else:
            ns.Relax(['time', time + time_step])
        
        for i, output_filex in enumerate(output_files):
            # if axis == 'x':
            ns.dp_getexactprofile((np.array([magnonDispersion.cellsize/2, y_vals[i], magnonDispersion.meshdims[2]-magnonDispersion.cellsize])*1e-9), (np.array([magnonDispersion.meshdims[0] - magnonDispersion.cellsize/2, y_vals[i], magnonDispersion.meshdims[2]])*1e-9), magnonDispersion.cellsize*1e-9, 0)
            if magnonDispersion.hard_axis and i % 2 != 0:
                ns.dp_div(3, Ms)
                ns.dp_saveappendasrow(output_filex, 3)
            else:
                ns.dp_div(int_dir, Ms)
                ns.dp_saveappendasrow(output_filex, int_dir)
        time += time_step

    # plotting.plot_magnon_dispersion(magnonDispersion)

def phonon_dispersion(meshdims, cellsize, t, damping, x_start, x_stop, MEC, ani, dir):

    dir1 = 0

    if dir == 'x':
        dir1 = 1
    elif dir == 'y':
        dir1 = 2
    elif dir == 'z':
        dir1 = 3
    else:
        print('Choose direction')
        exit()

    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/cache/' + mec_folder + 'dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    time_step = 0.01e-12 
    total_time = t*1e-12

    # sim_name = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/sims/V' + str(V) + '_damping' + str(damping) + '_steady_state.bsm'
    sim_name = 'C:/Users/mathimyh/Documents/Boris Data/Simulations/boris_fordypningsoppgave/' + ani + '/sims/' + mec_folder + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/ground_state.bsm'

    ns = NSClient(); ns.configure(True, False)
    
    ns.loadsim(sim_name)
    ns.reset()

    ns.editstagestop(0, 'time', 10e-12)
    ns.Run()
    ns.reset()

    time = 0.0

    output_file = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 'dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_phonon_dispersion.txt'
    ns.dp_newfile(output_file)

    while time < total_time:
        ns.editstagestop(0, 'time', time + time_step)
        ns.Run()
        ns.dp_getexactprofile('u', [x_start * 1e-9 + cellsize*1e-9/2, 50e-9/2 + cellsize*1e-9/2, 0], [x_stop * 1e-9 - cellsize*1e-9/2, 50e-9/2 + cellsize*1e-9/2, 0], cellsize*1e-9, 0)
        ns.dp_div(dir1, 1e-13)
        ns.dp_saveappendasrow(output_file, dir1)
        time += time_step

    plotting.plot_phonon_dispersion(meshdims, damping, MEC, ani, dir, time_step)

def trajectory(meshdims, t, damping, x_start, x_stop, MEC, ani, dir):

    dir1 = 0

    if dir == 'x':
        dir1 = 1
    elif dir == 'y':
        dir1 = 2
    elif dir == 'z':
        dir1 = 3
    else:
        print('Choose direction')
        exit()

    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/cache/' + mec_folder + 'trajectory/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    time_step = 0.1e-12

    Ms = 2.1e3

    # sim_name = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/sims/V' + str(V) + '_damping' + str(damping) + '_steady_state.bsm'
    sim_name = 'C:/Users/mathimyh/Documents/Boris Data/Simulations/boris_fordypningsoppgave/' + ani + '/sims/' + mec_folder + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/ground_state.bsm'

    ns = NSClient(); ns.configure(True, False)
    
    ns.loadsim(sim_name)
    ns.reset()

    time = 0.0
    ns.cuda(1)

    output_file1 = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 'trajectory/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + dir +  '_trajectory_M1.txt'
    output_file2 = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 'trajectory/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + dir +  '_trajectory_M2.txt'

    ns.dp_newfile(output_file1)
    ns.dp_newfile(output_file2)

    while time < t*1e-12:
        ns.editstagestop(0, 'time', time + time_step)
        ns.Run()
        ns.dp_getexactprofile('M', [x_start * 1e-9 + 5e-9/2, meshdims[1]*1e-9/2 + 5e-9/2, 0], [x_stop * 1e-9 - 5e-9/2, meshdims[1]*1e-9/2 + 5e-9/2, 0], 5e-9, 0)
        ns.dp_getexactprofile('M2', [x_start * 1e-9 + 5e-9/2, meshdims[1]*1e-9/2 + 5e-9/2, 0], [x_stop * 1e-9 - 5e-9/2, meshdims[1]*1e-9/2 + 5e-9/2, 0], 5e-9, 4)
        ns.dp_div(dir1, Ms)
        ns.dp_div(dir1+4, Ms)
        ns.dp_saveappendasrow(output_file1, dir1)
        ns.dp_saveappendasrow(output_file2, dir1+4)
        time += time_step

    plotting.plot_trajectory(meshdims, damping, MEC, ani, dir)

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

def magnon_dispersion_sinc(ns, magnonDispersionSinc):

    int_dir = 0

    if magnonDispersionSinc.component == 'x':
        int_dir = 1
    elif magnonDispersionSinc.component == 'y':
        int_dir = 2
    elif magnonDispersionSinc.component == 'z':
        int_dir = 3
    else:
        print('Choose direction')
        exit()

    ns.cuda(1)    # ns.selectcudadevice([0,1])
    ns.reset()
    
    modules = ['exchange', 'aniuni', 'Zeeman']
    if magnonDispersionSinc.MEC:
        modules.append('melastic')

    # Set up the antiferromagnet
    AFM = ns.AntiFerromagnet(np.array(magnonDispersionSinc.meshdims)*1e-9, [magnonDispersionSinc.cellsize*1e-9])
    AFM.modules(modules)
    ns.setode('LLG', 'RK4') # No temperature
    ns.setdt(1e-15)

    # Set parameters. Should possibly just save these in the database really    
    AFM.param.grel_AFM = 1
    AFM.param.damping_AFM = magnonDispersionSinc.damping
    AFM.param.Ms_AFM = 2.1e3
    AFM.param.Nxy = 0
    AFM.param.A_AFM = 76e-15 # J/m
    AFM.param.Ah = -460e3 # J/m^3
    AFM.param.Anh = 0.0
    AFM.param.J1 = 0
    AFM.param.J2 = 0
    if magnonDispersionSinc.hard_axis:
        AFM.param.K1_AFM = -21e-3 # J/m^3
        AFM.param.K2_AFM = 21 # J/m^3
    else:
        AFM.param.K1_AFM = 21 # J/m^3
        AFM.param.K2_AFM = 0
    AFM.param.K3_AFM = 0
    AFM.param.cHa = 1
    # Different anisotropies
    if magnonDispersionSinc.ani == 'OOP':
        AFM.param.ea1 = (0,0,1) # Set it z-direction
        AFM.setangle(0,90) # Just move m to z-direction, not necessary to wait every time
    elif magnonDispersionSinc.ani == 'IP':
        AFM.param.ea1 = (1,0,0)
        AFM.setangle(90,0)
    else:
        print('Choose anisotropy direction')
        exit()


    Ms = 2.1e3

    output_files = magnonDispersionSinc.cachename()
    params.make_folder(output_files[0])
    ns.reset()

    time = 0.0
    ns.cuda(1)

    for output_file in output_files:
        ns.dp_newfile(output_file)

    AFM.pbc('x', 10)

    equiltime = 0
    total_time = 2 *magnonDispersionSinc.t*1e-12
    H0 = 0
    He = 500e3

    ns.setstage('Hequation')
    ns.editstagevalue(0, 'H0, He *sinc(kc*(x-Lx/2))*sinc(kc*(y-Ly/2))*sinc(2*PI*fc*(t-t0)),He *sinc(kc*(x-Lx/2))*sinc(kc*(y-Ly/2))*sinc(2*PI*fc*(t-t0))')

    N = 600
    L = magnonDispersionSinc.meshdims[0]
    kc = 2*np.pi*N/(2*L)
    fc = 5e12
    time_step = 0.1e-12

    ns.equationconstants('H0', H0)
    ns.equationconstants('He', He)
    ns.equationconstants('kc', kc)
    ns.equationconstants('fc', fc)
    ns.equationconstants('t0', magnonDispersionSinc.t)

    #setup data extraction in command buffer : need to extract a normalized magnetization profile at each time_step and append it to output file
    ns.setdata('commbuf')
    ns.editstagestop(0, 'time', total_time)
    ns.editdatasave(0, 'time', time_step)

    ns.clearcommbuffer()
    #the following 3 commands are not executed now, but buffered so they can be executed at runtime every time_step ('commbuf' output data set)
    #get magnetisation profile along length through center


    #get magnetisation profile along length through center
    ns.dp_getexactprofile((np.array([magnonDispersionSinc.cellsize/2, magnonDispersionSinc.meshdims[1]/2, magnonDispersionSinc.meshdims[2]-magnonDispersionSinc.cellsize])*1e-9), 
                        (np.array([magnonDispersionSinc.meshdims[0] - magnonDispersionSinc.cellsize/2, magnonDispersionSinc.meshdims[1]/2, magnonDispersionSinc.meshdims[2]])*1e-9), 
                            magnonDispersionSinc.cellsize*1e-9, 0)
    
    #save only the y component of magnetisation at time_ste
    ns.dp_div(2, Ms, bufferCommand=True)
    ns.dp_saveappendasrow(output_files[0], 2, bufferCommand = True)

    ns.Run()

    # plotting.plot_magnon_dispersion_with_zoom(magnonDispersion)