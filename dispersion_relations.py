import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

from NetSocks import NSClient   # type: ignore
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import plotting
import transport

path = 'C:/Users/mathimyh/master/master-boris/'

def magnon_dispersion(ns, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, dir, axis, steadystate = False):

    # # ns = NSClient(); ns.configure(True, False)

    int_dir = 0

    if dir == 'x':
        int_dir = 1
    elif dir == 'y':
        int_dir = 2
    elif dir == 'z':
        int_dir = 3
    else:
        print('Choose direction')
        exit()

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += 'mec'
    if hard_axis:
        modules_folder += '+hard_axis'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    if type == 'AFM':
        time_step = 0.1e-12
        total_time = t*1e-12
    elif type == 'FM':
        time_step = 1e-12
        total_time = t*1e-12

    Ms = 2.1e3

    y_vals = [25, 5, 45]
    
    if steadystate:
        sim_name = path + type + '/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K_steady_state.bsm'
        ns.loadsim(sim_name)
        output_file1 = path + type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_axis' + axis + 'V' + str(V) + '_damping' + str(damping) + '_T' + str(T) + '_y=' + str(y_vals[0]) + '_dispersion.txt'
        output_file2 = path + type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_axis' + axis + 'V' + str(V) + '_damping' + str(damping) + '_T' + str(T) + '_y=' + str(y_vals[1]) + '_dispersion.txt'
        output_file3 = path + type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_axis' + axis + 'V' + str(V) + '_damping' + str(damping) + '_T' + str(T) + '_y=' + str(y_vals[2]) + '_dispersion.txt'
        output_files = [output_file1, output_file2, output_file3]

    else:
        if type == 'AFM':
            M = transport.Init_AFM(meshdims, cellsize, damping, MEC, ani, T, hard_axis)
        elif type == 'FM':
            M = transport.Init_FM(meshdims, cellsize, damping, MEC, ani, T, hard_axis)
        else:
            print('Choose type!')
            exit()
        output_file = path + type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_axis' + axis + 'groundstate' + '_damping' + str(damping) + '_T' + str(T) +  '_dispersion.txt'
        output_files = [output_file]

    ns.reset()

    time = 0.0
    ns.cuda(1)

    for output_filex in output_files:
        ns.dp_newfile(output_filex)

    while time < total_time:
        
        if steadystate:
            ns.V([0.001*V, 'time', time + time_step])
        else:
            ns.Relax(['time', time + time_step])
        
        for i, output_filex in enumerate(output_files):
            # if axis == 'x':
            ns.dp_getexactprofile((np.array([cellsize/2, y_vals[i], meshdims[2]-cellsize])*1e-9), (np.array([meshdims[0] - cellsize/2, y_vals[i], meshdims[2]])*1e-9), cellsize*1e-9, 0)
            # elif axis == 'z':
            #     ns.dp_getexactprofile((np.array([meshdims[0]/2-cellsize/2, meshdims[1]/2-cellsize/2, meshdims[2]-cellsize/2])*1e-9), (np.array([meshdims[0]/2-cellsize/2, meshdims[1]/2-cellsize/2, 40+cellsize/2])*1e-9), cellsize*1e-9, 0)
            ns.dp_div(int_dir, Ms)
            ns.dp_saveappendasrow(output_filex, int_dir)
        time += time_step

    plotting.plot_magnon_dispersion(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, dir, axis, steadystate)

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

def critical_T(ns, meshdims, cellsize, t, damping, MEC, ani, type, max_T):
    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += 'mec'
    modules_folder += '/'
    
    folder_name = type + '/' + modules_folder + ani + '/cache/critical_T/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    if type == 'AFM':
        M = transport.Init_AFM(meshdims, cellsize, damping, MEC, ani, 0)
        measuring_t = 10e-12
        step_t = 40e-12
    elif type == 'FM':
        M = transport.Init_FM(meshdims, cellsize, damping, MEC, ani, 0)
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

    output_file = path + type + '/' + modules_folder + ani + '/cache/critical_T/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/critical_T.txt'

    data = []

    avs = 20

    # Increase temperature over time 
    for i in range(max_T):
        ns.temperature(i)
        ns.Relax(['time', step_t])
        m = [0,0,0]
        for j in range(avs):
            temp = ns.showdata('<M>', type, np.array([0,0,0,meshdims[0],meshdims[1],meshdims[2]])*1e-9)
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

    plotting.plot_critical_T(meshdims, damping, MEC, ani, type)

def magnon_dispersion_sinc(ns, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, dir, axis):
    int_dir = 0

    if dir == 'x':
        int_dir = 1
    elif dir == 'y':
        int_dir = 2
    elif dir == 'z':
        int_dir = 3
    else:
        print('Choose direction')
        exit()

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += 'mec'
    if hard_axis:
        modules_folder += '+hard_axis'
    modules_folder += '/'

    ns.cuda(1)    # ns.selectcudadevice([0,1])
    ns.reset()
    
    modules = ['exchange', 'aniuni']
    if MEC:
        modules.append('melastic')

    # Set up the antiferromagnet
    AFM = ns.AntiFerromagnet(np.array(meshdims)*1e-9, [cellsize*1e-9])
    AFM.modules(modules)
    ns.setode('LLG', 'RK4') # No temperature
    ns.setdt(1e-15)

    # Set parameters. Should possibly just save these in the database really    
    AFM.param.grel_AFM = 1
    AFM.param.damping_AFM =  damping
    AFM.param.Ms_AFM = 2.1e3
    AFM.param.Nxy = 0
    AFM.param.A_AFM = 76e-15 # J/m
    AFM.param.Ah = -460e3 # J/m^3
    AFM.param.Anh = 0.0
    AFM.param.J1 = 0
    AFM.param.J2 = 0
    if hard_axis:
        AFM.param.K1_AFM = -21e-3 # J/m^3
        AFM.param.K2_AFM = 21 # J/m^3
    else:
        AFM.param.K1_AFM = 21 # J/m^3
        AFM.param.K2_AFM = 0
    AFM.param.K3_AFM = 0
    AFM.param.cHa = 1
    # Different anisotropies
    if ani == 'OOP':
        AFM.param.ea1 = (0,0,1) # Set it z-direction
        AFM.setangle(0,90) # Just move m to z-direction, not necessary to wait every time
    elif ani == 'IP':
        AFM.param.ea1 = (1,0,0)
    else:
        print('Choose anisotropy direction')
        exit()

    folder_name = type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    Ms = 2.1e3

    output_file = path + type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/sinc_' + 'dir' + dir + '_axis' + axis + 'groundstate' + '_damping' + str(damping) + '_T' + str(T) +  '_dispersion.txt'

    ns.reset()

    time = 0.0
    ns.cuda(1)

    ns.dp_newfile(output_file)

    AFM.pbc('x', 5)

    equiltime = 0
    total_time = 2 *t*1e-12
    H0 = 0
    He = 500e3

    ns.setstage('Hequation')
    ns.editstagevalue(0, 'H0', 'He *sinc(kc*(x-Lx/2))*sinc(kc*(y-Ly/2))*sinc(2*PI*fc*(t-t0)),He *sinc(kc*(x-Lx/2))*sinc(kc*(y-Ly/2))*sinc(2*PI*fc*(t-t0))')

    N = 600
    L = meshdims[0]
    kc = 2*np.pi*N/(2*L)
    fc = 5e12
    time_step = 0.1e-12

    ns.equationconstants('H0', H0)
    ns.equationconstants('He', He)
    ns.equationconstants('kc', kc)
    ns.equationconstants('fc', fc)
    ns.equationconstants('t0', t)

    while time < total_time:
        
        ns.editstagestop(0, 'time', time + time_step)
        ns.Run()

        # if axis == 'x':
        ns.dp_getexactprofile((np.array([cellsize/2, 25, meshdims[2]-cellsize])*1e-9), (np.array([meshdims[0] - cellsize/2, 25, meshdims[2]])*1e-9), cellsize*1e-9, 0)
        # elif axis == 'z':
        #     ns.dp_getexactprofile((np.array([meshdims[0]/2-cellsize/2, meshdims[1]/2-cellsize/2, meshdims[2]-cellsize/2])*1e-9), (np.array([meshdims[0]/2-cellsize/2, meshdims[1]/2-cellsize/2, 40+cellsize/2])*1e-9), cellsize*1e-9, 0)
        ns.dp_div(int_dir, Ms)
        ns.dp_saveappendasrow(output_file, int_dir)
        time += time_step

    plotting.plot_magnon_dispersion_with_zoom(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, dir, axis, False, True)