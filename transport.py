import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

from NetSocks import NSClient   # type: ignore
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import plotting

path = 'C:/Users/mathimyh/master/master-boris/'

# Initializes an AFM mesh with hematite parameters. Returns the mesh
def Init_AFM(meshdims, cellsize, damping, MEC, ani, T):

    ns = NSClient(); ns.configure(True, False)
    ns.cuda(1)
    ns.selectcudadevice([0,1])
    ns.reset()
    
    modules = ['exchange', 'aniuni']
    if MEC:
        modules.append('melastic')

    # Set up the antiferromagnet
    AFM = ns.AntiFerromagnet(np.array(meshdims)*1e-9, [cellsize*1e-9])
    AFM.modules(modules)
    temp = str(T) + 'K'
    ns.temperature(temp)
    ns.setode('sLLG', 'RK4') # Stochastic LLG for temperature effects
    ns.setdt(1e-15) # 1fs is good time step

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

    # Add increased damping at edges along x-axis to prevent reflections
    AFM.param.damping_AFM.setparamvar('abl_tanh', [300/meshdims[0], 300/meshdims[0], 0, 0, 0, 0, 1, 10, 300]) 
    
    # Add the magnetoelastic parameters if necessary
    if MEC:
        AFM.surfacefix('-z')
        AFM.seteldt(1e-15)
        AFM.mcellsize([cellsize*1e-9]) 
        AFM.param.cC = (36e10, 17e10, 8.86e10) # N/m^2       A. Yu. Lebedev et al (1989)
        AFM.param.density = 5250 #kg/m^3       found this on google
        AFM.param.MEc = (-3.44e6, 7.5e6) #J/m^3  (Original B2 = 7.5e6)   G. Wedler et al (1999) 
        AFM.param.mdamping = 1e15 # Should probably be lower than this 

    # Relax for 1 ps to get some fluctuations
    ns.Relax(['time', 1e-12])

    # Return the mesh, ready for simulations
    return AFM

# Sets up a simulation with a virtual current
def virtual_current(meshdims, cellsize, damping, MEC, ani, T, type):

    ns = NSClient(); ns.configure(True, False)
    
    # Retrieve the desired mesh
    if type == 'AFM':
        M = Init_AFM(meshdims, cellsize, damping, MEC, ani, T)
    elif type == 'FM':
        M = 0

    # Add the transport modules
    modules = ['exchange', 'aniuni', 'SOTfield', 'transport']
    if MEC:
        modules.append('melastic')
    M.modules(modules)
    ns.clearelectrodes()
    ns.reset()

    # Set spesific params for torque
    M.param.SHA = 1
    M.param.flST = 1

    # If OOP ani, we need to change the spin current direction
    if ani == 'OOP':
        M.param.STp = (1,0,0) # x-dir spin current and y-dir electrode gives z-dir torque

    # Current along y-direction
    ns.addelectrode(np.array([(meshdims[0]/2 - 100), 0, (meshdims[2]-cellsize), (meshdims[0]/2 + 100), 0, meshdims[2]])* 1e-9)
    ns.addelectrode(np.array([(meshdims[0]/2 - 100), meshdims[1], (meshdims[2]-cellsize), (meshdims[0]/2 + 100), meshdims[1], meshdims[2]]) * 1e-9)
    ns.designateground('1')

    # # This is current along x-direction
    # # elif ani == 'OOP':
    # ns.addelectrode(np.array([0, 0, 0, 0, meshdims[1], meshdims[2]]) * 1e-9)
    # ns.addelectrode(np.array([meshdims[0], 0, 0, meshdims[0], meshdims[1], meshdims[2]]) * 1e-9)
    # ns.designateground('1')

    # # Electrode along z-direction
    # ns.addelectrode(np.array([meshdims[0]/2 - 100, 0, 0, meshdims[0]/2 + 100, meshdims[1], 0]) * 1e-9)
    # ns.addelectrode(np.array([meshdims[0]/2 - 100, 0, meshdims[2], meshdims[0]/2 + 100, meshdims[1], meshdims[2]]) * 1e-9)
    # ns.designateground('1')
    
    # else:
    #     print('Which anisotropy?')
    #     exit(1)
    
    # Add step function so that torque only acts on region in the injector
    width = 40
    func = '(step(x-' + str(meshdims[0]/2 - width/2) + 'e-9)-step(x-' + str(meshdims[0]/2 + width/2) + 'e-9)) * (step(z-' + str(meshdims[2]-cellsize) + 'e-9)-step(z-' + str(meshdims[2]) + 'e-9))'
    M.param.SHA.setparamvar('equation', func)
    M.param.flST.setparamvar('equation',func)

    # # Maybe try periodic boundary conditions for the large one, instead of damping equation?
    # ns.pbc('base', 'x')

    return M

# Runs a simulation from ground state for a given time
# Saves the simulation after
# Can plot magnetization at x_vals to find plateau 
def save_steadystate(meshdims, cellsize, t, V, damping, MEC, ani, T, type, x_vals=False):
    
    # Folder system
    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    modules_folder += '/'

    ns = NSClient(); ns.configure(True, False)
    ns.reset()

    M = virtual_current(meshdims, cellsize, damping, MEC, ani, T, type)
    ns.iterupdate(200)

    savename = path + type + '/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K_steady_state.bsm'
    folder_name2 = type + '/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name2):
            os.makedirs(folder_name2)

    # Run the simulation while also saving <mxdmdt> at x_vals every 5ps
    if x_vals != False:
        
        x_vals_string = 'nm_'.join(str(x_val) for x_val in x_vals)

        filename = 'C:/Users/mathimyh/master/master-boris/' + type + '/' + modules_folder + ani + '/cache/plateau/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/plateau_V'  + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K.txt'
        
        folder_name = type + '/' + modules_folder + ani + '/cache/plateau/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        data = ['time']

        for x_val in x_vals:
            data.append(['<mxdmdt>', M, np.array([x_val, 0, meshdims[2], x_val + 1, meshdims[1], meshdims[2]]) * 1e-9])
            
        ns.setsavedata(filename, *data)

        ns.V([0.001*V, 'time', t*1e-12, 'time', 5e-12])
        ns.savesim(savename)
        plotting.plot_plateau(meshdims, V, damping, x_vals, MEC, ani, T, type)

    # Just run the simulation
    else:
        ns.V([0.001*V, 'time', t*1e-12])
        ns.savesim(savename)

# Loads a simulation in steady state, runs the simulation and saves time and <mxdmdt> along the x-axis
def time_avg_SA(meshdims, cellsize, t, V, damping, MEC, ani, T, type, x_start, x_stop):

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder = '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


    # sim_name = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/sims/' + mec_folder + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_steady_state.bsm'
    sim_name = path + type + '/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K_steady_state.bsm'

    ns = NSClient(); ns.configure(True, False)
    ns.reset()
    
    # Loading the sim. All the parameters and parameters variation is still there so don't need to add back
    ns.loadsim(sim_name)
    ns.reset()

    ns.setdata('time')
    for i in range(int((x_stop - x_start)/cellsize)):
        temp = np.array([x_start + (1*i*cellsize), 0, meshdims[2], x_start + (1 + i)*cellsize, meshdims[1], meshdims[2]]) * 1e-9 # Only measure at the top
        ns.adddata('<mxdmdt>', type, temp)

    savename = path + type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'

    ns.savedatafile(savename)

    # ns.cuda(1)

    # Voltage stage
    ns.V([0.001*V, 'time', t*1e-12, 'time', t*1e-12 / 200])

    # ns.Run()

    plotting.plot_tAvg_SA(meshdims, cellsize, t, V, damping, MEC, ani, T, type, x_start, x_stop)
 
# Save 2D magnetization
def time_avg_SA_2D(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani):
    savedata = data[1:-1]
    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/cache/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


    sim_name = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/sims/' + mec_folder + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_steady_state.bsm'
    
    ns = NSClient(); ns.configure(True, False)
    ns.reset()
    
    ns.loadsim(sim_name)
    ns.reset()

    # Voltage stage
    ns.setstage('V')

    ns.editstagevalue('0', str(0.001*V))
    
    ns.editstagestop(0, 'time', t * 1e-12)

    ns.editdatasave(0, 'time', t * 1e-12 /200)

    ns.setdata('time')


    # To not have an enormous amount of data, x-direction will only sample every 10th cellsize. 
    # z-direction will sample every nm.
    # At least for now, I will see if the resolution is fine or not. 
    for j in range(meshdims[2]):
        for i in range(int((x_stop - x_start)/cellsize*0.1)): 
            temp = np.array([x_start + i*cellsize*10, 0, j, x_start + (i+1)*cellsize*10, meshdims[1], j]) * 1e-9 # Average over y direction
            ns.adddata(data, "base", temp)


    savename = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + savedata  + '.txt'

    ns.savedatafile(savename)

    ns.cuda(1)
    # ns.selectcudadevice([0,1])

    ns.Run()

    plotting.plot_tAvg_SA_2D(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani)

def time_avg_SA_2D_y(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani):
    savedata = data[1:-1]
    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/cache/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


    sim_name = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/sims/' + mec_folder + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_steady_state.bsm'
    
    ns = NSClient(); ns.configure(True, False)
    ns.reset()
    
    ns.loadsim(sim_name)
    ns.reset()

    # Voltage stage
    ns.setstage('V')

    ns.editstagevalue('0', str(0.001*V))
    
    ns.editstagestop(0, 'time', t * 1e-12)

    ns.editdatasave(0, 'time', t * 1e-12 /200)

    ns.setdata('time')


    # To not have an enormous amount of data, x-direction will only sample every 10th cellsize. 
    # z-direction will sample every nm.
    # At least for now, I will see if the resolution is fine or not. 
    for j in range(meshdims[1]):
        for i in range(int((x_stop - x_start)/cellsize*0.1)): 
            temp = np.array([x_start + i*cellsize*10, j, meshdims[2], x_start + (i+1)*cellsize*10, j, meshdims[2]]) * 1e-9 # Average over y direction
            ns.adddata(data, "base", temp)


    savename = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/ydir_2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + savedata  + '.txt'

    ns.savedatafile(savename)

    ns.cuda(1)
    # ns.selectcudadevice([0,1])

    ns.Run()

    plotting.plot_tAvg_SA_2D_y(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani)

def time_avg_SA_z(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani):
    savedata = data[1:-1]
    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/cache/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


    sim_name = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/sims/' + mec_folder + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_steady_state.bsm'
    
    ns = NSClient(); ns.configure(True, False)
    ns.reset()
    
    ns.loadsim(sim_name)
    ns.reset()

    # Voltage stage
    ns.setstage('V')

    ns.editstagevalue('0', str(0.001*V))
    
    ns.editstagestop(0, 'time', t * 1e-12)

    ns.editdatasave(0, 'time', t * 1e-12 /200)

    ns.setdata('time')
    
    for p in range(meshdims[2]):
        temp = np.array([meshdims[0]/2 - cellsize, 0, meshdims[2]-p, meshdims[0]/2 + cellsize, meshdims[1], meshdims[2]-p]) * 1e-9
        ns.adddata(data, "base", temp)


    savename = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + savedata  + '_zdir.txt'

    ns.savedatafile(savename)

    ns.cuda(1)
    # ns.selectcudadevice([0,1])

    ns.Run()

    plotting.plot_tAvg_SA_z(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani)

# Get a profile of the magnetization
def profile_from_sim(t, V, damping, sim_name, x_start, x_stop):

    ns = NSClient(); ns.configure(True)
    ns.reset()
    
    ns.loadsim(sim_name)
    ns.reset()

    # Voltage stage
    ns.setstage('V')

    ns.editstagevalue('0', str(0.001*V))
    
    ns.editstagestop(0, 'time', t * 1e-12)

    ns.setdata("commbuf")
    ns.adddata("time")

    start = str(x_start) + 'e-9, 10e-9, 0'
    end = str(x_stop) + 'e-9, 10e-9, 0'

    savedt = 1e-12

    for i in range(0, 6):
        ns.editdatasave(i, "time", savedt)

    ns.dp_getexactprofile(start = start, end = end, step = '4e-9', dp_index = '0', bufferCommand = True)
    ns.dp_save("C:/Users/mathimyh/Documents/Boris data/Simulations/boris_fordypningsoppgave/cache/profile_test.txt", dp_indexes = 1, bufferCommand = True)

    ns.cuda(1)

    ns.Run()