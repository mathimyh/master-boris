import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import plotting
import params

def Init_AFM(ns, params):
    
    '''

    Initializes an AFM mesh with the given parameters and modules. 
    Standard is to start with magnetization along easy axis and relax for 10ps to get fluctuations.
    
    Returns the mesh object
    
    '''

    # Different modules depending on params
    modules = ['exchange']
    if params.hard_axis:
        modules.append('anitens')
    else:
        modules.append('aniuni')
    if params.MEC:
        modules.append('melastic')
    if params.Hfield:
        modules.append('Zeeman')

    # Set up the antiferromagnet
    AFM = ns.AntiFerromagnet(np.array(params.meshdims)*1e-9, [params.cellsize*1e-9])
    AFM.modules(modules)
    temp = str(params.T) + 'K'
    ns.temperature(temp)
    ns.setode('sLLG', 'RK4') # Stochastic LLG for temperature effects
    ns.setdt(1e-15)

    # Set material parameters    
    AFM.param.grel_AFM = 1
    AFM.param.damping_AFM =  params.damping
    AFM.param.Ms_AFM = 2.1e3
    AFM.param.Nxy = 0
    AFM.param.A_AFM = 76e-15 # J/m
    AFM.param.Ah = -460e3 # J/m^3
    AFM.param.Anh = 0.0
    AFM.param.J1 = 0
    AFM.param.J2 = 0
    if params.hard_axis:
        AFM.param.K1_AFM = 21 # J/m^3
        ns.setktens('-0.001x2', '1z2')
    else:
        AFM.param.K1_AFM = 21 # J/m^3
        AFM.param.K2_AFM = 0
    AFM.param.K3_AFM = 0
    AFM.param.cHa = 1
    # Different anisotropies
    if params.ani == 'OOP':
        AFM.param.ea1 = (0,0,1) # Set it z-direction
        AFM.setangle(0,90) # Just move m to z-direction, not necessary to wait every time
    elif params.ani == 'IP':
        AFM.param.ea1 = (1,0,0)
        # AFM.setangle(90,0)
    else:
        print('Choose anisotropy direction')
        exit()

    # Add increased damping at edges along x-axis to prevent reflections. For smaller meshes this needs to be shortened
    damping_x = 300
    if params.meshdims[0] == 1000:
        damping_x = 50
    AFM.param.damping_AFM.setparamvar('abl_tanh', [damping_x/params.meshdims[0], damping_x/params.meshdims[0], 0, 0, 0, 0, 1, 10, damping_x]) 
    
    # Add the magnetoelastic parameters if necessary
    if params.MEC:
        AFM.surfacefix('-z')
        AFM.seteldt(1e-15)
        AFM.mcellsize([params.cellsize*1e-9]) 
        AFM.param.cC = (36e10, 17e10, 8.86e10) # N/m^2       A. Yu. Lebedev et al (1989)
        AFM.param.density = 5250 #kg/m^3       found this on google
        AFM.param.MEc = (-3.44e6, 7.5e6) #J/m^3  (Original B2 = 7.5e6)   G. Wedler et al (1999) 
        AFM.param.mdamping = 1e15 # Should probably be lower than this 

    # ns.random()

    # Critical field, found from calculations
    AFM.setfield(params.Hfield, 90, 90)
    
    ns.cuda(1)

    # Relax to reach plateau of fluctuations
    relax_time = 10e-12
    if params.hard_axis and params.Hfield > 0: # Needs longer time with these conditions
        relax_time = 50e-12
    ns.Relax(['time', relax_time])
    

    # Return the mesh, ready for simulations
    return AFM

def Init_FM(ns, params):

    '''
    
    Initializes a FM mesh with the given parameters and modules. 
    Standard is to start with magnetization along easy axis and relax for 1000ps to get fluctuations.
    
    Returns the mesh object and saves the simulation.
    
    '''

    ns.iterupdate(200)
    
    modules = ['exchange', 'aniuni']
    if params.MEC:
        modules.append('melastic')

    # Set up the antiferromagnet
    FM = ns.Ferromagnet(np.array(params.meshdims)*1e-9, [params.cellsize*1e-9])
    FM.modules(modules)
    temp = str(params.T) + 'K'
    ns.temperature(temp)
    ns.setode('sLLG', 'RK4') # Stochastic LLG for temperature effects
    ns.setdt(1e-15) 

    # Set parameters. Should possibly just save these in the database really    
    FM.param.grel = 1
    FM.param.damping =  params.damping
    FM.param.Ms = 2.1e3
    FM.param.Nxy = 0
    FM.param.A = 76e-15 # J/m
    FM.param.J1 = 0
    FM.param.J2 = 0
    FM.param.K1 = 2100 # J/m^3
    FM.param.K2 = 0
    FM.param.K3 = 0
    FM.param.cHa = 1
    # Different anisotropies
    if params.ani == 'OOP':
        FM.param.ea1 = (0,0,1) # Set it z-direction
        FM.setangle(0,90) # Just move m to z-direction, not necessary to wait every time
    elif params.ani == 'IP':
        FM.param.ea1 = (1,0,0)
    else:
        print('Choose anisotropy direction')
        exit()

    # Add increased damping at edges along x-axis to prevent reflections
    FM.param.damping.setparamvar('abl_tanh', [300/params.meshdims[0], 300/params.meshdims[0], 0, 0, 0, 0, 1, 10, 300]) 
    
    # Add the magnetoelastic parameters if necessary
    if params.MEC:
        FM.surfacefix('-z')
        ns.seteldt(1e-15)
        FM.mcellsize([params.cellsize*1e-9]) 
        FM.param.cC = (36e10, 17e10, 8.86e10) # N/m^2       A. Yu. Lebedev et al (1989)
        FM.param.density = 5250 #kg/m^3       found this on google
        FM.param.MEc = (-3.44e6, 7.5e6) #J/m^3   G. Wedler et al (1999) 
        FM.param.mdamping = 1e15 # Should probably be lower than this 

    # FM.pbc(['x', 10])

    ns.random()

    # Relax for 1000 ps to thermally excite magnons
    ns.Relax(['time', 1000e-12])

    # Return the mesh, ready for simulations
    # sim_name = path + 'FM/' + modules_folder + params.ani + '/sims/' +  str(params.meshdims[0]) + 'x' + str(params.meshdims[1]) + 'x' + str(params.meshdims[2]) + '/ground_state.bsm'
    # ns.savesim(sim_name)
    return FM

def virtual_current(ns, params):
    
    '''
    
    Retrieves a ground state mesh from the Init function and adds all necessary modules and parameters
    for a virtual current simulation. 

    Returns the mesh.
    
    '''


    # Retrieve the desired mesh
    if params.type == 'AFM':
        M = Init_AFM(ns, params)
    elif params.type == 'FM':
        M = Init_FM(params)
    else:
        print('Choose type!')
        exit()

    # Add the transport modules
    modules = ['exchange', 'aniuni', 'SOTfield', 'transport']
    if params.MEC:
        modules.append('melastic')
    M.modules(modules)
    ns.clearelectrodes()
    ns.reset()

    # Set spesific params for torque
    M.param.SHA = 1
    if params.type == 'FM':
        M.param.flST = 0
    else:
        M.param.flST = 1

    # If OOP ani, we need to change the spin current direction
    if params.ani == 'OOP':
        M.param.STp = (1,0,0) # x-dir spin current and y-dir electrode gives z-dir torque

    meshdims = params.meshdims
    cellsize = params.cellsize

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

def save_steadystate(ns, steadystate, sim_num=False):

    '''
    
    Start with a system in its ground state and find steady state for transport. 
    Can save <mxdmdt> as a function of time to make a plateau plot.
    Saves the system in steady state.
    
    '''


    M = virtual_current(ns, steadystate)
    # ns.iterupdate(200)

    if not sim_num:
        savename = steadystate.simname()
    else:
       savename = steadystate.simname()[:-4]
       savename += '_sim' + str(sim_num) + '.bsm' 
    params.make_folder(savename)

    # Run the simulation while also saving <mxdmdt> at x_vals every 5ps
    if steadystate.x_vals != False:

        filename = steadystate.cachename()
        params.make_folder(filename)
        
        data = ['time']

        for x_val in steadystate.x_vals:
            data.append(['<mxdmdt>', M, np.array([x_val, 0, steadystate.meshdims[2], x_val + 1, steadystate.meshdims[1], steadystate.meshdims[2]]) * 1e-9])
            
        ns.setsavedata(filename, *data)

        ns.V([0.001*steadystate.V, 'time', steadystate.t*1e-12, 'time', 5e-12])
        ns.savesim(savename)
        plotting.plot_plateau(steadystate)

    # Just run the simulation
    else:
        ns.V([0.001*steadystate.V, 'time', steadystate.t*1e-12])
        ns.savesim(savename)

def save_steadystate2(ns, steadystate):
    
    # Folder system
    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    if hard_axis:
        modules_folder += '+hard_axis'
    modules_folder += '/'

    sim_name = path + type + '/' + modules_folder + ani + '/sims/' +  str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/ground_state.bsm'
    ns.loadsim(sim_name)    

    ns.addmodule("FM", "SOTfield")
    ns.addmodule("FM", "transport")
    ns.setparam("FM", "SHA", '1')
    ns.setparam("FM", "flST", '0')
    ns.clearelectrodes()
    ns.reset()

    # # If OOP ani, we need to change the spin current direction
    # if ani == 'OOP':
    #     M.param.STp = (1,0,0) # x-dir spin current and y-dir electrode gives z-dir torque

    # Current along y-direction
    ns.addelectrode(np.array([(meshdims[0]/2 - 100), 0, (meshdims[2]-cellsize), (meshdims[0]/2 + 100), 0, meshdims[2]])* 1e-9)
    ns.addelectrode(np.array([(meshdims[0]/2 - 100), meshdims[1], (meshdims[2]-cellsize), (meshdims[0]/2 + 100), meshdims[1], meshdims[2]]) * 1e-9)
    ns.designateground('1')

    ns.iterupdate(200)

    # Add step function so that torque only acts on region in the injector
    width = 20
    func = '(step(x-' + str(meshdims[0]/2 - width/2) + 'e-9)-step(x-' + str(meshdims[0]/2 + width/2) + 'e-9)) * (step(z-' + str(meshdims[2]-cellsize) + 'e-9)-step(z-' + str(meshdims[2]) + 'e-9))'
    ns.setparamvar('SHA','equation', func)
    ns.setparamvar('flST','equation',func)

    savename = path + type + '/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K_steady_state.bsm'
    folder_name2 = type + '/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name2):
            os.makedirs(folder_name2)

    # Run the simulation while also saving <mxdmdt> at x_vals every 5ps
    if x_vals != False:

        filename = 'C:/Users/mathimyh/master/master-boris/' + type + '/' + modules_folder + ani + '/cache/plateau/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/plateau_V'  + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K.txt'
        
        folder_name = type + '/' + modules_folder + ani + '/cache/plateau/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        ns.setdata('time')

        for x_val in x_vals:
            temp =  np.array([x_val, 0, meshdims[2], x_val + 1, meshdims[1], meshdims[2]]) * 1e-9
            ns.adddata('<mxdmdt>', 'FM', temp)

        ns.savedatafile(filename)
        ns.V([0.001*V, 'time', t*1e-12, 'time', 5e-12])
        ns.savesim(savename)
        plotting.plot_plateau(meshdims, V, damping, x_vals, MEC, ani, T, type, hard_axis)

    # Just run the simulation
    else:
        ns.V([0.001*V, 'time', t*1e-12])
        ns.savesim(savename)

# Run transport simulation and save current density over time (all 3 components)
def current_density(ns, currentDensity):

    '''
    
    Find the current density (at the position of the injector) for a system. 
    
    
    '''

    M = virtual_current(ns, currentDensity)


    output_file = currentDensity.cachename()
    print(output_file)
    params.make_folder(output_file)

    # See if there is difference in spin current within the injector. Shouldnt be
    one = np.array([currentDensity.meshdims[0]/2 - 20, 0, currentDensity.meshdims[2], currentDensity.meshdims[0]-20, currentDensity.meshdims[1], currentDensity.meshdims[2]])*1e-9
    two = np.array([currentDensity.meshdims[0]/2 - 20, 0, currentDensity.meshdims[2], currentDensity.meshdims[0]-20, currentDensity.meshdims[1], currentDensity.meshdims[2]])*1e-9
    three = np.array([currentDensity.meshdims[0]/2 - 20, 0, currentDensity.meshdims[2], currentDensity.meshdims[0]-20, currentDensity.meshdims[1], currentDensity.meshdims[2]])*1e-9
        
    ns.setsavedata(output_file, ['<Jc>', M, one], ['<Jc>', M, two], ['<Jc>', M, three])

    ns.V([0.001*currentDensity.V, 'time', currentDensity.t*1e-12, 'time', 1e-12])

 
def time_avg_SA(ns, timeAvgSA,sim_num=False):

    '''

    Loads a system in steady state and saves <mxdmdt> over time along the x-axis
    given by x_start and x_stop

    '''


    if not sim_num:
        sim_name = timeAvgSA.simname()
    else:
       sim_name = timeAvgSA.simname()[:-4]
       sim_name += '_sim' + str(sim_num) + '.bsm' 
    
    # Loading the sim. All the parameters and parameters variation is still there so don't need to add back
    ns.loadsim(sim_name)
    ns.reset()

    ns.setdata('time')
    
    if timeAvgSA.direction == 'x':
        for i in range(int((timeAvgSA.x_stop - timeAvgSA.x_start)/timeAvgSA.cellsize)):
            # Only measure at the top
            temp = np.array([timeAvgSA.x_start + (i*timeAvgSA.cellsize), 0, timeAvgSA.meshdims[2], 
                            timeAvgSA.x_start + (1 + i)*timeAvgSA.cellsize, timeAvgSA.meshdims[1], 
                            timeAvgSA.meshdims[2]]) * 1e-9 
            ns.adddata('<mxdmdt>', timeAvgSA.type, temp)
    
    elif timeAvgSA.direction == 'z':
        for i in range(int(timeAvgSA.meshdims[2]/timeAvgSA.cellsize)):
            # Measure through middle of sample
            temp = np.array([timeAvgSA.meshdims[0]/2, 0, timeAvgSA.cellsize*i, 
                             timeAvgSA.meshdims[0]/2, timeAvgSA.meshdims[1], 
                             timeAvgSA.cellsize*(i+1)]) * 1e-9 
            ns.adddata('<mxdmdt>', timeAvgSA.type, temp)
            
    elif timeAvgSA.direction == 'y':
        # Only measure at the top
        for i in range(int(timeAvgSA.meshdims[1]/timeAvgSA.cellsize)):
            temp = np.array([timeAvgSA.meshdims[0]/2, timeAvgSA.cellsize*i, timeAvgSA.meshdims[2],
                             timeAvgSA.meshdims[0]/2], timeAvgSA.cellsize*(i+1), 
                            timeAvgSA.meshdims[2])
            ns.adddata('<mxdmdt>', timeAvgSA.type, temp)
            
    if not sim_num:
        savename = timeAvgSA.cachename()
    else:
        savename = timeAvgSA.cachename()[:-4]
        savename += '_sim' + str(sim_num) + '.txt'
    params.make_folder(savename)
    ns.savedatafile(savename)

    ns.cuda(1)

    # Voltage stage
    ns.V([0.001*timeAvgSA.V, 'time', timeAvgSA.t*1e-12, 'time', timeAvgSA.t*1e-12 / 200])

    plotting.plot_tAvg_SA(timeAvgSA, sim_num)
 
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

def time_avg_SA_underneath(ns, meshdims, cellsize, t, V, damping, MEC, ani, T, type, x_start, x_stop):

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder = '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/underneath/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    # sim_name = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/sims/' + mec_folder + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_steady_state.bsm'
    sim_name = path + type + '/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K_steady_state.bsm'

    # Loading the sim. All the parameters and parameters variation is still there so don't need to add back
    ns.loadsim(sim_name)
    ns.reset()

    ns.setdata('time')
    ns.adddata('<mxdmdt>', type, np.array([x_start, 0, meshdims[2], x_start + 1, meshdims[1], meshdims[2]]) * 1e-9)
    for i in range(int((x_stop - x_start)/cellsize)):
        temp = np.array([x_start + (1*i*cellsize), 0, 0, x_start + (1 + i)*cellsize, meshdims[1], cellsize]) * 1e-9 # Only measure at the bottom
        ns.adddata('<mxdmdt>', type, temp)

    savename = path + type + '/' + modules_folder + ani + '/underneath/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'

    ns.savedatafile(savename)

    # ns.cuda(1)

    # Voltage stage
    ns.V([0.001*V, 'time', t*1e-12, 'time', t*1e-12 / 200])

    plotting.plot_tAvg_SA_underneath(meshdims, cellsize, t, V, damping, MEC, ani, T, type, x_start,x_stop)

def complete_simulation_AFM(ns, meshdims, cellsize, t_steady, t_saving, V, damping, MEC, ani, T, hard_axis):
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
    if hard_axis:
        AFM.param.K1_AFM = -21e-3 # J/m^3
        AFM.param.K2_AFM = 21
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

    # Add increased damping at edges along x-axis to prevent reflections. For short meshes this needs to be shortened
    damping_x = 300
    if meshdims[0] == 1000:
        damping_x = 50
    AFM.param.damping_AFM.setparamvar('abl_tanh', [damping_x/meshdims[0], damping_x/meshdims[0], 0, 0, 0, 0, 1, 10, damping_x]) 
    
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
    ns.cuda(1)
    ns.Relax(['time', 1e-12])

    # Add the transport modules
    modules = ['exchange', 'aniuni', 'SOTfield', 'transport']
    if MEC:
        modules.append('melastic')
    AFM.modules(modules)
    ns.clearelectrodes()
    ns.reset()

    # Set spesific params for torque
    AFM.param.SHA = 1
    AFM.param.flST = 1

    # If OOP ani, we need to change the spin current direction
    if ani == 'OOP':
        AFM.param.STp = (1,0,0) # x-dir spin current and y-dir electrode gives z-dir torque

    # Current along y-direction
    ns.addelectrode(np.array([(meshdims[0]/2 - 100), 0, (meshdims[2]-cellsize), (meshdims[0]/2 + 100), 0, meshdims[2]])* 1e-9)
    ns.addelectrode(np.array([(meshdims[0]/2 - 100), meshdims[1], (meshdims[2]-cellsize), (meshdims[0]/2 + 100), meshdims[1], meshdims[2]]) * 1e-9)
    ns.designateground('1')

    # Add step function so that torque only acts on region in the injector
    width = 40
    func = '(step(x-' + str(meshdims[0]/2 - width/2) + 'e-9)-step(x-' + str(meshdims[0]/2 + width/2) + 'e-9)) * (step(z-' + str(meshdims[2]-cellsize) + 'e-9)-step(z-' + str(meshdims[2]) + 'e-9))'
    AFM.param.SHA.setparamvar('equation', func)
    AFM.param.flST.setparamvar('equation',func)

    # Folder system
    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    if hard_axis:
        modules_folder += '+hard_axis'
    modules_folder += '/'

    ns.reset()
    ns.iterupdate(200)

    savename = path + 'AFM/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K_steady_state.bsm'
    folder_name2 = 'AFM/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name2):
            os.makedirs(folder_name2)

    ns.V([0.001*V, 'time', t_steady*1e-12])
    ns.savesim(savename)

    folder_name = 'AFM/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    ns.reset()

    x_start = int(meshdims[0]/2 + 20)
    x_stop = meshdims[0]

    ns.setdata('time')
    for i in range(int((x_stop - x_start)/cellsize)):
        temp = np.array([x_start + (1*i*cellsize), 0, meshdims[2], x_start + (1 + i)*cellsize, meshdims[1], meshdims[2]]) * 1e-9 # Only measure at the top
        ns.adddata('<mxdmdt>', 'AFM', temp)

    savename = path + 'AFM/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'

    ns.savedatafile(savename)

    # Voltage stage
    ns.V([0.001*V, 'time', t_saving*1e-12, 'time', t_saving*1e-12 / 200])

    plotting.plot_tAvg_SA(meshdims, cellsize, t_saving, V, damping, MEC, ani, T, 'AFM', hard_axis, x_start, x_stop)

def complete_simulation_AFM_z(ns, meshdims, cellsize, t_steady, t_saving, V, damping, MEC, ani, T, hard_axis):
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
    if hard_axis:
        AFM.param.K1_AFM = -21e-3 # J/m^3
        AFM.param.K2_AFM = 21
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

    # Add increased damping at edges along x-axis to prevent reflections. For short meshes this needs to be shortened
    damping_x = 300
    if meshdims[0] == 1000:
        damping_x = 50
    AFM.param.damping_AFM.setparamvar('abl_tanh', [damping_x/meshdims[0], damping_x/meshdims[0], 0, 0, 0, 0, 1, 10, damping_x]) 
    
    # Add the magnetoelastic parameters if necessary
    if MEC:
        AFM.surfacefix('-z')
        AFM.seteldt(1e-15)
        AFM.mcellsize([cellsize*1e-9]) 
        AFM.param.cC = (36e10, 17e10, 8.86e10) # N/m^2       A. Yu. Lebedev et al (1989)
        AFM.param.density = 5250 #kg/m^3       found this on google
        AFM.param.MEc = (-3.44e6, 7.5e6) #J/m^3  (Original B2 = 7.5e6)   G. Wedler et al (1999) 
        AFM.param.mdamping = 1e15 # Should probably be lower than this 

    # Relax for 10 ps to get some fluctuations
    ns.cuda(1)
    ns.Relax(['time', 10e-12])

    # Add the transport modules
    modules = ['exchange', 'aniuni', 'SOTfield', 'transport']
    if MEC:
        modules.append('melastic')
    AFM.modules(modules)
    ns.clearelectrodes()
    ns.reset()

    # Set spesific params for torque
    AFM.param.SHA = 1
    AFM.param.flST = 1

    # If OOP ani, we need to change the spin current direction
    if ani == 'OOP':
        AFM.param.STp = (1,0,0) # x-dir spin current and y-dir electrode gives z-dir torque

    # Current along y-direction
    ns.addelectrode(np.array([(meshdims[0]/2 - 100), 0, (meshdims[2]-cellsize), (meshdims[0]/2 + 100), 0, meshdims[2]])* 1e-9)
    ns.addelectrode(np.array([(meshdims[0]/2 - 100), meshdims[1], (meshdims[2]-cellsize), (meshdims[0]/2 + 100), meshdims[1], meshdims[2]]) * 1e-9)
    ns.designateground('1')

    # Add step function so that torque only acts on region in the injector
    width = 40
    func = '(step(x-' + str(meshdims[0]/2 - width/2) + 'e-9)-step(x-' + str(meshdims[0]/2 + width/2) + 'e-9)) * (step(z-' + str(meshdims[2]-cellsize) + 'e-9)-step(z-' + str(meshdims[2]) + 'e-9))'
    AFM.param.SHA.setparamvar('equation', func)
    AFM.param.flST.setparamvar('equation',func)

    # Folder system
    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    if hard_axis:
        modules_folder += '+hard_axis'
    modules_folder += '/'

    ns.reset()
    ns.iterupdate(200)

    savename = path + 'AFM/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V' + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K_steady_state.bsm'
    folder_name2 = 'AFM/' + modules_folder + ani + '/sims/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name2):
            os.makedirs(folder_name2)

    ns.V([0.001*V, 'time', t_steady*1e-12])
    ns.savesim(savename)

    folder_name = 'AFM/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    ns.reset()

    x_start = int(meshdims[0]/2 + 20)
    x_stop = meshdims[0]

    ns.setdata('time')
    for i in range(int(meshdims[2]/cellsize)):
        temp = np.array([meshdims[0]/2, cellsize*i, 0, meshdims[0]/2, meshdims[1], cellsize*(i+1)]) * 1e-9 # Measure through middle of sample
        ns.adddata('<mxdmdt>', 'AFM', temp)

    savename = path + 'AFM/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/z_dir_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'

    ns.savedatafile(savename)

    # Voltage stage
    ns.V([0.001*V, 'time', t_saving*1e-12, 'time', t_saving*1e-12 / 200])

    plotting.plot_tAvg_SA_z(meshdims, cellsize, t_saving, V, damping, MEC, ani, T, 'AFM')