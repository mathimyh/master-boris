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
    
    try:
        ns.cuda(1)
    except:
        print('No cuda!')

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
    Saves the system in steady state as a .bsm file, so that we dont need to do it more 
    than once. 
    
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

def current_density(ns, currentDensity):

    '''
    
    Find the current density (at the position of the injector) for a system
    during spin pumping. 
    Current density is proportional to applied voltage so not very useful
    
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

def time_avg_SA(ns, timeAvgSA, sim_num=False, show_plot=False):

    '''

    Loads a system in steady state and saves <mxdmdt> over time along the x-axis

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

    try:
        ns.cuda(1)
    except:
        print('No cuda!')

    # Voltage stage
    ns.V([0.001*timeAvgSA.V, 'time', timeAvgSA.t*1e-12, 'time', timeAvgSA.t*1e-12 / 200])

    # Plot the result
    plotting.plot_tAvg_SA(timeAvgSA, sim_num,show_plot)