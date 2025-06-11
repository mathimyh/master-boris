import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import to_rgb, to_hex
import matplotlib.colors as mcolors
from colorspacious import cspace_convert
import colorsys
import seaborn as sns
import numpy as np
import os
import params
import ruptures as rpt
from scipy.optimize import curve_fit
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import copy

plt.rcParams.update({'font.size': 26})


#### MAGNON TRANSPORT PLOTS ###

def plot_plateau(steadystate):

    '''
    
    Plots <mxdmdt> over time for a magnon transport simulation that starts in the ground state. This function
    is used to find the plateau of <mxdmdt>, i.e. when/if steady state is found. The result is several plots
    at different x-positions in the mesh, given by steadystate.x_vals

    '''

    filename = steadystate.cachename()
    f = open(filename, 'r')

    lines = f.readlines()
    skip = False
    if lines[7][0] == 's':
        skip = True
    lines = lines[11:]

    indexer = 0

    fig, ax = plt.subplots(2, 3)

    fig.set_figheight(10)
    fig.set_figwidth(14)

    direction = 1
    if steadystate.ani == 'OOP':
        direction = 3

    for i in range(direction, len(lines[0].split('\t'))-3, 3):
        ts = []
        vals = []

        for line in lines:
            vec = line.split('\t')
            if skip:
                vec = vec[3:]
            ts.append(float(vec[0])*1e12)
            vals.append(float(vec[i]))

        text = str(steadystate.x_vals[indexer]) + 'nm'
        if indexer == 0:
            ax[0][0].plot(ts, vals)
            ax[0][0].title.set_text(text)
        elif indexer == 1:
            ax[0][1].plot(ts, vals)
            ax[0][1].title.set_text(text)
        elif indexer == 2:
            ax[0][2].plot(ts, vals)
            ax[0][2].title.set_text(text)
        elif indexer == 3:
            ax[1][0].plot(ts, vals)
            ax[1][0].title.set_text(text)
        elif indexer == 4:
            ax[1][1].plot(ts, vals)
            ax[1][1].title.set_text(text)
        elif indexer == 5:
            ax[1][2].plot(ts, vals)
            ax[1][2].title.set_text(text)

        indexer += 1


    plotname = steadystate.plotname()
    params.make_folder(plotname)
    fig.suptitle('<mxdmdt> over time')
    fig.tight_layout()
    fig.savefig(plotname, dpi=500)
    plt.show()

def plot_tAvg_SA(timeAvgSA, sim_num=False):

    
    '''
    
    The standard function for plotting <mxdmdt> over distance in magnon transport simulations. 

    '''

    plt.figure(figsize=(10, 7))
    
    color = '#1B9E77'
    if timeAvgSA.hard_axis:
        color = '#D95F02'

    if not sim_num:
        filename = timeAvgSA.cachename()
    else:
        filename = timeAvgSA.cachename()[:-4]
        filename += '_sim' + str(sim_num) + '.txt'

    try:
        f = open(filename, 'r')

    except FileNotFoundError:
        print("No simulations have been done with these params... V = ", timeAvgSA.V)

    else:
        lines = f.readlines()
        lines = lines[10:]

        vals = []


        for line in lines:
            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if timeAvgSA.ani == 'OOP':
                direction = 2
            vec = line.split('\t')
            all_vals = vec[1:]
            temp = []
            while direction < len(all_vals):
                temp.append(float(all_vals[direction]))
                direction += 3
            vals.append(temp)

        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                try:
                    val += float(vals[j][i])
                except IndexError:
                    print('WHAT DA HELL')
                    break
            val /= len(vals)
            ys.append(val)
            
        if timeAvgSA.direction == 'x':
            xs = np.linspace(timeAvgSA.x_start/1000, timeAvgSA.x_stop/1000, len(ys))
            plt.xlabel('x (μm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=3, color=color)
        
        elif timeAvgSA.direction == 'z':
            xs = np.linspace(timeAvgSA.cellsize/2, timeAvgSA.meshdims[2]-timeAvgSA.cellsize/2, len(ys))
            plt.xlabel('z (nm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=2, color=color, marker='s', markersize=10)
            ticks = np.linspace(0, timeAvgSA.meshdims[2], int(timeAvgSA.meshdims[2]/10)+1)
            plt.xticks(ticks)
            plt.ylim(0, 30)
            
        elif timeAvgSA.direction == 'y':
            xs = np.linspace(timeAvgSA.cellsize/2, timeAvgSA.meshdims[1]-timeAvgSA.cellsize/2, len(ys))
            plt.xlabel('y (nm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=2, color=color, marker='s', markersize=10)
            ticks = np.linspace(0, timeAvgSA.meshdims[1], int(timeAvgSA.meshdims[1]/5)+1)
            plt.xticks(ticks)
        
        plt.ylabel(r'$\mu$')
        
            
        plt.tight_layout()
            

        if not sim_num:
            plotname = timeAvgSA.plotname()
        else:
            plotname = timeAvgSA.plotname()[:-4]
            plotname += '_sim' + str(sim_num) + '.png'
        params.make_folder(plotname)
        plt.savefig(plotname, dpi=500)
        
        plt.show()
        
        # plt.close()
 
def plot_tAvg_SA_both_systems(timeAvgSA, sim_num=False):

    
    '''
    
    The standard function for plotting <mxdmdt> over distance in magnon transport simulations. 
    Plots both biaxial and uniaxial systems together in the same plot.
    Possible to pass many timeAvgSA objects together in the same plot

    '''

    plt.figure(figsize=(10, 7))
    
    timeAvgSA.set_hard_axis(0,0)

    if not sim_num:
        filename = timeAvgSA.cachename()
    else:
        filename = timeAvgSA.cachename()[:-4]
        filename += '_sim' + str(sim_num) + '.txt'

    
    color = '#1B9E77'

    try:
        f = open(filename, 'r')

    except FileNotFoundError:
        print("No simulations have been done with these params... V = ", timeAvgSA.V)

    else:
        lines = f.readlines()
        lines = lines[10:]

        vals = []


        for line in lines:
            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if timeAvgSA.ani == 'OOP':
                direction = 2
            vec = line.split('\t')
            all_vals = vec[1:]
            temp = []
            while direction < len(all_vals):
                temp.append(float(all_vals[direction]))
                direction += 3
            vals.append(temp)

        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                try:
                    val += float(vals[j][i])
                except IndexError:
                    print('WHAT DA HELL')
                    break
            val /= len(vals)
            ys.append(val)
            
        if timeAvgSA.direction == 'x':
            xs = np.linspace(timeAvgSA.x_start/1000, timeAvgSA.x_stop/1000, len(ys))
            plt.xlabel('x (μm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=3, color=color)
        
        elif timeAvgSA.direction == 'z':
            xs = np.linspace(timeAvgSA.cellsize/2, timeAvgSA.meshdims[2]-timeAvgSA.cellsize/2, len(ys))
            plt.xlabel('z (nm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=2, color=color, marker='s', markersize=10)
            ticks = np.linspace(0, timeAvgSA.meshdims[2], int(timeAvgSA.meshdims[2]/15)+1)
            # plt.xticks(ticks)
            
        elif timeAvgSA.direction == 'y':
            xs = np.linspace(timeAvgSA.cellsize/2, timeAvgSA.meshdims[1]-timeAvgSA.cellsize/2, len(ys))
            plt.xlabel('y (nm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=2, color=color, marker='s', markersize=10)
            ticks = np.linspace(0, timeAvgSA.meshdims[1], int(timeAvgSA.meshdims[1]/5)+1)
            plt.xticks(ticks)
            
    
    critical_H = 4.711e6
    timeAvgSA.set_hard_axis(1, critical_H)
    
    if not sim_num:
        filename2 = timeAvgSA.cachename()
    else:
        filename2 = timeAvgSA.cachename()[:-4]
        filename2 += '_sim' + str(sim_num) + '.txt'
    
    color = '#D95F02'
    
    try:
        f = open(filename2, 'r')

    except FileNotFoundError:
        print("No simulations have been done with these params... V = ", timeAvgSA.V)

    else:
        lines = f.readlines()
        lines = lines[10:]

        vals = []


        for line in lines:
            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if timeAvgSA.ani == 'OOP':
                direction = 2
            vec = line.split('\t')
            all_vals = vec[1:]
            temp = []
            while direction < len(all_vals):
                temp.append(float(all_vals[direction]))
                direction += 3
            vals.append(temp)

        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                try:
                    val += float(vals[j][i])
                except IndexError:
                    print('WHAT DA HELL')
                    break
            val /= len(vals)
            ys.append(val)
            
        if timeAvgSA.direction == 'x':
            xs = np.linspace(timeAvgSA.x_start/1000, timeAvgSA.x_stop/1000, len(ys))
            plt.xlabel('x (μm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=3, color=color, linestyle=(0, (5, 10)))
        
        elif timeAvgSA.direction == 'z':
            xs = np.linspace(timeAvgSA.cellsize/2, timeAvgSA.meshdims[2]-timeAvgSA.cellsize/2, len(ys))
            plt.xlabel('z (nm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=2, color=color, marker='v', markersize=10)
            ticks = np.linspace(0, timeAvgSA.meshdims[2], int(timeAvgSA.meshdims[2]/10)+1)
            plt.xticks(ticks)
            
        elif timeAvgSA.direction == 'y':
            xs = np.linspace(timeAvgSA.cellsize/2, timeAvgSA.meshdims[1]-timeAvgSA.cellsize/2, len(ys))
            plt.xlabel('y (nm)')
            plt.plot(xs, np.array(ys)/1e11, linewidth=2, color=color, marker='s', markersize=10)
            ticks = np.linspace(0, timeAvgSA.meshdims[1], int(timeAvgSA.meshdims[1]/5)+1)
            plt.xticks(ticks)
        
    plt.ylabel(r'$\mu$')
    
    plt.legend(['Uniaxial system', 'Biaxial system'])    
    
    plt.tight_layout()
        

    if not sim_num:
        plotname = timeAvgSA.plotname()
    else:
        plotname = timeAvgSA.plotname()
        plotname += '_sim' + str(sim_num) + '.png'
    params.make_folder(plotname)
    plt.savefig(plotname, dpi=500)
    
    plt.show()
        
        # plt.close()

def plot_tAvg_SA_underneath(meshdims, cellsize, t, V, damping, MEC, ani, T, type, x_start, x_stop):

    '''
    
    Similar function to the standard, but plots <mxdmdt> at the bottom of the system instead of the top.

    NOTE: This function should be updated to the object-oriented version.
    
    '''

    plt.figure(figsize=(10, 7))

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/underneath/plots/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = type + '/' + modules_folder + ani + '/underneath/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'
    filename_2D = type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'

    if os.path.isfile(filename):
        f = open(filename, 'r')

        lines = f.readlines()
        lines = lines[10:]

        xs = np.linspace(x_start/1000, x_stop/1000, int((x_stop - x_start)/cellsize))

        vals = []

        for line in lines:
            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if ani == 'OOP':
                direction = 2
            vec = line.split('\t')
            all_vals = vec[1:]
            temp = []
            while direction < len(all_vals):
                temp.append(float(all_vals[direction]))
                direction += 3
            vals.append(temp)


        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)

    elif os.path.isfile(filename_2D):
        f = open(filename_2D, 'r')

        lines = f.readlines()
        lines = lines[10:]

        # Turn the raw data into a list of numpy arrays. Every entry in the arrays are floats.
        raw_data = []
        for line in lines:

            # Make a list of all entries and an empty array to fill only the component we want in
            vec = line.strip().split('\t')[1:]
            temp = np.empty(int(len(vec)/3))

            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if ani == 'OOP':
                direction = 2

            # Iterate over all and add only the component we want. Convert to float
            indexer = 0
            while direction < len(vec):
                temp[indexer] = float(vec[direction])
                indexer += 1
                direction += 3

            # Reshape to 2D array and add to the data list
            raw_data.append(temp.reshape(meshdims[2], int(len(temp)/(meshdims[2]))))

        # Now find the time averages for all the data
        tAvg_data = np.zeros_like(raw_data[0])

        for k, matrix in enumerate(raw_data):
            for i, row in enumerate(matrix):
                for j, col in enumerate(row):
                    tAvg_data[i][j] += col
                    if k == len(raw_data)-1:
                        tAvg_data[i][j] /= len(raw_data)


        
        ys = [tAvg_data[0][i] for i in range(len(tAvg_data[0]))]
        xs = np.linspace(x_start, x_stop/1000, len(ys))

    else:
        print("No simulations have been done with these params")
        exit()

    ys2 = [y/ys[0] for y in ys[1:]] # Normalize to the injector

    plt.plot(xs, ys, linewidth=2)
    plt.xlabel(r'$x$ $(\mu m)$')
    plt.ylabel(r'$\mu_x$')
    plt.tight_layout()

    plotname = type + '/' + modules_folder + ani + '/underneath/plots/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.png'
    plt.savefig(plotname, dpi=500)
    # plt.show()

def plot_tAvg_SA_2D(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani):

    '''
    
    Similar function to the standard, but plots <mxdmdt> as a 2D heatmap of the system (along x and z) instead of the top.

    NOTE: This function should be updated to the object-oriented version.
    
    '''

    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/plots/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = ani + '/cache/' + mec_folder + 't_avg/'  + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '.txt'
    f = open(filename, 'r')

    lines = f.readlines()
    lines = lines[10:]

    # Turn the raw data into a list of numpy arrays. Every entry in the arrays are floats.
    raw_data = []
    for line in lines:

        # Make a list of all entries and an empty array to fill only the component we want in
        vec = line.strip().split('\t')[1:]
        temp = np.empty(int(len(vec)/3))

         # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
        direction = 0
        if ani == 'OOP':
            direction = 2

        # Iterate over all and add only the component we want. Convert to float
        indexer = 0
        while direction < len(vec):
            temp[indexer] = float(vec[direction])
            indexer += 1
            direction += 3

        # Reshape to 2D array and add to the data list
        raw_data.append(temp.reshape(meshdims[2], int(len(temp)/(meshdims[2]))))

    # Now find the time averages for all the data
    tAvg_data = np.zeros_like(raw_data[0]) # Haven't tried this before, but should work, right?

    for k, matrix in enumerate(raw_data):
        for i, row in enumerate(matrix):
            for j, col in enumerate(row):
                tAvg_data[i][j] += col
                if k == len(raw_data)-1:
                    tAvg_data[i][j] /= len(raw_data)
    plt.figure(figsize=(12, 4))
    plt.imshow((np.flip(tAvg_data)), extent=[0, (x_stop-x_start)/1000, 0, meshdims[2]], aspect='auto', interpolation='bilinear', cmap='cividis')
    plt.colorbar(label=r'$\mu$')
    plt.xlabel(r'x ($\mu$m)')
    plt.ylabel(r'z (nm)')
    plt.tight_layout()

    plotname = ani + '/plots/' + mec_folder + 't_avg/' +  str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '_t' + str(t) + 'ps.png'
    plt.savefig(plotname, dpi=500)
    plt.show()

def plot_tAvg_SA_2D_y(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani):

    '''
    
    Similar function to the standard 2D, but plots a 2D heatmap along x and y direction.

    NOTE: This function should be updated to the object-oriented version.
    
    '''

    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/plots/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = ani + '/cache/' + mec_folder + 't_avg/'  + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/ydir_2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '.txt'
    f = open(filename, 'r')

    lines = f.readlines()
    lines = lines[10:]

    # Turn the raw data into a list of numpy arrays. Every entry in the arrays are floats.
    raw_data = []
    for line in lines:

        # Make a list of all entries and an empty array to fill only the component we want in
        vec = line.strip().split('\t')[1:]
        temp = np.empty(int(len(vec)/3))

         # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
        direction = 0
        if ani == 'OOP':
            direction = 2

        # Iterate over all and add only the component we want. Convert to float
        indexer = 0
        while direction < len(vec):
            temp[indexer] = float(vec[direction])
            indexer += 1
            direction += 3

        # Reshape to 2D array and add to the data list
        raw_data.append(temp.reshape(meshdims[1], int(len(temp)/(meshdims[1]))))

    # Now find the time averages for all the data
    tAvg_data = np.zeros_like(raw_data[0]) # Haven't tried this before, but should work, right?

    for k, matrix in enumerate(raw_data):
        for i, row in enumerate(matrix):
            for j, col in enumerate(row):
                tAvg_data[i][j] += col
                if k == len(raw_data)-1:
                    tAvg_data[i][j] /= len(raw_data)
    plt.figure(figsize=(12, 4))
    plt.imshow((np.flip(tAvg_data)), extent=[0, (x_stop-x_start)/1000, 0, meshdims[1]], aspect='auto', interpolation='bilinear', cmap='cividis')
    plt.colorbar(label=r'$\mu$')
    plt.xlabel(r'x ($\mu$m)')
    plt.ylabel(r'y (nm)')
    plt.tight_layout()

    plotname = ani + '/plots/' + mec_folder + 't_avg/' +  str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/ydir_2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '_t' + str(t) + 'ps.png'
    plt.savefig(plotname, dpi=500)
    plt.show()

def plot_tAvg_SA_2D_subplots(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani):
    
    '''
    
    Plots several 2D magnon transport systems in the same figure.

    NOTE: This function should be updated to the object-oriented version.
    
    '''
    
    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/plots/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = ani + '/cache/' + mec_folder + 't_avg/'  + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '.txt'
    f = open(filename, 'r')

    lines = f.readlines()
    lines = lines[10:]

    # Turn the raw data into a list of numpy arrays. Every entry in the arrays are floats.
    raw_data = []
    for line in lines:

        # Make a list of all entries and an empty array to fill only the component we want in
        vec = line.strip().split('\t')[1:]
        temp = np.empty(int(len(vec)/3))

         # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
        direction = 0
        if ani == 'OOP':
            direction = 2

        # Iterate over all and add only the component we want. Convert to float
        indexer = 0
        while direction < len(vec):
            temp[indexer] = float(vec[direction])
            indexer += 1
            direction += 3

        # Reshape to 2D array and add to the data list
        raw_data.append(temp.reshape(meshdims[1], int(len(temp)/(meshdims[1]))))

    # Now find the time averages for all the data
    tAvg_data = np.zeros_like(raw_data[0]) # Haven't tried this before, but should work, right?

    for k, matrix in enumerate(raw_data):
        for i, row in enumerate(matrix):
            for j, col in enumerate(row):
                tAvg_data[i][j] += col
                if k == len(raw_data)-1:
                    tAvg_data[i][j] /= len(raw_data)
    plt.figure(figsize=(12, 4))
    plt.imshow((np.flip(tAvg_data)), extent=[0, (x_stop-x_start)/1000, 0, meshdims[1]], aspect='auto', interpolation='bilinear', cmap='cividis')
    plt.colorbar()
    plt.xlabel(r'x ($\mu$m)')
    plt.ylabel(r'y (nm)')
    plt.tight_layout()

    plotname = ani + '/plots/' + mec_folder + 't_avg/' +  str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '_t' + str(t) + 'ps.png'
    plt.savefig(plotname, dpi=500)
    plt.show()

def plot_tAvg_comparison_with_reference(plots, legends, savename, ani):

    '''
    
    Function for a very specific plot.

    NOTE: This function should be updated to the object-oriented version.
    
    '''

    indexer = 0
    plt.figure(figsize=(13,8))


    N = len(plots)
    colors = plt.cm.viridis(np.linspace(0,1,N+1))
    # colors = ['tab:blue', 'tab:red']

    f1 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x5/tAvg_damping0.0004_V-0.06_0.3K.txt'
    f2 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x10/tAvg_damping0.0004_V-0.12_0.3K.txt'
    f3 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x15/tAvg_damping0.0004_V-0.18_0.3K.txt'
    f4 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x20/tAvg_damping0.0004_V-0.25_0.3K.txt'
    f5 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x25/tAvg_damping0.0004_V-0.32_0.3K.txt'
    f6 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x30/tAvg_damping0.0004_V-0.38_0.3K.txt'
    f7 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x35/tAvg_damping0.0004_V-0.45_0.3K.txt'
    f8 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x40/tAvg_damping0.0004_V-0.52_0.3K.txt'
    f9 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x45/tAvg_damping0.0004_V-0.58_0.3K.txt'

    plts1 = [f1,f2,f3,f4,f5,f6,f7,f8,f9]
    inj_vals = []

    for k, plot in enumerate(plts1):

        f = open(plot, 'r')

        plat = plot.split('/')[-1]
        temps = plat.split('_')
        temp = temps[2]
        V = float(temp[1:])

        lines = f.readlines()
        lines = lines[10:]

        vals = []

        for i in range(len(lines)):
            vec1 = lines[i].split('\t')
            all_vals = vec1[1:]
            ani_int = 0
            if ani == 'OOP':
                ani_int = 2
            temp = []
            while ani_int < len(all_vals):
                temp.append(float(all_vals[ani_int]))
                ani_int += 3
            vals.append(temp)

        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)


        # if len(ys) > 1980:
        #     ys = ys[len(ys)-1980:]

        # ys = [(y) / ys[int(len(ys)/2)]  for y in ys]
        # if k == 0:
        #     ys = ys[int(len(ys)/2):]
        # else:
        #     ys = ys[int(len(ys)/2):]
        #     xs = np.linspace(2.02, 2.5, len(ys))

        # If large temperatures we have to adjust for background noise
        if k > 2:
            background = 4*sum([y for y in ys[3*int(len(ys)/4):]]) / len(ys)
            ys = [y - background for y in ys]

        inj_vals.append(ys[0])

    for k, plot in enumerate(plots):

        f = open(plot, 'r')

        plat = plot.split('/')[-1]
        temps = plat.split('_')
        temp = temps[2]
        V = float(temp[1:])

        lines = f.readlines()
        lines = lines[10:]

        vals = []

        for i in range(len(lines)):
            vec1 = lines[i].split('\t')
            all_vals = vec1[1:]
            ani_int = 0
            if ani == 'OOP':
                ani_int = 2
            temp = []
            while ani_int < len(all_vals):
                temp.append(float(all_vals[ani_int]))
                ani_int += 3
            vals.append(temp)

        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)


        # if len(ys) > 1980:
        #     ys = ys[len(ys)-1980:]

        # ys = [(y) / ys[int(len(ys)/2)]  for y in ys]
        # if k == 0:
        #     ys = ys[int(len(ys)/2):]
        # else:
        #     ys = ys[int(len(ys)/2):]
        #     xs = np.linspace(2.02, 2.5, len(ys))

        # If large temperatures we have to adjust for background noise
        if k > 2:
            background = 4*sum([y for y in ys[3*int(len(ys)/4):]]) / len(ys)
            ys = [y - background for y in ys]

        leg = round(ys[0] * 1e-11, 1)

        ys = [(y )/inj_vals[k] for y in ys]
        xs = np.linspace(2.02, 4, len(ys))
        # xs = np.linspace(2.02, 4, len(ys))

        

        plt.plot(xs, ys, label=legends[k], linewidth=3, color=colors[k])
        # plt.plot(xs, ys, label=str(leg), linewidth=3, color=colors[k])

        indexer += 1

    plt.xlabel('x (μm)')
    plt.ylabel(r'$\mu$ (normalized)')

    # plt.legend(title= r'$\mu$ value underneath injector ($10^{11}$)')
    plt.legend(title='Thickness (layers)')

    plt.savefig(savename, dpi=500)
    plt.show()

def plot_tAvg_comparison(plots, legends, savename, ani):

    '''
    
    Plots several systems in the same plot for comparison of spin diffusion lengths.
    
    '''

    plt.figure(figsize=(14,7))
    ax = plt.subplot(111)

    N = len(plots)
    colors = plt.cm.viridis(np.linspace(0,1,N+1))
    
    def get_lab_ramp(base_hex, N, lightness_range=(85, 40)):
        # Convert base color to LAB
        base_rgb = np.array(to_rgb(base_hex))[None, :]
        base_lab = cspace_convert(base_rgb, "sRGB1", "CIELab")[0]

        # Vary lightness, keep a and b (chroma) fixed
        L_vals = np.linspace(lightness_range[0], lightness_range[1], N)
        lab_colors = np.array([[L, base_lab[1], base_lab[2]] for L in L_vals])
        rgb_colors = cspace_convert(lab_colors, "CIELab", "sRGB1")

        # Clip to valid RGB
        rgb_colors = np.clip(rgb_colors, 0, 1)
        return [to_hex(c) for c in rgb_colors]
    
    def generate_hue_variants(base_hex, N, hue_shift=0.3, lightness_range=(0.35, 0.7)):
        """Generate N colors around a base hue using HSL variation."""
        base_rgb = np.array(to_rgb(base_hex))
        base_hls = colorsys.rgb_to_hls(*base_rgb)

        hues = np.linspace(base_hls[0] - hue_shift, base_hls[0] + hue_shift, N)
        lights = np.linspace(lightness_range[0], lightness_range[1], N)

        colors = [
            to_hex(colorsys.hls_to_rgb(h, l, base_hls[2]))
            for h, l in zip(hues, lights)
        ]
        return colors
    
    green_color = '#1B9E77'
    orange_color = '#D95F02'

    # Farger: dus grønn → hovedgrønn → blågrønn
    cmap_colors_green = [
    (0.75, 0.93, 0.87),              # dus grønn (ikke for lys)
    mcolors.to_rgb(green_color),     # hovedfarge
    (0.0, 0.55, 0.65)               # dyp blågrønn
    ]
    
    cmap_colors_orange = [
        (1.0, 0.60, 0.35),              # mettet dus rød-oransje
        mcolors.to_rgb(orange_color),   # hovedoransje
        (0.4, 0.15, 0.05)              # mørk rødbrun oransje
    ]
    
    custom_cmap = mcolors.LinearSegmentedColormap.from_list("CustomGn", cmap_colors_green)
    
    # colors = [custom_cmap(i / (N - 1)) for i in range(N)]
    # colors = get_lab_ramp(green_color, N)
    # colors = generate_hue_variants(green_color, N)
    
    
    x0 = 0.5
    y0 = 0.3
    dx = 0.2
    dy = 0.33
    
    w = 0.2
    h = 0.33
    
    x1, x2, y1, y2 = 3.1 , 3.43, 0.0 , 0.33
    axins = ax.inset_axes(
        [x0, y0, w*2, h*2],
        xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
    ax.indicate_inset_zoom(axins, edgecolor='black')
    axins.set_yticks([])


    for k, plot in enumerate(plots):

        f = open(plot, 'r')

        plat = plot.split('/')[-1]
        temps = plat.split('_')
        temp = temps[2]
        V = float(temp[1:])

        lines = f.readlines()
        lines = lines[10:]

        vals = []

        for i in range(len(lines)):
            vec1 = lines[i].split('\t')
            all_vals = vec1[1:]
            ani_int = 0
            if ani == 'OOP':
                ani_int = 2
            temp = []
            while ani_int < len(all_vals):
                temp.append(float(all_vals[ani_int]))
                ani_int += 3
            vals.append(temp)

        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)


        # if len(ys) > 1980:
        #     ys = ys[len(ys)-1980:]

        # ys = [(y) / ys[int(len(ys)/2)]  for y in ys]
        # if k == 0:
        #     ys = ys[int(len(ys)/2):]
        # else:
        #     ys = ys[int(len(ys)/2):]
        #     xs = np.linspace(2.02, 2.5, len(ys))

        # Adjust for background noise
        if k > 2:
            background = 4*sum([y for y in ys[3*int(len(ys)/4):]]) / len(ys)
            ys = [y - background for y in ys]

        leg = round(ys[0] * 1e-11, 1)

        ys = [(y )/ys[0] for y in ys]
        xs = np.linspace(3.02, 6, len(ys))

        ax.plot(xs, ys, label=legends[k], linewidth=3, color=colors[k])
        # plt.plot(xs, ys, label=str(leg), linewidth=3, color=colors[k])
        
        
        if k == 4:
            x1, x2, y1, y2 = 3.1 , 3.43, 0.0 , 0.33
            axins = ax.inset_axes(
                [x0, y0+dy, w, h],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            
        elif k == 5:
            x1, x2, y1, y2 = 3.1 , 3.43, 0.0 , 0.33
            axins = ax.inset_axes(
                [x0+dx, y0+dy, w, h],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            
        elif k == 6:
            x1, x2, y1, y2 = 3.1 , 3.43, 0.0 , 0.33
            axins = ax.inset_axes(
                [x0, y0, w, h],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            
        elif k == 7:
            x1, x2, y1, y2 = 3.1 , 3.43, 0.0 , 0.33
            axins = ax.inset_axes(
                [x0+dx, y0, w, h],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
        

        if k + factor > 3:
            axins.plot(xs, ys, color=colors[k], linewidth=3)
            
            axins.set_yticks([])

    plt.xlabel('x (μm)')
    plt.ylabel(r'$\mu$/$\mu_0$')
    
    

    

    # plt.legend(title= r'$\mu$ value underneath injector ($10^{11}$)')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Thickness (layers)')

    plt.tight_layout()
    plt.savefig(savename, dpi=500)
    plt.show()
    
def plot_tAvg_comparison_thick(plots, legends, savename, ani):

    '''
    
    Plots several systems in the same plot for comparison of spin diffusion lengths.
    Used for the very thick systems
    
    '''

    plt.figure(figsize=(12,7))
    ax = plt.subplot(111)

    N = len(plots)
    # colors = plt.cm.viridis(np.linspace(0,1,N+1))
    
    def get_lab_ramp(base_hex, N, lightness_range=(85, 40)):
        # Convert base color to LAB
        base_rgb = np.array(to_rgb(base_hex))[None, :]
        base_lab = cspace_convert(base_rgb, "sRGB1", "CIELab")[0]

        # Vary lightness, keep a and b (chroma) fixed
        L_vals = np.linspace(lightness_range[0], lightness_range[1], N)
        lab_colors = np.array([[L, base_lab[1], base_lab[2]] for L in L_vals])
        rgb_colors = cspace_convert(lab_colors, "CIELab", "sRGB1")

        # Clip to valid RGB
        rgb_colors = np.clip(rgb_colors, 0, 1)
        return [to_hex(c) for c in rgb_colors]
    
    def generate_hue_variants(base_hex, N, hue_shift=0.3, lightness_range=(0.35, 0.7)):
        """Generate N colors around a base hue using HSL variation."""
        base_rgb = np.array(to_rgb(base_hex))
        base_hls = colorsys.rgb_to_hls(*base_rgb)

        hues = np.linspace(base_hls[0] - hue_shift, base_hls[0] + hue_shift, N)
        lights = np.linspace(lightness_range[0], lightness_range[1], N)

        colors = [
            to_hex(colorsys.hls_to_rgb(h, l, base_hls[2]))
            for h, l in zip(hues, lights)
        ]
        return colors
    
    green_color = '#1B9E77'
    orange_color = '#D95F02'

    # Farger: dus grønn → hovedgrønn → blågrønn
    cmap_colors_green = [
    (0.75, 0.93, 0.87),              # dus grønn (ikke for lys)
    mcolors.to_rgb(green_color),     # hovedfarge
    (0.0, 0.55, 0.65)               # dyp blågrønn
    ]
    
    cmap_colors_orange = [
        (1.0, 0.60, 0.35),              # mettet dus rød-oransje
        mcolors.to_rgb(orange_color),   # hovedoransje
        (0.4, 0.15, 0.05)              # mørk rødbrun oransje
    ]
    
    custom_cmap = mcolors.LinearSegmentedColormap.from_list("CustomGn", cmap_colors_green)
    
    colors = [custom_cmap(i / (N - 1)) for i in range(N)]
    # colors = get_lab_ramp(green_color, N)
    # colors = generate_hue_variants(green_color, N)

    for k, plot in enumerate(plots):

        f = open(plot, 'r')

        plat = plot.split('/')[-1]
        temps = plat.split('_')
        temp = temps[2]
        V = float(temp[1:])

        lines = f.readlines()
        lines = lines[10:]

        vals = []

        for i in range(len(lines)):
            vec1 = lines[i].split('\t')
            all_vals = vec1[1:]
            ani_int = 0
            if ani == 'OOP':
                ani_int = 2
            temp = []
            while ani_int < len(all_vals):
                temp.append(float(all_vals[ani_int]))
                ani_int += 3
            vals.append(temp)

        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)


        # if len(ys) > 1980:
        #     ys = ys[len(ys)-1980:]

        # ys = [(y) / ys[int(len(ys)/2)]  for y in ys]
        # if k == 0:
        #     ys = ys[int(len(ys)/2):]
        # else:
        #     ys = ys[int(len(ys)/2):]
        #     xs = np.linspace(2.02, 2.5, len(ys))

        # Adjust for background noise
        if k > 2:
            background = 4*sum([y for y in ys[3*int(len(ys)/4):]]) / len(ys)
            ys = [y - background for y in ys]

        leg = round(ys[0] * 1e-11, 1)

        ys = [(y )/ys[0] for y in ys]
        xs = np.linspace(3.02, 6, len(ys))

        ax.plot(xs, ys, label=legends[k], linewidth=3, color=colors[k])
        # plt.plot(xs, ys, label=str(leg), linewidth=3, color=colors[k])
        

    plt.xlabel('x (μm)')
    plt.ylabel(r'$\mu$/$\mu_0$')
    
    

    

    # plt.legend(title= r'$\mu$ value underneath injector ($10^{11}$)')
    ax.legend(title='Thickness (layers)')

    plt.tight_layout()
    plt.savefig(savename, dpi=500)
    plt.show()


#### MAGNON DISPERSION PLOTS ####

def plot_analytical_AFM_uni(ax):
    
    '''
    
    Plot the analytical dispersion relation for uniaxial anisotropy on a given axis. 
    Used for plotting on top of numerically found dispersions
    
    '''
    
    # Define the constants
    Ah = -460e3 # J/m^3
    K_easy = 21 # J/m^3
    A = 76e-15 # J/m
    a = 5e-9 # m
    hbar = 1.0546e-34 # m^2 kg/s
    gamma_e_over_2pi = 2.802e10 # s^-1 T^-1
    Ms = 2.1e3 # A/m
    C = gamma_e_over_2pi/(Ms)
    
    # Define the dispersion relation
    def f(q_a):
        q = q_a / 5e-9
        return np.sqrt(C**2 * 16 * np.abs(Ah)*K_easy + 16*A*np.abs(Ah)*q**2*C**2)/1e12
    
    x = np.linspace(-0.75, 0.75, 1000)
    
    ax.plot(x, f(x), color='red', linestyle='dashed', linewidth=2.5)
    
def plot_analytical_AFM_bi(axes):
    
    '''
    
    Plot the analytical dispersion relation for biaxial anisotropy on two given axes. 
    Used for plotting on top of numerically found dispersions
    
    '''
    
    # Define the constants
    Ah = 460e3 # J/m^3
    K_hard = 21
    K_easy = 21e-3 # J/m^3
    A = 76e-15 # J/m
    a = 5e-9 # m
    hbar = 1.0546e-34 # m^2 kg/s
    gamma_e_over_2pi = 2.802e10 # s^-1 T^-1
    Ms = 2.1e3 # A/m
    C = gamma_e_over_2pi/(Ms)
    
    # Define the dispersion relations
    def high(q_a):
        q = q_a / 5e-9
        return C*np.sqrt(16*A*Ah*q**2 + 16*Ah*(K_easy+K_hard))/1e12
    
    def low(q_a):
        q = q_a /5e-9
        return C*np.sqrt(16*A*Ah*q**2 + 16*Ah*K_easy)/1e12
    
    x = np.linspace(-0.75, 0.75, 1000)
    
    axes[0].plot(x, low(x), color='red', linestyle = 'dashed', linewidth=2.5)
    axes[1].plot(x, high(x), color='red', linestyle = 'dashed', linewidth=2.5)
    
def plot_analytical_FM_uni(ax):
    
    '''
    
    Plot the FM analytical dispersion relation for uniaxial anisotropy on a given axis. 
    Used for plotting on top of numerically found dispersions
    
    '''
    
    # Define the constants
    Ah = -460e3 # J/m^3
    K_easy = 2100 # J/m^3
    A = 76e-15 # J/m
    a = 5e-9 # m
    hbar = 1.0546e-34 # m^2 kg/s
    gamma_e_over_2pi = 2.802e10 # s^-1 T^-1
    Ms = 2.1e3 # A/m
    C = gamma_e_over_2pi/(Ms)
    
    # Define the dispersion relation
    def f(q_a):
        q = q_a / 5e-9
        return C * (2*A*q**2 + 2*K_easy) / 1e9
    
    x = np.linspace(-0.75, 0.75, 1000)
    
    ax.plot(x, f(x), color='red', linestyle='dashed', linewidth=2.5)

def plot_magnon_dispersion(magnonDispersion, zoom, clim_max = 1000, analytical=False, sim_num=False):

    '''
    
    Plot the magnon dispersion of a given system, including a zoom-in on an area if given
    In the case of easy-plane plot the two bands in separate subplots (NOTE: Maybe change?)
    Can plot the analytical solutions on top if analytical == True
    
    '''


    if magnonDispersion.type == 'AFM':
        time_step = 0.1e-12
        ylabel = 'f (THz)'
        divisor = 1e12
    elif magnonDispersion.type == 'FM':
        time_step = 1e-12
        ylabel = 'f (GHz)'
        divisor = 1e9
    else:
        raise Exception('Choose type!')
    
    
    output_files = magnonDispersion.cachename()
    if not sim_num:
        output_files = magnonDispersion.cachename()
    else:
       raw_output_files = magnonDispersion.cachename()
       output_files = []
       for output_file in raw_output_files:
            temp = output_file[:-4]
            temp += '_sim' + str(sim_num) + '.txt' 
            output_files.append(temp) 
            
            
    if not sim_num:
        savename = magnonDispersion.plotname()
    else:
        savename = magnonDispersion.plotname()[:-4]
        savename += '_sim' + str(sim_num) + '.png'
    
    params.make_folder(savename)

    if magnonDispersion.hard_axis:
        fig1,ax1 = plt.subplots(1,2)
        fig1.set_figheight(6)
        fig1.set_figwidth(14)
    else:
        fig1,ax1 = plt.subplots()

    titles = ['y-component', 'z-component']
    
    axes_list = []

    for i, output_file in enumerate(output_files):

        pos_time = np.loadtxt(output_file)

        fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

        freq_len = len(fourier_data)
        k_len = len(fourier_data[0])
        freq = np.fft.fftfreq(freq_len, time_step)
        kvector = np.fft.fftfreq(k_len, 5e-9)

        k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
        f_min = np.abs(freq[0])
        f_max = np.abs(freq[int(0.5 * len(freq))])/divisor # to make it THz
        f_points = int(0.5 * freq_len)

        final_result = [fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]

        label = r'$q_{' + magnonDispersion.axis + r'}$a' 

        if magnonDispersion.hard_axis:
            ax1[i].imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))
            ax1[i].set_xlabel(label)
            ax1[i].set_ylabel(ylabel)

            if zoom:
                x1, x2, y1, y2 = -0.3, 0.3, 0, 0.4
                axins = ax1[i].inset_axes(
                    [0.4, 0.52, 0.6, 0.4],
                    xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
                if i == 0:
                    axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,20000))
                else:
                    axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower')
                
                ax1[i].indicate_inset_zoom(axins, edgecolor='white', alpha=1)
            
                for spine in axins.spines.values():
                    spine.set_edgecolor('white')
                    
                axins.tick_params(colors='white')
                axins.set_yticks([0.0, 0.2, 0.4])
                
                ax1[i].title.set_text(titles[i])
                
                if analytical:
                    axes_list.append(axins)

        else:
            ax1.imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))
            ax1.set_xlabel(label)
            ax1.set_ylabel(ylabel)

            if zoom:
                x1, x2, y1, y2 = -0.3, 0.3, 30, 0.4
                axins = ax1.inset_axes(
                    [0.4, 0.52, 0.6, 0.4],
                    xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
                axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,5000))
                axins.set_yticks([0.0, 0.2, 0.4])
                # axins.set_yticks([0.0, 0.3, 0.6])

                ax1.indicate_inset_zoom(axins, edgecolor='white', alpha=1)
            
                for spine in axins.spines.values():
                    spine.set_edgecolor('white')
                    
                axins.tick_params(colors='white')
                
                if analytical and magnonDispersion.type == 'AFM':
                    plot_analytical_AFM_uni(axins)
            if analytical and magnonDispersion.type == 'FM':
                plot_analytical_FM_uni(ax1)
                plt.xlim(-2,2)

        # ax1.set_ylim(0, 0.1)
    
    if analytical and magnonDispersion.hard_axis:
        plot_analytical_AFM_bi(axes_list)
    
    plt.tight_layout()

    plt.savefig(savename, dpi=500)

    plt.show()

def plot_magnon_dispersion_separate(magnonDispersion, clim_max = 1000, analytical=False):

    '''
    
    Plot the magnon dispersion of a given system. In the case of easy-plane, the y and z components
    are plotted next to each other. 

    '''

    if magnonDispersion.type == 'AFM':
        time_step = 0.1e-12
        ylabel = 'f (THz)'
        divisor = 1e12
    elif magnonDispersion.type == 'FM':
        time_step = 1e-12
        ylabel = 'f (GHz)'
        divisor = 1e9
        
    axins_list = []

    if magnonDispersion.triple:
        output_files = magnonDispersion.cachename()
        savename = magnonDispersion.plotname()
        
        fig, ax = plt.subplots(1,3)

        fig.set_figheight(5)
        fig.set_figwidth(15)
        yvals =[5,25,45]
        
        for i, output_filex in enumerate(output_files):
            
            pos_time = np.loadtxt(output_filex)

            fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

            freq_len = len(fourier_data)
            k_len = len(fourier_data[0])
            freq = np.fft.fftfreq(freq_len, time_step)
            kvector = np.fft.fftfreq(k_len, 5e-9)

            k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
            f_min = np.abs(freq[0])
            f_max = np.abs(freq[int(0.5 * len(freq))])/divisor # THz for FM and GHz for FM
            f_points = int(0.5 * freq_len)

            result = [fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]

            ax[i].imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))

            label = r'$q_{' + magnonDispersion.axis + r'}$a' 

            ax[i].set_xlabel(label)
            ax[i].set_ylabel(ylabel)
            # ax1.set_ylim(0, 0.1)
            title = r'y = ' + str(yvals[i]*1e-3) + ' $\mu$m'
            ax[i].title.set_text(title)

        folder_name = '/'.join(savename.split('/')[:-2])
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        plt.tight_layout()

    elif magnonDispersion.hard_axis:
        
        titles = ['y component', 'z component']
        
        output_files = magnonDispersion.cachename()
        savename = magnonDispersion.plotname()
        params.make_folder(savename)
        
        fig, ax = plt.subplots(1,2)

        fig.set_figheight(5)
        fig.set_figwidth(12)
        yvals =[5,25,45]
        
        for i, output_filex in enumerate(output_files):
            
            pos_time = np.loadtxt(output_filex)

            fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

            freq_len = len(fourier_data)
            k_len = len(fourier_data[0])
            freq = np.fft.fftfreq(freq_len, time_step)
            kvector = np.fft.fftfreq(k_len, 5e-9)

            k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
            f_min = np.abs(freq[0])
            f_max = np.abs(freq[int(0.5 * len(freq))])/divisor # THz for FM and GHz for FM
            f_points = int(0.5 * freq_len)

            result = [fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]

            ax[i].imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))

            label = r'$q_{' + magnonDispersion.axis + r'}$a'

            ax[i].set_xlabel(label)
            ax[i].set_ylabel(ylabel)
            ax[i].set_title(titles[i])

            x1, x2, y1, y2 = -0.3, 0.3, 0, 0.4
            axins = ax[i].inset_axes(
                [0.4, 0.52, 0.6, 0.4],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            axins.imshow(result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,15000))
            axins.set_yticks([0.0, 0.2, 0.4])

            ax[i].indicate_inset_zoom(axins, edgecolor='white', alpha=1)
        
            for spine in axins.spines.values():
                spine.set_edgecolor('white')
                
            axins.tick_params(colors='white')
            
            if analytical:
                axins_list.append(axins)

        if analytical and magnonDispersion.hard_axis:
            plot_analytical_AFM_bi(axins_list)
        
        plt.tight_layout()
        
    
    
    else:
        output_file = magnonDispersion.cachename()[0]
        savename = magnonDispersion.plotname()

        pos_time = np.loadtxt(output_file)

        fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

        freq_len = len(fourier_data)
        k_len = len(fourier_data[0])
        freq = np.fft.fftfreq(freq_len, time_step)
        kvector = np.fft.fftfreq(k_len, 5e-9)

        k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
        f_min = np.abs(freq[0])
        f_max = np.abs(freq[int(0.5 * len(freq))])/divisor # to make it THz
        f_points = int(0.5 * freq_len)

        result = [fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]

        fig1,ax1 = plt.subplots()

        ax1.imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))

        label = r'$q_{' + magnonDispersion.axis + r'}$a'

        ax1.set_xlabel(label)
        ax1.set_ylabel(ylabel)
        # ax1.set_ylim(0, 0.1)

        plt.tight_layout()

        folder_name = '/'.join(savename.split('/')[:-1])
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

    plt.savefig(savename, dpi=500)

    plt.show()

def plot_magnon_dispersion_triple(magnonDispersion, zoom, clim_max = 1000):

    '''
    
    Plot 3 magnon dispersion plots horizontally in the same figure, including a zoom-in on an area.
    Mostly used for comparison between magnon dispersion at different y-values in the system.
    
    '''

    # FM and AFM range in different frequencies
    if magnonDispersion.type == 'AFM':
        time_step = 0.1e-12
        ylabel = 'f (THz)'
        divisor = 1e12
    elif magnonDispersion.type == 'FM':
        time_step = 1e-12
        ylabel = 'f (GHz)'
        divisor = 1e9

    output_files = magnonDispersion.cachename()
    savename = magnonDispersion.plotname()
    params.make_folder(savename)
    
    fig, ax = plt.subplots(1,3, sharey=True)

    fig.set_figheight(5)
    fig.set_figwidth(16)
    yvals =[5,25,45]

    # If hard axis then two consecutive files are y and x components
    # Then they need to be superimposed on the same plot
    loops = 1
    if magnonDispersion.hard_axis:
        loops = 2
    
    # Loop over the output files 
    for i in range(0, int(len(output_files)), loops):

        results = []

        # Add both y and z component to the results list if hard axis
        # If easy axis then just the given component
        for j in range(loops):
            pos_time = np.loadtxt(output_files[i+j])

            fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

            freq_len = len(fourier_data)
            k_len = len(fourier_data[0])
            freq = np.fft.fftfreq(freq_len, time_step)
            kvector = np.fft.fftfreq(k_len, 5e-9)

            k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
            f_min = np.abs(freq[0])
            f_max = np.abs(freq[int(0.5 * len(freq))])/divisor # THz for FM and GHz for FM
            f_points = int(0.5 * freq_len)

            result = np.array([fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)])
            results.append(result)

        if loops == 2:
            final_result = results[0] + results[1]
        else:
            final_result = results[0]

        index = int(i/loops)

        ax[index].imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))

        label = r'$q_{' + magnonDispersion.axis + r'}$a' 

        ax[index].set_xlabel(label)
        if index == 0:
            ax[index].set_ylabel(ylabel)
        title = 'y = ' + str(yvals[index]) + 'nm'
        ax[index].title.set_text(title)

        if zoom:
            x1, x2, y1, y2 = -0.35, 0.35, 0, 0.5
            axins = ax[index].inset_axes(
                [0.35, 0.45, 0.7, 0.5],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,5000))
            axins.set_yticks([0.0, 0.2, 0.4])

            ax[index].indicate_inset_zoom(axins, edgecolor='white', alpha=1)
        
            for spine in axins.spines.values():
                spine.set_edgecolor('white')
                
            axins.tick_params(colors='white')

    plt.tight_layout()

    plt.savefig(savename, dpi=500)

    plt.show()

def plot_magnon_dispersion_double(magnonDispersion, zoom, clim_max = 1000, sim_num=False):
    
    '''
    
    Plot 2 magnon dispersion plots horizontally in the same figure, including a zoom-in on an area if desired.
    Used in the thesis since y=5 and y=45 results are the same, so only show y=5 and y=25.
    
    '''

    # FM and AFM range in different frequencies
    if magnonDispersion.type == 'AFM':
        time_step = 0.1e-12
        ylabel = 'f (THz)'
        divisor = 1e12
    elif magnonDispersion.type == 'FM':
        time_step = 1e-12
        ylabel = 'f (GHz)'
        divisor = 1e9

    output_files = magnonDispersion.cachename()
    if not sim_num:
        output_files = magnonDispersion.cachename()
    else:
       raw_output_files = magnonDispersion.cachename()
       output_files = []
       for output_file in raw_output_files:
            temp = output_file[:-4]
            temp += '_sim' + str(sim_num) + '.txt' 
            output_files.append(temp) 
            
            
    if not sim_num:
        savename = magnonDispersion.plotname()
    else:
        savename = magnonDispersion.plotname()[:-4]
        savename += '_sim' + str(sim_num) + '.png'
    params.make_folder(savename)
    
    fig, ax = plt.subplots(1,2, sharey=True)

    fig.set_figheight(5)
    fig.set_figwidth(12)
    yvals =[5,25]

    # If hard axis then two consecutive files are y and x components
    # Then they need to be superimposed on the same plot
    loops = 1
    if magnonDispersion.hard_axis:
        loops = 2
    
    # Loop over the output files 
    for i in range(0, int(len(output_files))-loops, loops):

        results = []

        # Add both y and z component to the results list if hard axis
        # If easy axis then just the given component
        for j in range(loops):
            pos_time = np.loadtxt(output_files[i+j])

            fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

            freq_len = len(fourier_data)
            k_len = len(fourier_data[0])
            freq = np.fft.fftfreq(freq_len, time_step)
            kvector = np.fft.fftfreq(k_len, 5e-9)

            k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
            f_min = np.abs(freq[0])
            f_max = np.abs(freq[int(0.5 * len(freq))])/divisor # THz for FM and GHz for FM
            f_points = int(0.5 * freq_len)

            result = np.array([fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)])
            results.append(result)

        if loops == 2:
            final_result = results[0] + results[1]
        else:
            final_result = results[0]

        index = int(i/loops)

        ax[index].imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))

        label = r'$q_{' + magnonDispersion.axis + r'}$a' 

        ax[index].set_xlabel(label)
        if index == 0:
            ax[index].set_ylabel(ylabel)
            title = 'Edge'
        else:
            title = 'Middle'
        ax[index].title.set_text(title)

        if zoom:
            x1, x2, y1, y2 = -0.8, 0.8, 0.7, 1.5
            axins = ax[index].inset_axes(
                [0.35, 0.45, 0.7, 0.5],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,5000))
            # axins.set_yticks([0.0, 0.2, 0.4])

            ax[index].indicate_inset_zoom(axins, edgecolor='white', alpha=1)
        
            for spine in axins.spines.values():
                spine.set_edgecolor('white')
                
            axins.tick_params(colors='white')

    plt.tight_layout()

    plt.savefig(savename, dpi=500)

    plt.show()
    
def plot_magnon_dispersion_together(magnonDispersion, clim_max = 1000):

    '''

    Plot both y and z component of a magnon dispersion in the same plot. 
    
    '''

    if magnonDispersion.type == 'AFM':
        time_step = 0.1e-12
        ylabel = 'f (THz)'
        divisor = 1e12
    elif magnonDispersion.type == 'FM':
        time_step = 1e-12
        ylabel = 'f (GHz)'
        divisor = 1e9
    else:
        raise Exception('Choose type!')
    
    
    output_files = magnonDispersion.cachename()
    savename = magnonDispersion.plotname()
    params.make_folder(savename)

    fig1,ax1 = plt.subplots()
    fig1.set_figheight(6)
    fig1.set_figwidth(8)

    results = []

    for i, output_file in enumerate(output_files):

        pos_time = np.loadtxt(output_file)

        fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

        freq_len = len(fourier_data)
        k_len = len(fourier_data[0])
        freq = np.fft.fftfreq(freq_len, time_step)
        kvector = np.fft.fftfreq(k_len, 5e-9)

        k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
        f_min = np.abs(freq[0])
        f_max = np.abs(freq[int(0.5 * len(freq))])/divisor # to make it THz
        f_points = int(0.5 * freq_len)

        results.append(np.array([fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]))

    final_result = results[0] + results[1]

    label = 'q' + magnonDispersion.axis

    ax1.imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))
    ax1.set_xlabel(label)
    ax1.set_ylabel(ylabel)

    x1, x2, y1, y2 = -0.25, 0.25, 0, 0.5
    axins = ax1.inset_axes(
        [0.5, 0.45, 0.47, 0.47],
        xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
    axins.set_yticks([0.0,0.2,0.4])
    axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', interpolation='bilinear', clim=(0, 15000))

    ax1.indicate_inset_zoom(axins, edgecolor='black')
        

        # ax1.set_ylim(0, 0.1)

    plt.tight_layout()

    plt.savefig(savename, dpi=500)

    plt.show()

def plot_dispersion_field_comparisons(plots, c_max1, c_max2):

    '''
    
    Plots several magnon dispersion relations in the same figure (subplots). 
    This is a specific function for comparison of external magnetic field in the biaxial system. 
    
    '''

    fig, axs = plt.subplots(4, 1)

    fig.set_figheight(20)
    fig.set_figwidth(7)

    annotations = [r'$\mu_0$H = 0.0T', r'$\mu_0$H = 3.0T', r'$\mu_0$H = 5.9T', r'$\mu_0$H = 11.8T']
    annotations2 = [r'$(a)$', r'$(b)$', r'$(c)$', r'$(d)$']

    indexer = 0

    for i in range(8):
        output_file = plots[i]

        time_step = 0.1e-12
        pos_time = np.loadtxt(output_file)

        fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

        freq_len = len(fourier_data)
        k_len = len(fourier_data[0])
        freq = np.fft.fftfreq(freq_len, time_step)
        kvector = np.fft.fftfreq(k_len, 5e-9)

        k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
        f_min = np.abs(freq[0])
        f_max = np.abs(freq[int(0.5 * len(freq))])/1e12 # to make it THz
        f_points = int(0.5 * freq_len)
        
        if i % 2 == 0:
            result = np.array([fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)])
        else:
            result += np.array([fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)])
        
            label = r'$q_{x}$a' 

            x1, x2, y1, y2 = -0.3, 0.3, 0, 0.4
            axins = axs[indexer].inset_axes(
                [0.4, 0.52, 0.6, 0.4],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])

            axs[indexer].imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0,c_max1[indexer]))
            axins.imshow(result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,c_max2[indexer]))
            axs[indexer].indicate_inset_zoom(axins, edgecolor='white', alpha=1)
            axs[indexer].annotate(annotations[indexer], xy=(0.025,0.05  ), xycoords=('axes fraction'), color='white')
            axs[indexer].annotate(annotations2[indexer], xy=(0.9,0.05), xycoords=('axes fraction'), color='white')

            for spine in axins.spines.values():
                spine.set_edgecolor('white')
                
            axins.tick_params(colors='white')

            axs[indexer].set_xlabel(label)
            axs[indexer].set_ylabel('f (THz)')
            
            indexer += 1

    fig.tight_layout()
    
    savename = 'AFM/custom/plots/biaxial_singlelayer_dispersions_magnetic_field_comparison.png'

    plt.savefig(savename, dpi=500)

    # plt.show()

def plot_phonon_dispersion(meshdims, damping, MEC, ani, dir,time_step):

    '''
    
    Plots the phonon-dispersion of a system. Since phonons are only excited by MEC, the results are not reliable.

    NOTE: This function should be updated to the object-oriented version.
    
    '''

    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'


    output_file = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 'dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/dir' + dir + '_phonon_dispersion.txt'

    pos_time = np.loadtxt(output_file)

    fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

    freq_len = len(fourier_data)
    k_len = len(fourier_data[0])
    freq = np.fft.fftfreq(freq_len, time_step)
    kvector = np.fft.fftfreq(k_len, 5e-9)

    k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
    f_min = np.abs(freq[0])
    f_max = np.abs(freq[int(0.5 * len(freq))])/1e12 # to make it THz
    f_points = int(0.5 * freq_len)

    result = [fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]

    fig1,ax1 = plt.subplots()

    ax1.imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0,1200), vmax=20)

    ax1.set_xlabel(r'$q_x$')
    ax1.set_ylabel('f (THz)')
    # ax1.set_ylim(0, 0.1)

    plt.tight_layout()

    folder_name = ani + '/plots/' + mec_folder + 'dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    savename = ani + '/plots/' + mec_folder + 'dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/damping' + str(damping) + 'dir' + dir + '_phonon_dispersion.png'

    plt.savefig(savename, dpi=500)

    plt.show()

def plot_phonon_dispersion_specific(output_file, savename, time_step):


    pos_time = np.loadtxt(output_file)

    fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

    freq_len = len(fourier_data)
    k_len = len(fourier_data[0])
    freq = np.fft.fftfreq(freq_len, time_step)
    kvector = np.fft.fftfreq(k_len, 5e-9)

    k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
    f_min = np.abs(freq[0])
    f_max = np.abs(freq[int(0.5 * len(freq))])/1e12 # to make it THz
    f_points = int(0.5 * freq_len)

    result = [fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]

    fig1,ax1 = plt.subplots()

    ax1.imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0,1000), vmax=20)

    ax1.set_xlabel(r'$q_x$')
    ax1.set_ylabel('f (THz)')
    # ax1.set_ylim(0, 0.1)

    plt.tight_layout()

    plt.savefig(savename, dpi=500)

    plt.show()

def plot_dispersions(plots, savename):

    '''
    
    Plots several magnon dispersion relations in the same figure (subplots). Titles and such must be specified. 
    
    '''

    dim1 = 0
    dim2 = 0

    if len(plots) == 2:
        dim1 = 1
        dim2 = 2

    elif len(plots) == 3:
        dim1 = 1
        dim2 = 3


    fig, axs = plt.subplots(dim1, dim2, sharey=True)

    fig.set_figheight(5)
    fig.set_figwidth(14)

    annotations = [r'$\frac{|V|}{L_z} = 19μV$ $nm^{-1}$', r'$\frac{|V|}{L_z} = 30.0μV$ $nm^{-1}$']
    # titles = ['1 layer', '50 layers', '100 layers', '150 layers', '150 layers \n + damping layer']
    # titles = ['T = 0.3K', 'T = 0.8K', 'T = 3.0K']

    clim_max = [5000, 15000, 25000]

    for i, ax in enumerate(list((axs))):

        output_file = plots[i]

        time_step = 0.1e-12
        pos_time = np.loadtxt(output_file)

        fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

        freq_len = len(fourier_data)
        k_len = len(fourier_data[0])
        freq = np.fft.fftfreq(freq_len, time_step)
        kvector = np.fft.fftfreq(k_len, 5e-9)

        k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
        f_min = np.abs(freq[0])
        f_max = np.abs(freq[int(0.5 * len(freq))])/1e12 # to make it THz
        f_points = int(0.5 * freq_len)

        result = [fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]

        # clim_max = 1000
        # if i == 0:
        #     clim_max = 1000

        label = r'$q_x a$'

        x1, x2, y1, y2 = -0.25, 0.25, 0, 0.5
        axins = ax.inset_axes(
            [0.5, 0.45, 0.47, 0.47],
            xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])

        ax.imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max[i]))
        axins.imshow(result, extent = [-k_max, k_max,f_min, f_max], origin='lower')
        ax.indicate_inset_zoom(axins, edgecolor='black')

        ax.annotate(annotations[i], (0.05, 0.85), xycoords = 'axes fraction', color='white', fontsize=32)

        # if i == 0:
        #     ax.title.set_text('Without MEC')
        # elif i == 1:
        #     ax.title.set_text('With MEC')

        # ax.title.set_text(titles[i])
        # ax.set_xticks([-2,2])
        ax.set_xlabel(label)
        if i == 0:
            ax.set_ylabel('f (THz)')
    # ax1.set_ylim(0, 0.1)

    fig.tight_layout()
    # fig.subplots_adjust(wspace=0.03)
    

    plt.savefig(savename, dpi=500)


    plt.show()

#### MISCELLANEOUS ####

def plot_trajectory(meshdims, damping, MEC, ani, dir):

    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    output_file1 = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 'trajectory/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' +  dir + '_trajectory_M1.txt'
    output_file2 = 'C:/Users/mathimyh/documents/boris data/simulations/boris_fordypningsoppgave/' + ani + '/cache/' + mec_folder + 'trajectory/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' +  dir + '_trajectory_M2.txt'

    f1 = open(output_file1, 'r')

    lines1 = f1.readlines()[:-1]

    ys1 = np.zeros(len(lines1[0].strip().split('\t')))

    for line1 in lines1:
        vals = line1.strip().split('\t')
        for i, val in enumerate(vals):
            ys1[i] += float(val)

    ys1 = [y/len(ys1) for y in ys1]

    f2 = open(output_file2, 'r')

    lines2 = f2.readlines()[:-1]

    ys2 = np.zeros(len(lines2[0].strip().split('\t')))

    for line2 in lines2:
        vals = line2.strip().split('\t')
        for i, val in enumerate(vals):
            ys2[i] += float(val)

    ys2 = [y/len(ys2) for y in ys2]

    ys = []

    for i in range(len(ys1)):
        ys.append(1/2 * (ys1[i] - ys2[i]))

    xs = np.linspace(0, len(ys), len(ys))

    plt.plot(xs, ys)

    folder_name = ani + '/plots/' + mec_folder + 'trajectory/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    savename = ani + '/plots/' + mec_folder + 'trajectory/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/damping' + str(damping) + '_' + dir + '_trajectory.png'
    plt.savefig(savename, dpi=500)
    plt.show()

def plot_critical_T(criticalT):

    plt.figure(figsize=(10,6))

    output_file = criticalT.cachename()

    f = open(output_file, 'r')

    lines = f.readlines()

    xs = []
    ys = []
    

    for i, line in enumerate(lines):
        temp = line[1:-2]
        vals = temp.strip().split(', ')
        T = int(vals[0])
        m = float(vals[1])
        xs.append(T)
        ys.append(m)
        
    ys = [abs(y) for y in ys]
    max_y = max(ys)
    ys = [y/max_y for y in ys]

    plt.plot(xs, ys, linewidth=3)

    savename = criticalT.plotname()
    params.make_folder(savename)
    plt.xlabel(r'Temperature (K)')
    plt.ylabel(r'|${{m}}_{x}$|')
    plt.tight_layout()
    plt.savefig(savename, dpi=500)
    plt.show()

def plot_diffusion_length_all_systems(all_plots, all_ts, savename, ani):
    '''
    
    Plots the spin diffusion lengths of given magnon transport simulations. 
    Used for all systems in the same figure
    ts are the values that is compared between them
    
    '''
    
    
    plt.figure(figsize=(10,7))
    
    markers = ['s', 'v', '.', 's', 'v', '.']
    thicknesses = [10,20,40,10,20,40]
    colors = ['#2874a6', '#3498db', '#85c1e9', '#af601a', '#e67e22', '#f0b27a']

    for m, plots in enumerate(all_plots):
        
        for k, plot in enumerate(plots):

            try:
                f = open(plot, 'r')
            except:
                print('No file for V = ', -all_ts[m][k])
            else:
                plat = plot.split('/')[-1]
                temps = plat.split('_')
                temp = temps[2]

                lines = f.readlines()
                lines = lines[10:]

                vals = []

                for i in range(len(lines)):
                    vec1 = lines[i].split('\t')
                    all_vals = vec1[1:]
                    ani_int = 0
                    if ani == 'OOP':
                        ani_int = 2
                    temp = []
                    while ani_int < len(all_vals):
                        temp.append(float(all_vals[ani_int]))
                        ani_int += 3
                    vals.append(temp)

                ys = []

                for i in range(len(vals[0])):
                    val = 0
                    for j in range(len(vals)):
                        val += float(vals[j][i])
                    val /= len(vals)
                    ys.append(val)


                # When <mxdmdt> flattens out we get NaN and inf after doing the logarithm. 
                # This needs to be removed
                # From looking at the data, it seems to be fine to just find the 
                # slope from the first half of ys, as the end point introduces lots of noise
                try: 
                    min_y = min(ys)
                except ValueError:
                    print('min(ys) is empty for V = ', all_ts[m][k])
                else:
                    # ys =[y + abs(min_y) for y in ys]
                    ys = np.array([np.log(p) for p in ys])
                    ys = ys[np.isfinite(ys)]
                    xs = np.linspace(0, 2, len(ys))
                        
                    # Find the cut-off before fitting the function
                    algo = rpt.Dynp(model="l2").fit(ys)
                    try:
                        result = algo.predict(n_bkps=1)
                    except:
                        print('Too noisy systesm for V = ', all_ts[m][k])
                    else:
                    
                        # If the signal travels almost to the edges the cut-off is not needed
                        # If the cut-off is larger than 1/2 of the distance, this is the case
                        if result[0] > len(ys*0.5):
                        
                            # Then cut off the list and create x-values list
                            ys = ys[:result[0]]
                            xs = xs[:result[0]]
                        
                        # In this case we just remove the very edges
                        else:
                            ys = ys[:int(3*len(ys)/4)]
                            xs = xs[:int(3*len(xs)/4)]

                        def f(x, a, b):
                            return a - b*x   

                        try:
                            params, params_cov = curve_fit(f, xs, ys)
                            if 1/params[1] < 0.6 and 1/params[1] > 0:
                                plt.plot(10e3*all_ts[m][k]/thicknesses[m], 1/params[1], marker=markers[m], markersize = 10, color=colors[m])
                        except TypeError:
                            print('Could not find a curve for this voltage: ', all_ts[m][k])
                        except ValueError:
                            print('ydata was empty for this voltage: ', all_ts[m][k])
                        
                        
    plt.xlabel(r'Voltage/thickness (μV/nm)')
    plt.ylabel(r'$l_d$')
    
    handles = []
    
    aniuni_2_layer = mlines.Line2D([], [], color=colors[0], marker=markers[0], linestyle='None',
                          markersize=10, label='Uniaxial system')
    
    labels = ['2 layer', '4 layer', '8 layer', '2 layer', '4 layer', '8 layer']
    
    for n in range(len(all_plots)):
        handles.append(mlines.Line2D([], [], color=colors[n], marker=markers[n], linestyle='None',
                          markersize=10, label=labels[n]))
    
    
    plt.legend(handles=handles)
    
    plt.tight_layout()
    plt.savefig(savename, dpi=500)
    

    plt.show()
    
def plot_diffusion_length_both_systems(plots1, plots2, ts, savename, ani, thickness):
    
    '''
    
    Plots the spin diffusion lengths of given magnon transport simulations (two systems in the same plot). 
    ts are the values that is compared between them
    
    '''
    
    
    plt.figure(figsize=(10,7))
    
    both = [plots1, plots2]
    markers = ['s', 'v']
    colors = ['#1B9E77', '#D95F02']

    flips_2layer = {0.0875 : [], 0.100 : [4,5,6,7,9], 0.125 : [3,4,6,7,9], 0.1375: [1,2,3,4,5,7,8],
                    0.1625 : [2,3,5,6,7,8], 0.175 : [2,4,6,7,8], 0.200 : [1,2,3,4,5,6,7,9],
                    0.2125 : [1,2,3,4,5,6,7,8], 0.2375 : [1,2,3,4,5,6,7,8,9], 
                    0.250 : [1,2,3,4,5,6,7,8,9], 0.275 : [1,5,6,9], 0.2875 : [2,7,9]}

    ensembles = 9
    
    # Loop over biaxial and uniaxial
    for m, plots in enumerate(both):    
            
        # Loop over voltages
        for k, raw_plot in enumerate(plots):

            results = []
            
            failed = False

            # Loop over ensembles
            for e in range(1,ensembles+1):

                # For now only 4 layers have many ensembles, update when sims are done for the rest
                if ensembles == 1:
                    plot = raw_plot
                else:
                    plot = raw_plot[:-4]
                    plot += '_sim' + str(e) + '.txt' 
                
                try:
                    f = open(plot, 'r')
                except:
                    print('No file found', plot)
                else:
                    plat = plot.split('/')[-1]
                    temps = plat.split('_')
                    temp = temps[2]

                    lines = f.readlines()
                    lines = lines[10:]

                    vals = []

                    for i in range(len(lines)):
                        vec1 = lines[i].split('\t')
                        all_vals = vec1[1:]
                        ani_int = 0
                        if ani == 'OOP':
                            ani_int = 2
                        temp = []
                        while ani_int < len(all_vals):
                            temp.append(float(all_vals[ani_int]))
                            ani_int += 3
                        vals.append(temp)

                    ys = []

                    for i in range(len(vals[0])):
                        val = 0
                        for j in range(len(vals)):
                            val += float(vals[j][i])
                        val /= len(vals)
                        ys.append(val)


                    # We need to find only the data in the exponential decay using a cut-off
                    # Then we fit this data to a function
                    # Either linear for exponential decay, or exponential with harmonic wave for
                    # resonance plots
                    try: 
                        min_y = min(ys)
                    except ValueError:
                        print('min(ys) is empty for V = ', -ts[k])
                        failed = True
                    else:
                        
                        original_ys = copy.deepcopy(ys)
                        
                        # This is the data set we use for curve fitting
                        fit_ys = [y + abs(min_y) for y in ys]
                        fit_ys = np.array([np.log(p) for p in fit_ys])
                        
                        # This is the data set we use for determining the cut-off
                        cutoff_ys = np.array([np.log(p) for p in ys])
                        cutoff_ys  = cutoff_ys[np.isfinite(cutoff_ys)]
                        
                        xs = np.linspace(0, 2.98, len(fit_ys))
                            
                        # Find the cut-off before fitting the function
                        algo = rpt.Dynp(model="l2").fit(cutoff_ys)
                        try:
                            result = algo.predict(n_bkps=1)
                        except:
                            print('Data was too noisy for V = ', -ts[k])
                            failed = True
                        else:

                            # Then cut off the list ands create new x-values list
                            # One for linear and one for exponential function fitting
                            lin_fit_ys = fit_ys[:result[0]]
                            exp_fit_ys = original_ys[:result[0]]
                            fit_xs = xs[:result[0]]

                            # Linear
                            def f(x, a, b):
                                return a - b*x   
                            
                            # Exponential with harmonic wave
                            def g(x, a, b, c, d, e, f):
                                return a*np.exp(-b*x) + c*np.sin(d*x + e)*np.exp(-f*x)
                            
                            # Exponential
                            def h(x, a, b):
                                return a*np.exp(-b*x)
                            
                            
                            if thickness == 40:
                                try:
                                    params, params_cov = curve_fit(g, fit_xs, exp_fit_ys, p0=[1e12, 2, 5e11, 20, 0, 4])
                                except TypeError:
                                    print('Could not find curve for this voltage')
                                
                                # If we get runtime error then there is no harmonic wave, just exponential decay
                                except RuntimeError:
                                    params, params_cov = curve_fit(h, fit_xs, lin_fit_ys, p0=[1e12, 2])
                                    if 1/params[1] < 2:
                                        results.append(1/params[1])
                                else:
                                    results.append(1/params[1])
                            else:
                                params, params_cov = curve_fit(f, fit_xs, lin_fit_ys)
                                if thickness != 10 or (e not in flips_2layer[ts[k]] or m == 0):
                                    results.append(1/params[1])
                                    
            array = np.array([x for x in results if 0 <= x < 1])
            standard_deviation = np.std(array)
            average = np.average(array)
            if average < 1: 
                plt.errorbar([1e3*ts[k]/thickness], average, yerr=standard_deviation, marker=markers[m], markersize=10, color=colors[m], capsize=3)
            
    plt.xlabel(r'Voltage/thickness (μV/nm)')
    plt.ylabel(r'$L_D (μm)$')
    
    
    green_square = mlines.Line2D([], [], color=colors[0], marker=markers[0], linestyle='None',
                          markersize=10, label='Uniaxial system')
    orange_triangle = mlines.Line2D([], [], color=colors[1], marker=markers[1], linestyle='None',
                          markersize=10, label='Biaxial system')
    
    plt.legend(handles=[green_square, orange_triangle])
    
    plt.tight_layout()
    plt.savefig(savename, dpi=500)
    

    plt.show()

def plot_diffusion_length(plots, ts, savename, ani, thickness):
    
    '''
    
    Plots the spin diffusion lengths of given magnon transport simulations. 
    ts are different voltages
    
    NOTE: Needs to be updated
    
    '''

    plt.figure(figsize=(10,7))

    for k, plot in enumerate(plots):

        try:
            f = open(plot, 'r')
        except:
            print('No file for V = ', -ts[k])
        else:
            plat = plot.split('/')[-1]
            temps = plat.split('_')
            temp = temps[2]

            lines = f.readlines()
            lines = lines[10:]

            vals = []

            for i in range(len(lines)):
                vec1 = lines[i].split('\t')
                all_vals = vec1[1:]
                ani_int = 0
                if ani == 'OOP':
                    ani_int = 2
                temp = []
                while ani_int < len(all_vals):
                    temp.append(float(all_vals[ani_int]))
                    ani_int += 3
                vals.append(temp)

            ys = []

            for i in range(len(vals[0])):
                val = 0
                for j in range(len(vals)):
                    val += float(vals[j][i])
                val /= len(vals)
                ys.append(val)

            try: 
                min_y = min(ys)
            except ValueError:
                print('min(ys) is empty for V = ', ts[k])
            else:
                
                # This is the data set we use for curve fitting
                ys1 = [y + abs(min_y) for y in ys]
                ys1 = np.array([np.log(p) for p in ys1])
                
                # This is the data set we use for determining the cut-off
                ys = np.array([np.log(p) for p in ys])
                ys  = ys[np.isfinite(ys)]
                
                xs = np.linspace(0, 2.98, len(ys1))
                
                # Apply to cut-off data
                algo = rpt.Dynp(model="l2").fit(ys)
                result = algo.predict(n_bkps=1)
                
                # Then cut off the list ands create new x-values list
                ys1 = ys1[:result[0]]
                xs = xs[:result[0]]
                

                def f(x, a, b):
                    return a - b*x   

                try:
                    params, params_cov = curve_fit(f, xs, ys1)
                    plt.plot(10e3*ts[k]/thickness, 1/params[1], marker='v', markersize = 10, color='black')
                    # plt.plot(xs, ys1)
                except TypeError:
                    print('Could not find a curve for this voltage: ', ts[k])
                except ValueError:
                    print('ydata was empty for this voltage: ', ts[k])

        
    plt.xlabel(r'Voltage/thickness (μV/nm)')
    plt.ylabel(r'$l_d$')
    # plt.yticks([0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50])

    plt.tight_layout()
    plt.savefig(savename, dpi=500)

    plt.show()

def plot_diffusion_length_across_thicknesses(plots, ts, savename, ani):
    
    '''
    
    Plots the spin diffusion lengths of given magnon transport simulations. 
    ts are the different thicknesses
    
    '''

    plt.figure(figsize=(10,7))
    
    uniax_color = '#1B9E77'
    uniax_marker = 's'
    biax_color = '#D95F02'
    biax_marker = 'v'

    for k, plot in enumerate(plots):

        try:
            f = open(plot, 'r')
        except:
            print('No file for V = ', -ts[k])
        else:
            plat = plot.split('/')[-1]
            temps = plat.split('_')
            temp = temps[2]

            lines = f.readlines()
            lines = lines[10:]

            vals = []

            for i in range(len(lines)):
                vec1 = lines[i].split('\t')
                all_vals = vec1[1:]
                ani_int = 0
                if ani == 'OOP':
                    ani_int = 2
                temp = []
                while ani_int < len(all_vals):
                    temp.append(float(all_vals[ani_int]))
                    ani_int += 3
                vals.append(temp)

            ys = []

            for i in range(len(vals[0])):
                val = 0
                for j in range(len(vals)):
                    val += float(vals[j][i])
                val /= len(vals)
                ys.append(val)


            try: 
                min_y = min(ys)
            except ValueError:
                print('min(ys) is empty for thickness = ', ts[k])
            else:
                
                original_ys = copy.deepcopy(ys)
                    
                # This is the data set we use for curve fitting
                fit_ys = [y + abs(min_y) for y in ys]
                fit_ys = np.array([np.log(p) for p in fit_ys])
                
                # This is the data set we use for determining the cut-off
                cutoff_ys = np.array([np.log(p) for p in ys])
                cutoff_ys  = cutoff_ys[np.isfinite(cutoff_ys)]
                
                xs = np.linspace(0, 2.98, len(fit_ys))
                    
                # Find the cut-off before fitting the function
                algo = rpt.Dynp(model="l2").fit(cutoff_ys)
                try:
                    result = algo.predict(n_bkps=1)
                except:
                    print('Data was too noisy for V = ', -ts[k])
                else:

                    # Then cut off the list ands create new x-values list
                    # One for linear and one for exponential function fitting
                    lin_fit_ys = fit_ys[:result[0]]
                    exp_fit_ys = original_ys[:result[0]]
                    fit_xs = xs[:result[0]]

                    # Linear fit
                    def f(x, a, b):
                        return a - b*x   
                    
                    # Exponential with harmonic wave
                    def g(x, a, b, c, d, e, f):
                        return a*np.exp(-b*x) + c*np.sin(d*x + e)*np.exp(-f*x)
                    
                    # Exponential
                    def h(x, a, b):
                        return a*np.exp(-b*x)
                    
                    if k < 8:
                        color = uniax_color
                        marker = uniax_marker
                    else:
                        color = biax_color
                        marker = biax_marker
                    
                    if (k > 3 and k < 8) or (k > 11):
                        try:
                            params, params_cov = curve_fit(g, fit_xs, exp_fit_ys, p0=[1e12, 2, 5e11, 20, 0, 4])
                        except TypeError:
                            print('Could not find curve for this voltage')
                        
                        # If we get runtime error then there is no harmonic wave, just exponential decay
                        except RuntimeError:
                            params, params_cov = curve_fit(h, fit_xs, lin_fit_ys, p0=[1e12, 2])
                            if 1/params[1] < 2:
                                plt.plot(ts[k]/5, 1/params[1], marker=marker, markersize = 10, color=color)
                        else:
                            print(1/params[5])
                            plt.plot(ts[k]/5, 1/params[1], marker=marker, markersize = 10, color=color)
                            
                    else:
                        params, params_cov = curve_fit(f, fit_xs, lin_fit_ys)
                        plt.plot(ts[k]/5, 1/params[1], marker=marker, markersize = 10, color=color)

                    # try:
                    #     params, params_cov = curve_fit(f, xs, ys1)
                    #     if k < 8:
                    #         color = uniax_color
                    #         marker = uniax_marker
                    #     else:
                    #         color = biax_color
                    #         marker = biax_marker
                    #     plt.plot(ts[k]/5, 1/params[1], marker=marker, markersize = 10, color=color)
                    #     # plt.plot(xs, ys1)
                    # except TypeError:
                    #     print('Could not find a curve for this thickness: ', ts[k])
                    # except ValueError:
                    #     print('ydata was empty for this thickness: ', ts[k])

        
    plt.xlabel(r'Thickness (layers)')
    plt.ylabel(r'$L_D (μm)$')
    # plt.yticks([0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50])
    plt.xticks([1,2,3,4,5,6,7,8])
    
    green_square = mlines.Line2D([], [], color=uniax_color, marker=uniax_marker, linestyle='None',
                          markersize=10, label='Uniaxial system')
    orange_triangle = mlines.Line2D([], [], color=biax_color, marker=biax_marker, linestyle='None',
                          markersize=10, label='Biaxial system')
    
    plt.legend(handles=[green_square, orange_triangle])

    plt.tight_layout()
    plt.savefig(savename, dpi=500)

    plt.show()

def plot_cut_offs(plots, ts, ani):
    
    '''
    
    Plots the spin diffusion along their cut-off distances. This is used
    to check if the cut of distances are correct, as they are used for 
    the linear regression to find spin diffusion lengths.
    
    '''

    plt.figure(figsize=(10,7))


    N = len(ts)
    colors = plt.cm.viridis(np.linspace(0,1,N+1))
    # colors = ['tab:blue', 'tab:red']


    for k, plot in enumerate(plots):

        try:
            f = open(plot, 'r')
        except:
            print('No file for V = ', -ts[k])
        else:
            plat = plot.split('/')[-1]
            temps = plat.split('_')
            temp = temps[2]

            lines = f.readlines()
            lines = lines[10:]

            vals = []

            for i in range(len(lines)):
                vec1 = lines[i].split('\t')
                all_vals = vec1[1:]
                ani_int = 0
                if ani == 'OOP':
                    ani_int = 2
                temp = []
                while ani_int < len(all_vals):
                    temp.append(float(all_vals[ani_int]))
                    ani_int += 3
                vals.append(temp)

            ys = []

            for i in range(len(vals[0])):
                val = 0
                for j in range(len(vals)):
                    val += float(vals[j][i])
                val /= len(vals)
                ys.append(val)

            original_ys = copy.deepcopy(ys)
                
            try: 
                min_y = min(ys)
            except ValueError:
                print('min(ys) is empty for V = ', -ts[k])
            else:
                
                # This is the data set we use for curve fitting
                fit_ys = [y + abs(min_y) for y in ys]
                fit_ys = np.array([np.log(p) for p in fit_ys])
                
                # This is the data set we use for determining the cut-off
                cutoff_ys = np.array([np.log(p) for p in ys])
                cutoff_ys  = cutoff_ys[np.isfinite(cutoff_ys)]
                
                xs = np.linspace(0, 2.98, len(fit_ys))
                    
                # Find the cut-off before fitting the function
                algo = rpt.Dynp(model="l2").fit(cutoff_ys)
                try:
                    result = algo.predict(n_bkps=1)
                except:
                    print('Data was too noisy for V = ', -ts[k])
                else:
                    
                    plt.plot(xs, original_ys, color=colors[k])

                    plt.axvline(x = xs[result[0]], color = 'b')

                    plt.xlabel(r'x')
                    plt.ylabel(r'$\mu$')

                    plt.tight_layout()

                    plt.show()             
        
def plot_logarithms(plots, ts, ani):
    
    '''
    
    Plots the logarithms of given magnon transport simulations. 
    Used to check whether data can be used in plot_diffusion_lengths()
    
    '''

    plt.figure(figsize=(10,7))


    N = len(ts)
    colors = plt.cm.viridis(np.linspace(0,1,N+1))
    # colors = ['tab:blue', 'tab:red']


    for k, plot in enumerate(plots):

        try:
            f = open(plot, 'r')
        except:
            print('No file for V = ', -ts[k])
        else:
            plat = plot.split('/')[-1]
            temps = plat.split('_')
            temp = temps[2]

            lines = f.readlines()
            lines = lines[10:]

            vals = []

            for i in range(len(lines)):
                vec1 = lines[i].split('\t')
                all_vals = vec1[1:]
                ani_int = 0
                if ani == 'OOP':
                    ani_int = 2
                temp = []
                while ani_int < len(all_vals):
                    temp.append(float(all_vals[ani_int]))
                    ani_int += 3
                vals.append(temp)

            ys = []

            for i in range(len(vals[0])):
                val = 0
                for j in range(len(vals)):
                    val += float(vals[j][i])
                val /= len(vals)
                ys.append(val)

            min_y = min(ys)
            ys1 = [y + abs(min_y) for y in ys]
            ys1 = np.array([np.log(p) for p in ys1])
            ys = np.array([np.log(p) for p in ys])
            ys  = ys[np.isfinite(ys)]
            xs = np.linspace(0, 2.98, len(ys1))
            
            # Apply to non-transformed data
            algo = rpt.Dynp(model="l2").fit(ys)
            result = algo.predict(n_bkps=1)

            plt.plot(xs, ys1, color=colors[k])
            
            plt.axvline(x = xs[result[0]], color = 'b')
            

            plt.xlabel(r'x')
            plt.ylabel(r'$\mu$')

            plt.tight_layout()
            # print(result[0])
            plt.show()

def plot_data_with_fitted_functions(plots, ani, thickness):
    
    '''
    
    Plots the actual data along with the fitted function for many plots. Used to validate of the 
    calculations are correct.
    
    '''

    plt.figure(figsize=(10,7))


    # N = len(ts)
    # colors = plt.cm.viridis(np.linspace(0,1,N+1))

    for k, plot in enumerate(plots):

        try:
            f = open(plot, 'r')
        except:
            print('No file')
        else:
            plat = plot.split('/')[-1]
            temps = plat.split('_')
            temp = temps[2]

            lines = f.readlines()
            lines = lines[10:]

            vals = []

            for i in range(len(lines)):
                vec1 = lines[i].split('\t')
                all_vals = vec1[1:]
                ani_int = 0
                if ani == 'OOP':
                    ani_int = 2
                temp = []
                while ani_int < len(all_vals):
                    temp.append(float(all_vals[ani_int]))
                    ani_int += 3
                vals.append(temp)

            ys = []

            for i in range(len(vals[0])):
                val = 0
                for j in range(len(vals)):
                    val += float(vals[j][i])
                val /= len(vals)
                ys.append(val)

            original_ys = copy.deepcopy(ys)
            original_xs = np.linspace(0, 2.98, len(original_ys))
            try: 
                min_y = min(ys)
            except ValueError:
                print('min(ys) is empty')
            else:
                
                # This is the logarithm data set we use for curve fitting
                fit_ys = [y + abs(min_y) for y in ys]
                fit_ys = np.array([np.log(p) for p in fit_ys])
                
                # This is the data set we use for determining the cut-off
                cutoff_ys = np.array([np.log(p) for p in ys])
                cutoff_ys  = cutoff_ys[np.isfinite(cutoff_ys)]
                
                plot_ys = copy.deepcopy(fit_ys)
                xs = np.linspace(0, 2.98, len(fit_ys))
                    
                # Find the cut-off before fitting the function
                algo = rpt.Dynp(model="l2").fit(cutoff_ys)
                try:
                    result = algo.predict(n_bkps=1)
                except:
                    print('Data was too noisy')
                else:

                    # Then cut off the list ands create new x-values list
                    if k == 7:
                        result[0] += 20
                    fit_ys = fit_ys[:result[0]]
                    fit_ys2 = original_ys[:result[0]]
                    fit_xs = xs[:result[0]]

                    # Linear fit
                    def f(x, a, b):
                        return a - b*x   
                    
                    # Exponential with harmonic wave
                    def g(x, a, b, c, d, e, f):
                        return a*np.exp(-b*x) + c*np.sin(d*x + e)*np.exp(-f*x)
                    
                    # Trying new function
                    def i(x, a, b, c, d, e):
                        return np.exp(-b*x)*(a + c*np.sin(d*x + e))
                    
                    # Exponential
                    def h(x, a, b):
                        return a*np.exp(-b*x)

                    if thickness == 40 or (k > 3 and k < 8) or (k > 11):
                        # if k == 7:
                        try:
                            params, params_cov = curve_fit(i, fit_xs, fit_ys2, p0=[5e11, 2, 5e11, 20, 0])
                        except TypeError:
                            print('Could not find curve for this voltage')
                        # If we get runtime error then there is no harmonic wave, just exponential decay
                        except RuntimeError:
                            params, params_cov = curve_fit(h, fit_xs, fit_ys2, p0=[1e12, 2])
                            plt.plot(original_xs+3.02, np.array(original_ys)/1e11, color='tab:blue')
                            plt.plot(fit_xs+3.02, np.array(h(fit_xs, *params))/1e11, color='tab:orange')
                        else:
                            plt.plot(original_xs+3.02, np.array(original_ys)/1e11, color="#e76f8a", linewidth=2, marker='s', markevery=[0, 2, 5, 8, 13, 20, 27, 37, 42, 48,58])
                            plt.plot(fit_xs+3.02, np.array(i(fit_xs, *params))/1e11, color="#4575b4", linewidth=2)
                        plt.xlabel(r'x (μm)')
                        plt.ylabel(r'$\mu$')
                        plt.legend(['Measured data', 'Fitted function'])
                        plt.tight_layout()
                        # plt.savefig('AFM/custom/plots/function_fit_uniaxial_8layer_constJc_over_d.png', dpi=500)
                        plt.show()
                            
                    # if k == 2:
                    #     params, params_cov = curve_fit(f, fit_xs, fit_ys)
                    #     plt.plot(xs+3.02, np.array(plot_ys), color="#e76f8a", linewidth=2, marker='s', markevery=range(0, 310,30))
                    #     plt.plot(fit_xs+3.02, np.array(f(fit_xs, params[0], params[1])), color="#4575b4", linewidth=2)
                    #     plt.xlabel(r'x (μm)')
                    #     plt.ylabel(r'ln($\mu$)')
                    #     plt.legend(['Measured data', 'Fitted function'])
                    #     plt.tight_layout()
                    #     plt.savefig('AFM/custom/plots/function_fit_uniaxial_2layer_constJc_over_d.png', dpi=500)
                    #     plt.show()
                                   
def plot_current_density(meshdims, cellsize, t, V, damping, MEC, ani, T, type):

    plt.figure(figsize=(10, 7))

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/plots/' + 'current_density/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


    filename = type + '/' + modules_folder + ani + '/cache/current_density/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/Jc'  + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K.txt'

    f = open(filename, 'r')

    lines = f.readlines()
    lines = lines[10:]

    array = np.array([np.array([line.split('\t')]) for line in lines])

    transposed = array.transpose()

    for elem in transposed:
        ys = []
        for val in elem:
            ys.append(float(val))
        xs = np.linspace(0, t, len(ys))
        plt.plot(xs, ys)

    savename = type + '/' + modules_folder + ani + '/plots/current_density/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/Jc'  + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K.png'
    plt.savefig(savename, dpi=500)
    plt.show()

def fft_transport(meshdims, cellsize, t, V, damping, MEC, ani, T, type):

    plt.figure(figsize=(10, 7))

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/plots/' + 'frequency/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'
    filename_2D = type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'

    if os.path.isfile(filename):
        f = open(filename, 'r')

        lines = f.readlines()
        lines = lines[10:]

        # xs = np.linspace(x_start/1000, x_stop/1000, int((x_stop - x_start)/cellsize))

        vals = []


        for line in lines:
            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if ani == 'OOP':
                direction = 2
            vec = line.split('\t')
            all_vals = vec[1:]
            temp = []
            while direction < len(all_vals):
                temp.append(float(all_vals[direction]))
                direction += 3
            vals.append(temp)


        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)

    elif os.path.isfile(filename_2D):
        f = open(filename_2D, 'r')

        lines = f.readlines()
        lines = lines[10:]

        # Turn the raw data into a list of numpy arrays. Every entry in the arrays are floats.
        raw_data = []
        for line in lines:

            # Make a list of all entries and an empty array to fill only the component we want in
            vec = line.strip().split('\t')[1:]
            temp = np.empty(int(len(vec)/3))

            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if ani == 'OOP':
                direction = 2

            # Iterate over all and add only the component we want. Convert to float
            indexer = 0
            while direction < len(vec):
                temp[indexer] = float(vec[direction])
                indexer += 1
                direction += 3

            # Reshape to 2D array and add to the data list
            raw_data.append(temp.reshape(meshdims[2], int(len(temp)/(meshdims[2]))))

        # Now find the time averages for all the data
        tAvg_data = np.zeros_like(raw_data[0])

        for k, matrix in enumerate(raw_data):
            for i, row in enumerate(matrix):
                for j, col in enumerate(row):
                    tAvg_data[i][j] += col
                    if k == len(raw_data)-1:
                        tAvg_data[i][j] /= len(raw_data)


        
        ys = [tAvg_data[0][i] for i in range(len(tAvg_data[0]))]
        # xs = np.linspace(x_start, x_stop/1000, len(ys))

    else:
        print("No simulations have been done with these params")
        exit()


    # Take the fourier transform of the data
    fourier_data = np.abs(fft(ys))
    freq_len = len(fourier_data)
    freq = fftfreq(freq_len, cellsize*1e-9)

    result = [fourier_data[i] for i in range(0,int(0.5 *freq_len))]

    peaks = find_peaks(result, 0.3e12)
    print(peaks)
    peaks_str = ([freq[peak] for peak in peaks[0]])

    plt.plot(freq[:(freq_len)//2], result)
    plt.annotate(peaks_str, xy=(0.5, 0.8), xycoords='axes fraction')
    savename = type + '/' + modules_folder + ani + '/plots/frequency/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V_' + str(V) + '.png'
    plt.savefig(savename, dpi=500)
    plt.show()

def fft_transport_underneath(meshdims, cellsize, t, V, damping, MEC, ani, T, type):

    plt.figure(figsize=(10, 7))

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/underneath/plots/' + 'frequency/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = type + '/' + modules_folder + ani + '/underneath/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'
    filename_2D = type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'

    if os.path.isfile(filename):
        f = open(filename, 'r')

        lines = f.readlines()
        lines = lines[10:]

        # xs = np.linspace(x_start/1000, x_stop/1000, int((x_stop - x_start)/cellsize))

        vals = []


        for line in lines:
            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if ani == 'OOP':
                direction = 2
            vec = line.split('\t')
            all_vals = vec[1:]
            temp = []
            while direction < len(all_vals):
                temp.append(float(all_vals[direction]))
                direction += 3
            vals.append(temp)


        ys = []

        for i in range(len(vals[0])):
            val = 0
            for j in range(len(vals)):
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)

    elif os.path.isfile(filename_2D):
        f = open(filename_2D, 'r')

        lines = f.readlines()
        lines = lines[10:]

        # Turn the raw data into a list of numpy arrays. Every entry in the arrays are floats.
        raw_data = []
        for line in lines:

            # Make a list of all entries and an empty array to fill only the component we want in
            vec = line.strip().split('\t')[1:]
            temp = np.empty(int(len(vec)/3))

            # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
            direction = 0
            if ani == 'OOP':
                direction = 2

            # Iterate over all and add only the component we want. Convert to float
            indexer = 0
            while direction < len(vec):
                temp[indexer] = float(vec[direction])
                indexer += 1
                direction += 3

            # Reshape to 2D array and add to the data list
            raw_data.append(temp.reshape(meshdims[2], int(len(temp)/(meshdims[2]))))

        # Now find the time averages for all the data
        tAvg_data = np.zeros_like(raw_data[0])

        for k, matrix in enumerate(raw_data):
            for i, row in enumerate(matrix):
                for j, col in enumerate(row):
                    tAvg_data[i][j] += col
                    if k == len(raw_data)-1:
                        tAvg_data[i][j] /= len(raw_data)


        
        ys = [tAvg_data[0][i] for i in range(len(tAvg_data[0]))]
        # xs = np.linspace(x_start, x_stop/1000, len(ys))

    else:
        print("No simulations have been done with these params")
        exit()


    # Take the fourier transform of the data
    fourier_data = np.abs(fft(ys))
    freq_len = len(fourier_data)
    freq = fftfreq(freq_len, cellsize*1e-9)

    result = [fourier_data[i] for i in range(0,int(0.5 *freq_len))]

    peaks = find_peaks(result, 0.3e12)
    print(peaks)
    peaks_str = ([freq[peak] for peak in peaks[0]])

    plt.plot(freq[:(freq_len)//2], result)
    plt.annotate(peaks_str, xy=(0.5, 0.8), xycoords='axes fraction')
    savename = type + '/' + modules_folder + ani + '/underneath/plots/frequency/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/V_' + str(V) + '.png'
    plt.savefig(savename, dpi=500)
    plt.show()

def plot_uniaxial_analytical_dispersion():
    
    '''
    
    Plots the analytical disperison relation for a system with uniaxial anisotropy.
    Plots 2D and 3D in the same plot
    
    '''
    
    # Define the constants
    Ah = -460e3 # J/m^3
    K_easy = 21 # J/m^3
    A = 76e-15 # J/m
    a = 5e-9 # m
    hbar = 1.0546e-34 # m^2 kg/s
    gamma_e_over_2pi = 2.802e10 # s^-1 T^-1
    Ms = 2.1e3 # A/m
    C = gamma_e_over_2pi/(Ms)
    
    # Define the dispersion relation
    def f(q_a):
        q = q_a / 5e-9
        return np.sqrt(C**2 * 16 * np.abs(Ah)*K_easy + 8*A*np.abs(Ah)*q**2*C**2)/1e12
    
    x = np.linspace(-2, 2, 1000)
    
    plt.plot(x, f(x))
    plt.plot(x, f(x))
    
    plt.tight_layout()
    plt.show()
    
def plot_uniaxial_analytical_dispersion_width_modes(magnonDispersion):
    
    '''
    
    Plots the analytical disperison relation for a system with uniaxial anisotropy.
    Plots 2D and 3D in the same plot
    
    '''
    
    
    if magnonDispersion.type == 'AFM':
        time_step = 0.1e-12
        ylabel = 'f (THz)'
        divisor = 1e12
    elif magnonDispersion.type == 'FM':
        time_step = 1e-12
        ylabel = 'f (GHz)'
        divisor = 1e9
    else:
        raise Exception('Choose type!')
    
    
    output_files = magnonDispersion.cachename()
        

    if magnonDispersion.hard_axis:
        fig1,ax1 = plt.subplots(1,2)
        fig1.set_figheight(6)
        fig1.set_figwidth(14)
    else:
        fig1,ax1 = plt.subplots()

    titles = ['y-component', 'z-component']
    
    axes_list = []

    for i, output_file in enumerate(output_files):

        pos_time = np.loadtxt(output_file)

        fourier_data = np.fft.fftshift(np.abs(np.fft.fft2(pos_time)))

        freq_len = len(fourier_data)
        k_len = len(fourier_data[0])
        freq = np.fft.fftfreq(freq_len, time_step)
        kvector = np.fft.fftfreq(k_len, 5e-9)

        k_max = 2*np.pi*kvector[int(0.5 * len(kvector))]*5e-9
        f_min = np.abs(freq[0])
        f_max = np.abs(freq[int(0.5 * len(freq))])/divisor # to make it THz
        f_points = int(0.5 * freq_len)

        final_result = [fourier_data[i] for i in range(int(0.5 *freq_len),freq_len)]

        label = r'$q_{' + magnonDispersion.axis + r'}$a' 

        if magnonDispersion.hard_axis:
            ax1[i].imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0,1000))
            ax1[i].set_xlabel(label)
            ax1[i].set_ylabel(ylabel)


        else:
            ax1.imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0,1000))
            ax1.set_xlabel(label)
            ax1.set_ylabel(ylabel)
        # ax1.set_ylim(0, 0.1)
        
    # Define the constants
    Ah = -460e3 # J/m^3
    K_easy = 21 # J/m^3
    A = 76e-15 # J/m
    a = 5e-9 # m
    hbar = 1.0546e-34 # m^2 kg/s
    gamma_e_over_2pi = 2.802e10 # s^-1 T^-1
    Ms = 2.1e3 # A/m
    C = gamma_e_over_2pi/(Ms)
    
    # Define the dispersion relation
    def g(q_a, t, n, m):
        q = q_a / 5e-9
        return np.sqrt(C**2 * 16 * np.abs(Ah)*K_easy + 16*A*np.abs(Ah)*(q**2 + (n*np.pi/(50e-9)**2) + (m*np.pi/(t*1e-9)**2))*C**2)/1e12
    
    # No width dependence
    def f(q_a, n):
        q = q_a / 5e-9
        return np.sqrt(C**2 * 16 * np.abs(Ah)*K_easy + 16*A*np.abs(Ah)*(q**2 + (n*np.pi/(50e-9)**2))*C**2)/1e12
    
    t = 5
    
    x = np.linspace(-0.75, 0.75, 1000)
    for n in [0, 3, 12, 27, 45]:
        for m in range(1):
            plt.plot(x, g(x, t, n, m), linestyle='dashed', color='red')
    
    plt.tight_layout()
    plt.savefig('AFM/custom/plots/analytical_width_modes.png', dpi=500)
    plt.show()

def main():


    # plot_uniaxial_analytical_dispersion()
    # plot_uniaxial_analytical_dispersion_width_modes()


    f1 = 'AFM/ex+ani/IP/cache/t_avg/6000x50x20/tAvg_damping0.0004_V-0.625_0.3K_sim6.txt'
    f2 = 'AFM/ex+ani/IP/cache/t_avg/6000x50x20/tAvg_damping0.0004_V-0.625_0.3K_sim9.txt'

    savename = 'AFM/custom/plots/4layer_dispersion_comparison_V-0.625_sim6_vs_sim9.png'
    
    # plot_dispersions([f1, f2], savename)

    #### MAGNON TRANSPORT ### 

    def thickness_comparison(hard_axis, V0):


        Vs = []
        ts = []
        for i in range(1, 9):
            # Have to make an exception here because the biaxial system at thickness=10nm and V=-0.12 flips and
            # can not be included
            Vs.append(V0*i)
            ts.append(i)
            
        hard_axis_str = ''
        system = 'uniaxial'
        if hard_axis:
            hard_axis_str = '+hard_axis+Hfield'
            system = 'biaxial'
            
        fs = []
        for V, t in zip(Vs, ts):
            f = 'AFM/ex+ani' + hard_axis_str + f'/IP/cache/t_avg/6000x50x{t*5}/tAvg_damping0.0004_V{V:.3f}_0.3K.txt'
            fs.append(f)

        savename = 'AFM/custom/plots/' + system + f'_mxdmdt_over_thicknesses_V0={V0}.png'

        plot_tAvg_comparison(fs, ts, savename, 'IP')

    thickness_comparison(1, -0.06)
    
    ### VERY thick meshes
    
    def thickness_comparison_thick(hard_axis):

        Vs = [-1.45, -2.175, -2.9]
        ts = [50, 75, 100]
            
        hard_axis_str = ''
        system = 'uniaxial'
        if hard_axis:
            hard_axis_str = '+hard_axis+Hfield'
            system = 'biaxial'
            
        fs = []
        for V, t in zip(Vs, ts):
            f = 'AFM/ex+ani' + hard_axis_str + f'/IP/cache/t_avg/4000x50x{t}/tAvg_damping0.0004_V{V:.3f}_0.3K.txt'
            fs.append(f)

        savename = 'AFM/custom/plots/' + system + f'_mxdmdt_over_thicknesses_V0={Vs[0]}_thick_meshes.png'

        plot_tAvg_comparison_thick(fs, ts, savename, 'IP')

    # thickness_comparison_thick(1)

    ### DISPERSIONS ###

    def magnetic_field_comparison():
        f1 = 'AFM/ex+ani+hard_axis/IP/cache/dispersions/4000x50x5/ycomponent_axisxgroundstate_damping0.0004_T0.3_dispersion.txt'
        f2 = 'AFM/ex+ani+hard_axis/IP/cache/dispersions/4000x50x5/zcomponent_axisxgroundstate_damping0.0004_T0.3_dispersion.txt'
        f3 = 'AFM/ex+ani+hard_axis+Hfield/IP/cache/dispersions/4000x50x5/ycomponent_axisxgroundstate_damping0.0004_T0.3_H2355500.0_dispersion.txt'
        f4 = 'AFM/ex+ani+hard_axis+Hfield/IP/cache/dispersions/4000x50x5/zcomponent_axisxgroundstate_damping0.0004_T0.3_H2355500.0_dispersion.txt'
        f5 = 'AFM/ex+ani+hard_axis+Hfield/IP/cache/dispersions/4000x50x5/ycomponent_axisxgroundstate_damping0.0004_T0.3_H4711000.0_dispersion.txt'
        f6 = 'AFM/ex+ani+hard_axis+Hfield/IP/cache/dispersions/4000x50x5/zcomponent_axisxgroundstate_damping0.0004_T0.3_H4711000.0_dispersion.txt'
        f7 = 'AFM/ex+ani+hard_axis+Hfield/IP/cache/dispersions/4000x50x5/ycomponent_axisxgroundstate_damping0.0004_T0.3_H7066500.0_dispersion.txt'
        f8 = 'AFM/ex+ani+hard_axis+Hfield/IP/cache/dispersions/4000x50x5/zcomponent_axisxgroundstate_damping0.0004_T0.3_H7066500.0_dispersion.txt'

        plots = [f1,f2,f3,f4,f5,f6,f7,f8]
        
        c_max1 = [3000, 3500, 3500, 3500]
        c_max2 = [9000, 10000, 8000, 9000]

        plot_dispersion_field_comparisons(plots, c_max1, c_max2)
        
    # magnetic_field_comparison()

    # f1 = ''

    ### DIFFUSION LENGTHS ###
    
    ### Constant Jc over d ###
    
    def const_Jc_over_d():
        Vs = []
        thicknesses = []
        V0 = 0.06
        fs = []
        for i in range(1, 9):
            Vs.append(V0*i)
            thicknesses.append(5*i)
        
        thicknesses += thicknesses
            
        for k in range(2):
            if k == 0:
                easy_plane_string = ''
            else:
                easy_plane_string = '+hard_axis+Hfield'

            for V, thickness in zip(Vs, thicknesses):
                plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/6000x50x{thickness}/tAvg_damping0.0004_V-{round(V,3):.3f}_0.3K.txt'
                fs.append(plotname)
                
        
        savename_thickness = f'AFM/custom/plots/diffusion_length_across_thicknesses_V0={V0:.3f}.png'
        
        plot_diffusion_length_across_thicknesses(fs, thicknesses, savename_thickness, 'IP')
        # plot_cut_offs(fs, thicknesses, 'IP')
        # plot_data_with_fitted_functions(fs, 'IP', 0)
    
    const_Jc_over_d()
    
    ### Across voltages ###
    
    def across_voltages(thickness, easy_plane):

        fs = []
        ts = []
        factor = thickness/20
        easy_plane_string = ''
        if easy_plane:
            easy_plane_string = '+hard_axis+Hfield'
        for i in range(3,55):
            voltage = (i+1)*0.025
            ts.append(round(voltage*factor,3))
            plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/4000x50x{thickness}/tAvg_damping0.0004_V-{round(voltage*factor,3):.3f}_0.3K.txt'
            fs.append(plotname)

        savename = 'AFM/ex+ani' + easy_plane_string + f'/IP/plots/custom/6micro_voltage_dependence_diffusion_lengths_{int(thickness/5)}layer.png'

        # plot_diffusion_length(fs, ts, savename, 'IP', thickness)
        # plot_logarithms(fs, ts, 'IP')
        # plot_cut_offs(fs, ts, 'IP')
        plot_data_with_fitted_functions(fs, 'IP', thickness)
    
    # across_voltages(40, 0)
    
    ## Both systems across voltages ###
    
    def both_systems_across_voltages(thickness):
        
        start = 9
        length = 4000
        if thickness == 10:
            length = 6000
            start = 4
        elif thickness == 20:
            start = 7
            length = 6000
            
        
        easy_plane = 0
        
        fs1 = []
        ts = []
        factor = thickness/20
        easy_plane_string = ''
        if easy_plane:
            easy_plane_string = '+hard_axis+Hfield'
            
        # Vs1 = np.linspace(0.200, 0.575, 6)
        # Vs2 = np.linspace(0.175,0.550, 6)
        # Vs = np.concatenate((Vs1, Vs2))
        Vs = [0.175, 0.200, 0.250, 0.275, 0.325, 0.350, 0.400, 0.425, 0.475, 0.500, 0.550, 0.575]
        Vs = [0.175, 0.200, 0.250, 0.275, 0.325, 0.350, 0.400, 0.425, 0.500, 0.550, 0.575]
            
        for i in Vs:
            voltage = (i)
            # ts.append(round(voltage*factor,3))
            ts.append(factor*voltage)
            plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/{length}x50x{thickness}/tAvg_damping0.0004_V-{round(factor*voltage,3):.3f}_0.3K.txt'
            fs1.append(plotname)
        
        easy_plane = 1
        
        fs2 = []
        factor = thickness/20
        easy_plane_string = ''
        if easy_plane:
            easy_plane_string = '+hard_axis+Hfield'
        for i in Vs:
            voltage = (i)
            plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/{length}x50x{thickness}/tAvg_damping0.0004_V-{round(voltage*factor,3):.3f}_0.3K.txt'
            fs2.append(plotname)
            
        savename_both = f'AFM/custom/plots/voltage_dependence_diffusion_lengths_{int(thickness/5)}layer.png'
        
        plot_diffusion_length_both_systems(fs1, fs2, ts, savename_both, 'IP', thickness)
    
    # both_systems_across_voltages(10)
    
    ### All systems across voltages ###
    
    # all_fs = []
    # all_ts = []
    # for thickness in [10,20,40]:
        
    #     fs1 = []
    #     ts = []
    #     factor = thickness/20
    #     easy_plane_string = ''
    #     for i in range(4,55):
    #         voltage = (i+1)*0.025
    #         ts.append(round(voltage*factor,3))
    #         plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/4000x50x{thickness}/tAvg_damping0.0004_V-{round(voltage*factor,3):.3f}_0.3K.txt'
    #         fs1.append(plotname)
            
    #     all_fs.append(fs1)
    #     all_ts.append(ts)
        
    # for thickness in [10,20,40]:

    #     fs2 = []
    #     factor = thickness/20
    #     easy_plane_string = '+hard_axis+Hfield'
    #     for i in range(4,55):
    #         voltage = (i+1)*0.025
    #         ts.append(round(voltage*factor,3))
    #         plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/4000x50x{thickness}/tAvg_damping0.0004_V-{round(voltage*factor,3):.3f}_0.3K.txt'
    #         fs2.append(plotname)
            
    #     all_fs.append(fs2)
    #     all_ts.append(ts)
        
        
    # savename_all = f'AFM/custom/plots/voltage_dependence_diffusion_lengths_all_layers.png'
    
    # plot_diffusion_length_all_systems(all_fs, all_ts, savename_all, 'IP')

if __name__ == '__main__':
    main()