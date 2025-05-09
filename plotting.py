import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
import os
import params
from itertools import chain
from scipy.optimize import curve_fit
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
from pathlib import Path

plt.rcParams.update({'font.size': 26})

path = Path(__file__).resolve().parent

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
    fig.savefig(plotname, dpi=600)
    plt.show()

def plot_tAvg_SA(timeAvgSA):

    '''
    
    The standard function for plotting <mxdmdt> over distance in magnon transport simulations. 

    '''

    plt.figure(figsize=(10, 7))

    filename = timeAvgSA.cachename()

    if os.path.isfile(filename):
        f = open(filename, 'r')

        lines = f.readlines()
        lines = lines[10:]

        xs = np.linspace(timeAvgSA.x_start/1000, timeAvgSA.x_stop/1000, int((timeAvgSA.x_stop - timeAvgSA.x_start)/timeAvgSA.cellsize))

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
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)

    # elif os.path.isfile(filename_2D):
    #     f = open(filename_2D, 'r')

    #     lines = f.readlines()
    #     lines = lines[10:]

    #     # Turn the raw data into a list of numpy arrays. Every entry in the arrays are floats.
    #     raw_data = []
    #     for line in lines:

    #         # Make a list of all entries and an empty array to fill only the component we want in
    #         vec = line.strip().split('\t')[1:]
    #         temp = np.empty(int(len(vec)/3))

    #         # This is the component we look at. In plane means x-component (0) and out-of-plane means z (2)
    #         direction = 0
    #         if ani == 'OOP':
    #             direction = 2

    #         # Iterate over all and add only the component we want. Convert to float
    #         indexer = 0
    #         while direction < len(vec):
    #             temp[indexer] = float(vec[direction])
    #             indexer += 1
    #             direction += 3

    #         # Reshape to 2D array and add to the data list
    #         raw_data.append(temp.reshape(meshdims[2], int(len(temp)/(meshdims[2]))))

    #     # Now find the time averages for all the data
    #     tAvg_data = np.zeros_like(raw_data[0])

    #     for k, matrix in enumerate(raw_data):
    #         for i, row in enumerate(matrix):
    #             for j, col in enumerate(row):
    #                 tAvg_data[i][j] += col
    #                 if k == len(raw_data)-1:
    #                     tAvg_data[i][j] /= len(raw_data)


        
    #     ys = [tAvg_data[0][i] for i in range(len(tAvg_data[0]))]
    #     xs = np.linspace(x_start, x_stop/1000, len(ys))

    else:
        print("No simulations have been done with these params")
        exit()

    plt.plot(xs, ys, linewidth=2)
    plt.xlabel(r'$x$ $(\mu m)$')
    plt.ylabel(r'$\mu_x$')
    plt.tight_layout()

    plotname = timeAvgSA.plotname()
    params.make_folder(plotname)
    plt.savefig(plotname, dpi=600)
    plt.show()

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
    plt.savefig(plotname, dpi=600)
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

def plot_tAvg_SA_z(meshdims, cellsize, t, V, damping, MEC, ani, T, type):

    '''
    
    Similar function to the standard, but plots <mxdmdt> along z-direction directly underneath the injector.

    NOTE: This function should be updated to the object-oriented version.
    
    '''

    plt.figure(figsize=(10, 7))

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/plots/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/z_dir_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'
    filename_2D = type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'

    if os.path.isfile(filename):
        f = open(filename, 'r')

        lines = f.readlines()
        lines = lines[10:]

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

        # ys = [y/ys[-1] for y in ys]
        xs = np.linspace(0, meshdims[2]/1000, len(ys))
        ys.reverse()

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
        xs = np.linspace(0, meshdims[2]/1000, len(ys))

    else:
        print("No simulations have been done with these params")
        exit()

    plt.plot(xs, ys, linewidth=2)
    plt.xlabel(r'$z$ $(\mu m)$')
    plt.ylabel(r'$\mu_x$')
    plt.tight_layout()

    plotname = type + '/' + modules_folder + ani + '/plots/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/z_dir_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.png'
    plt.savefig(plotname, dpi=600)
    # plt.show()

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

    plt.savefig(savename, dpi=600)
    plt.show()

def plot_tAvg_comparison(plots, legends, savename, ani):

    '''
    
    Plots several systems in the same plot for comparison of spin diffusion lengths.
    
    '''

    plt.figure(figsize=(13,8))

    N = len(plots)
    # colors = plt.cm.viridis(np.linspace(0,1,N+1))
    # colors = ['teal', 'darkviolet']

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

        ys = [(y )/ys[0] for y in ys]
        xs = np.linspace(2.02, 4, len(ys))

        plt.plot(xs, ys, label=legends[k], linewidth=3)#, color=colors[k])
        # plt.plot(xs, ys, label=str(leg), linewidth=3, color=colors[k])

    plt.xlabel('x (μm)')
    plt.ylabel(r'$\mu$ (normalized)')

    # plt.legend(title= r'$\mu$ value underneath injector ($10^{11}$)')
    plt.legend(title='Thickness (layers)')

    plt.savefig(savename, dpi=600)
    plt.show()


#### MAGNON DISPERSION PLOTS ####

def plot_magnon_dispersion_with_zoom(magnonDispersion, clim_max = 1000):

    '''
    
    Plot the magnon dispersion of a given system, including a zoom-in on an area, usually
    the gap. In the case of easy-plane, y and z components are superimposed on the same plot
    
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

    if magnonDispersion.hard_axis:
        fig1,ax1 = plt.subplots(1,2)
        fig1.set_figheight(7)
        fig1.set_figwidth(16)
    else:
        fig1,ax1 = plt.subplots()

    titles = ['y-component', 'z-component']

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

        label = 'q' + magnonDispersion.axis

        if magnonDispersion.hard_axis:
            ax1[i].imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))
            ax1[i].set_xlabel(label)
            ax1[i].set_ylabel(ylabel)

            x1, x2, y1, y2 = -0.15, 0.15, 0, 0.3
            axins = ax1[i].inset_axes(
                [0.5, 0.4, 0.47, 0.47],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower')

            ax1[i].indicate_inset_zoom(axins, edgecolor='black')

            ax1[i].title.set_text(titles[i])

        else:
            ax1.imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))
            ax1.set_xlabel(label)
            ax1.set_ylabel(ylabel)

            x1, x2, y1, y2 = -0.25, 0.25, 0, 0.5
            axins = ax1.inset_axes(
                [0.5, 0.5, 0.47, 0.47],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,5000))
            axins.set_yticks([0.0, 0.2, 0.4])

            ax1.indicate_inset_zoom(axins, edgecolor='black')
        

        # ax1.set_ylim(0, 0.1)

        plt.tight_layout()

    plt.savefig(savename, dpi=600)

    plt.show()

def plot_magnon_dispersion(magnonDispersion, clim_max = 1000):

    '''
    
    Plot the magnon dispersion of a given system. In the case of easy-plane, the y and z components
    are superimposed on the same plot. 

    '''

    if magnonDispersion.type == 'AFM':
        time_step = 0.1e-12
        ylabel = 'f (THz)'
        divisor = 1e12
    elif magnonDispersion.type == 'FM':
        time_step = 1e-12
        ylabel = 'f (GHz)'
        divisor = 1e9


    if magnonDispersion.steadystate:
        output_files = magnonDispersion.cachename()
        savename = magnonDispersion.plotname()
        
        fig, ax = plt.subplots(1,3)

        fig.set_figheight(5)
        fig.set_figwidth(16)
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

            label = 'q' + magnonDispersion.axis

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
        
        output_files = magnonDispersion.cachename()
        savename = magnonDispersion.plotname()
        params.make_folder(savename)
        
        fig, ax = plt.subplots(1,2)

        fig.set_figheight(5)
        fig.set_figwidth(10)
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

            label = 'q' + magnonDispersion.axis

            ax[i].set_xlabel(label)
            ax[i].set_ylabel(ylabel)
            # ax1.set_ylim(0, 0.1)

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

        label = 'q' + magnonDispersion.axis

        ax1.set_xlabel(label)
        ax1.set_ylabel(ylabel)
        # ax1.set_ylim(0, 0.1)

        plt.tight_layout()

        folder_name = '/'.join(savename.split('/')[:-1])
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

    plt.savefig(savename, dpi=600)

    # plt.show()

def plot_magnon_dispersion_triple_with_zoom(magnonDispersion, clim_max = 1000):

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
    
    fig, ax = plt.subplots(1,3)

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
            pos_time = np.loadtxt(output_files[i+j-1])

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

        label = 'q' + magnonDispersion.axis

        ax[index].set_xlabel(label)
        ax[index].set_ylabel(ylabel)
        title = r'y = ' + str(yvals[index]*1e-3) + ' $\mu$m'
        ax[index].title.set_text(title)

        x1, x2, y1, y2 = -0.25, 0.25, 0, 0.5
        axins = ax[index].inset_axes(
            [0.5, 0.4, 0.47, 0.47],
            xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
        axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,5000))

        ax[index].indicate_inset_zoom(axins, edgecolor='black')


    folder_name = '/'.join(savename.split('/')[:-2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    plt.tight_layout()

    plt.savefig(savename, dpi=600)

    plt.show()

def plot_magnon_dispersion_separate(magnonDispersion, clim_max = 1000):

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

    plt.savefig(savename, dpi=600)

    plt.show()

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

    plt.savefig(savename, dpi=600)

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

    plt.savefig(savename, dpi=600)

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

    annotations = ['a', 'b', 'c', 'd']
    # titles = ['1 layer', '50 layers', '100 layers', '150 layers', '150 layers \n + damping layer']
    titles = ['T = 0.3K', 'T = 0.8K', 'T = 3.0K']

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

        label = 'qy'

        x1, x2, y1, y2 = -0.25, 0.25, 0, 0.5
        axins = ax.inset_axes(
            [0.5, 0.45, 0.47, 0.47],
            xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])

        ax.imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max[i]))
        axins.imshow(result, extent = [-k_max, k_max,f_min, f_max], origin='lower')
        ax.indicate_inset_zoom(axins, edgecolor='black')

        # ax.annotate(annotations[i], (0.05, 0.85), xycoords = 'axes fraction', color='white', fontsize=32)

        # if i == 0:
        #     ax.title.set_text('Without MEC')
        # elif i == 1:
        #     ax.title.set_text('With MEC')

        ax.title.set_text(titles[i])
        # ax.set_xticks([-2,2])
        ax.set_xlabel(label)
        if i == 0:
            ax.set_ylabel('f (THz)')
    # ax1.set_ylim(0, 0.1)

    fig.tight_layout()
    # fig.subplots_adjust(wspace=0.03)
    

    plt.savefig(savename, dpi=600)


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
    plt.savefig(savename, dpi=600)
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
    plt.xlabel(r'Temperature ($K$)')
    plt.ylabel(r'|${m}_{x}$|')
    plt.tight_layout()
    plt.savefig(savename, dpi=600)
    plt.show()

def plot_diffusion_length(plots, ts, savename, ani, thickness):
    
    '''
    
    Plots the spin diffusion lengths of given magnon transport simulations. ts are the values that is compared between them
    
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


            # ys = [(((y )/abs(ys[0]))) for y in ys]
            # min_y = min(ys)
            # ys =[y + abs(min_y) for y in ys]
            

            # When <mxdmdt> flattens out we get NaN and inf after doing the logarithm. 
            # This needs to be removed
            # From looking at the data, it seems to be fine to just find the 
            # slope from the first half of ys, as the end point introduces lots of noise
            ys = np.array([np.log(p) for p in ys])
            ys = ys[np.isfinite(ys)]
            ys = ys[:int(len(ys)/2)]
            xs = np.linspace(0, 2, len(ys))

            def f(x, a, b):
                return a - b*x   

            try:
                params, params_cov = curve_fit(f, xs, ys)
                plt.plot(10e3*ts[k]/thickness, 1/params[1], marker='v', markersize = 10, color='black')
            except TypeError:
                print('Could not find a curve for this voltage: ', ts[k])
            except ValueError:
                print('ydata was empty for this voltage: ', ts[k])

        # plt.plot(ts[k], params[1], marker='v', markersize = 10, color='black')
        # plt.plot(f(xs, params[0], params[1]), xs, color=colors[k])
        # plt.plot(xs, ys, color=colors[k])

        

    plt.xlabel(r'Voltage/thickness (V/nm)')
    plt.ylabel(r'$L_d$')

    plt.tight_layout()
    plt.savefig(savename, dpi=600)

    plt.show()

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

def main():


    #### MAGNON TRANSPORT ### 

    Jc_dict = {5 : '-0.06', 10: '-0.12', 15 : '-0.18', 20 : '-0.25', 25 : '-0.32',
                 30 : '-0.38', 35 : '-0.45', 40 : '-0.52', 45 : '-0.58'}
    

    this = 45

    f1 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x' + str(this) + '/tAvg_damping0.0004_V' + Jc_dict[this] + '_0.3K.txt'
    f2 = 'AFM/ex+ani+hard_axis+Hfield/IP/cache/t_avg/4000x50x' + str(this) + '/tAvg_damping0.0004_V' + Jc_dict[this] + '_0.3K.txt'
    # f3 = 'AFM/ex+ani/IP/cache/t_avg/3000x50x100/tAvg_damping0.0004_V-4.5_10K.txt'
    # f4 = 'AFM/ex+ani/IP/cache/t_avg/3000x50x100/tAvg_damping0.0004_V-4.5_20K.txt'

    l1 = 'Easy-axis'
    l2 = 'Easy-plane'
    # # l3 = '10'
    # l4 = '20'

    savename = 'AFM/ex+ani+hard_axis+Hfield/IP/plots/custom/easy_vs_hard_comparison_' + str(this/5) + 'layer.png'

    # plot_tAvg_comparison([f1,f2], [l1,l2], savename, 'IP')


    ### DISPERSIONS ###

    # f1 = 'AFM/ex+ani+hard_axis/IP/cache/dispersions/4000x50x5/diry_axisxgroundstate_damping0.0004_T0.3_dispersion.txt'
    # f2 = 'AFM/ex+ani+hard_axis/IP/cache/dispersions/4000x50x5/diry_axisxgroundstate_damping0.0004_T0.8_dispersion.txt'
    # f3 = 'AFM/ex+ani+hard_axis/IP/cache/dispersions/4000x50x5/diry_axisxgroundstate_damping0.0004_T3_dispersion.txt'
    # # f4 = 'IP/cache/dispersions/1000x50x150/steady/diry_axisx_dispersion.txt'
    # # f5 = 'IP/cache/dispersions/1000x50x190/steady/diry_axisx_dispersion.txt'


    # savename = 'AFM/ex+ani+hard_axis/IP/plots/custom/dispersion_temperature_comparison.png'

    # plot_dispersions((f1,f2,f3), savename)

    ### DIFFUSION LENGTHS ###

    thickness = 5
    easy_plane = 0

    fs = []
    ts = []
    factor = thickness/20
    easy_plane_string = ''
    if easy_plane:
        easy_plane_string = '+hard_axis+Hfield'
    for i in range(27):
        voltage = (i+1)*0.025
        ts.append(voltage*factor)
        plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/4000x50x{thickness}/tAvg_damping0.0004_V-{voltage*factor:.3f}_0.3K.txt'
        fs.append(plotname)

    savename = 'AFM/ex+ani' + easy_plane_string + f'/IP/plots/custom/voltage_dependence_diffusion_lengths_{int(thickness/5)}layer.png'

    plot_diffusion_length(fs, ts, savename, 'IP', thickness)


if __name__ == '__main__':
    main()