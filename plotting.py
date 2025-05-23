import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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

def plot_tAvg_SA(timeAvgSA):

    '''
    
    The standard function for plotting <mxdmdt> over distance in magnon transport simulations. 

    '''

    plt.figure(figsize=(10, 7))

    filename = timeAvgSA.cachename()

    try:
        f = open(filename, 'r')

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
                val += float(vals[j][i])
            val /= len(vals)
            ys.append(val)
            
        xs = np.linspace(timeAvgSA.x_start/1000, timeAvgSA.x_stop/1000, len(ys))
            
        plt.plot(xs, ys, linewidth=2)
        plt.xlabel(r'$x$ $(\mu m)$')
        plt.ylabel(r'$\mu_x$')
        plt.tight_layout()

        plotname = timeAvgSA.plotname()
        params.make_folder(plotname)
        plt.savefig(plotname, dpi=500)
        
        plt.close()
        # plt.show()
        
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
    plt.savefig(plotname, dpi=500)
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

    plt.savefig(savename, dpi=500)
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

    plt.savefig(savename, dpi=500)
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
                [0.5, 0.5, 0.47, 0.47],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower')

            ax1[i].indicate_inset_zoom(axins, edgecolor='black')

            ax1[i].title.set_text(titles[i])

        else:
            ax1.imshow(final_result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))
            ax1.set_xlabel(label)
            ax1.set_ylabel(ylabel)

            x1, x2, y1, y2 = -0.3, 0.3, 0, 0.4
            axins = ax1.inset_axes(
                [0.4, 0.52, 0.6, 0.4],
                xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
            axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,5000))
            axins.set_yticks([0.0, 0.2, 0.4])

            ax1.indicate_inset_zoom(axins, edgecolor='white', alpha=1)
        
            for spine in axins.spines.values():
                spine.set_edgecolor('white')
                
            axins.tick_params(colors='white')
        # ax1.set_ylim(0, 0.1)

        plt.tight_layout()

    plt.savefig(savename, dpi=500)

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

    plt.savefig(savename, dpi=500)

    plt.show()

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
    params.make_folder(savename)
    
    fig, ax = plt.subplots(1,3)

    fig.set_figheight(4)
    fig.set_figwidth(12)
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

        x1, x2, y1, y2 = -0.75, 0.75, 0, 1.5
        axins = ax[index].inset_axes(
            [0.5, 0.4, 0.47, 0.47],
            xlim=(x1,x2), ylim=(y1,y2), xticklabels=[])
        axins.imshow(final_result, extent = [-k_max, k_max,f_min, f_max], origin='lower', clim=(0,5000))

        ax[index].indicate_inset_zoom(axins, edgecolor='black')

    plt.tight_layout()

    plt.savefig(savename, dpi=500)

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
        
            label = 'qy'

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
    
    savename = 'AFM/custom/plots/dispersions_magnetic_field_comparison.png'

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
    plt.xlabel(r'Temperature ($K$)')
    plt.ylabel(r'|${m}_{x}$|')
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
    colors = ['tab:blue', 'tab:orange']

    for m, plots in enumerate(both):
        
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


                # When <mxdmdt> flattens out we get NaN and inf after doing the logarithm. 
                # This needs to be removed
                # From looking at the data, it seems to be fine to just find the 
                # slope from the first half of ys, as the end point introduces lots of noise
                try: 
                    min_y = min(ys)
                except ValueError:
                    print('min(ys) is empty for V = ', -ts[k])
                else:
                    
                    # This is the data set we use for curve fitting
                    ys1 = [y + abs(min_y) for y in ys]
                    ys1 = np.array([np.log(p) for p in ys1])
                    
                    # This is the data set we use for determining the cut-off
                    ys = np.array([np.log(p) for p in ys])
                    ys  = ys[np.isfinite(ys)]
                    
                    xs = np.linspace(0, 2.98, len(ys1))
                        
                    # Find the cut-off before fitting the function
                    algo = rpt.Dynp(model="l2").fit(ys)
                    try:
                        result = algo.predict(n_bkps=1)
                    except:
                        print('Data was too noisy for V = ', -ts[k])
                    else:

                        # Then cut off the list ands create new x-values list
                        ys1 = ys1[:result[0]]
                        xs = xs[:result[0]]

                        def f(x, a, b):
                            return a - b*x   
                        
                        def g(x, a, b, c, d):
                            return a*np.exp(-b*x) * np.sin(c*x + d)

                        try:
                            if thickness == 40:
                                params, params_cov = curve_fit(g, xs, ys1)
                                plt.plot(10e3*ts[k]/thickness, 1/params[1] / params[0], marker=markers[m], markersize = 10, color=colors[m])
                            else:
                                params, params_cov = curve_fit(f, xs, ys1)
                                plt.plot(10e3*ts[k]/thickness, 1/params[1], marker=markers[m], markersize = 10, color=colors[m])
                        except TypeError:
                            print('Could not find a curve for this voltage: ', ts[k])
                        except ValueError:
                            print('ydata was empty for this voltage: ', ts[k])
                        
                        
    plt.xlabel(r'Voltage/thickness (μV/nm)')
    plt.ylabel(r'$l_d$')
    
    
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

    for k, plot in enumerate(plots):

        print(plot)

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


            # When <mxdmdt> flattens out we get NaN and inf after doing the logarithm. 
            # This needs to be removed
            # From looking at the data, it seems to be fine to just find the 
            # slope from the first half of ys, as the end point introduces lots of noise
            try: 
                min_y = min(ys)
            except ValueError:
                print('min(ys) is empty for V = ', ts[k])
            else:
                ys = np.array([np.log(p) for p in ys])
                ys = ys[np.isfinite(ys)]
                xs = np.linspace(0, 2, len(ys))
                    
                # Find the cut-off before fitting the function
                algo = rpt.Dynp(model="l2").fit(ys)
                result = algo.predict(n_bkps=1)
                
                # If the signal travels almost to the edges the cut-off is not needed
                # If the cut-off is larger than 1/2 of the distance, this is the case
                print(result[0])
                if result[0] > 90:
                
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
                    plt.plot(ts[k]/5, 1/params[1], marker='v', markersize = 10, color='black')
                except TypeError:
                    print('Could not find a curve for this voltage: ', ts[k])
                except ValueError:
                    print('ydata was empty for this voltage: ', ts[k])

        
    plt.xlabel(r'Thickness (layers)')
    plt.ylabel(r'$l_d$')
    # plt.yticks([0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50])
    plt.xticks([1,2,3,4,5,6,7,8])

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

            
            # If the signal is small enough, it could be negative. It needs to all be positive 
            # for us to find a logarithm
            # When <mxdmdt> flattens out we get NaN and inf after doing the logarithm. 
            # This needs to be removed
            try:
                min_y = min(ys)
            except ValueError:
                print('min(ys) is empty for V = ', ts[k])
            else:
                ys1 = [y + abs(min_y) for y in ys]
                ys1 = np.array([np.log(p) for p in ys])
                # ys1 = ys1[np.isfinite(ys1)]
                
                # Apply to non-transformed data
                algo = rpt.Dynp(model="l2").fit(ys1)
                result = algo.predict(n_bkps=1)
                
                xs = np.linspace(0, 2.98, len(ys))

                plt.plot(xs, ys, color=colors[k])

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

def plot_data_with_fitted_functions(plots, ts, ani, thickness):
    
    '''
    
    Plots the actual data along with the fitted function for many plots. Used to validate of the 
    calculations are correct.
    
    '''

    plt.figure(figsize=(10,7))


    N = len(ts)
    colors = plt.cm.viridis(np.linspace(0,1,N+1))

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
            original_xs = np.linspace(0, 2.98, len(original_ys))
            # If the signal is small enough, it could be negative. It needs to all be positive 
            # for us to find a logarithm
            # When <mxdmdt> flattens out we get NaN and inf after doing the logarithm. 
            # This needs to be removed
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
                    
                # Find the cut-off before fitting the function
                algo = rpt.Dynp(model="l2").fit(ys)
                try:
                    result = algo.predict(n_bkps=1)
                except:
                    print('Data was too noisy for V = ', -ts[k])
                else:

                    # TFor logarithmic function fitting
                    ys1 = ys1[:result[0]]
                    xs1 = xs[:result[0]]
                    
                    # For exponential function fitting
                    ys2 = original_ys[:result[0]]
                    # ys2 = [y - abs(min_y) for y in ys2]
                    xs2 = original_xs[:result[0]]

                    def f(x, a, b):
                        return a - b*x   
                    
                    def g(x, a, b, c, d):
                        return a*np.exp(-b*x) * np.sin(c*x + d)

                    try:
                        if thickness == 40:
                            params, params_cov = curve_fit(g, xs2, ys2, p0=[1e12, 2.5, 10, 0])
                            plt.plot(xs2, ys2, color='tab:blue')
                            def h(x):
                                return params[0]*np.exp(-params[1]*x) * np.sin(params[2]*x+params[3])
                            plt.plot(xs2, h(xs2), color='tab:orange')
                        else:
                            params, params_cov = curve_fit(f, xs1, ys1)
                            plt.plot(xs1, ys1, color='tab:blue')
                            plt.plot(xs1, f(xs1, params[0], params[1]), color='tab:orange')
                        
                        plt.xlabel(r'x')
                        plt.ylabel(r'$\mu$')
                        plt.tight_layout()
                        plt.show()
                        
                    except TypeError:
                        print('Could not find a curve for this voltage: ', ts[k])
                    # except ValueError:
                    #     print('ydata was empty for this voltage: ', ts[k])

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

    # Jc_dict = {5 : '-0.06', 10: '-0.12', 15 : '-0.18', 20 : '-0.25', 25 : '-0.32',
    #              30 : '-0.38', 35 : '-0.45', 40 : '-0.52', 45 : '-0.58'}
    

    # this = 45

    # f1 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x' + str(this) + '/tAvg_damping0.0004_V' + Jc_dict[this] + '_0.3K.txt'
    # f2 = 'AFM/ex+ani+hard_axis+Hfield/IP/cache/t_avg/4000x50x' + str(this) + '/tAvg_damping0.0004_V' + Jc_dict[this] + '_0.3K.txt'
    # f3 = 'AFM/ex+ani/IP/cache/t_avg/3000x50x100/tAvg_damping0.0004_V-4.5_10K.txt'
    # f4 = 'AFM/ex+ani/IP/cache/t_avg/3000x50x100/tAvg_damping0.0004_V-4.5_20K.txt'

    # l1 = 'Easy-axis'
    # l2 = 'Easy-plane'
    # # l3 = '10'
    # l4 = '20'

    # savename = 'AFM/ex+ani+hard_axis+Hfield/IP/plots/custom/easy_vs_hard_comparison_' + str(this/5) + 'layer.png'

    # plot_tAvg_comparison([f1,f2], [l1,l2], savename, 'IP')


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
                plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/4000x50x{thickness}/tAvg_damping0.0004_V-{round(V,3):.3f}_0.3K.txt'
                fs.append(plotname)
                
        
        savename_thickness = f'AFM/custom/plots/diffusion_length_across_thicknesses_V0={V0:.3f}.png'
        
        plot_diffusion_length_across_thicknesses(fs, thicknesses, savename_thickness, 'IP')
        plot_cut_offs(fs, thicknesses, 'IP')
    
    ### Across voltages ###
    
    def across_voltages():
        thickness = 40
        easy_plane = 0

        fs = []
        ts = []
        factor = thickness/20
        easy_plane_string = ''
        if easy_plane:
            easy_plane_string = '+hard_axis+Hfield'
        for i in range(3,14):
            voltage = (i+1)*0.025
            ts.append(round(voltage*factor,3))
            plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/4000x50x{thickness}/tAvg_damping0.0004_V-{round(voltage*factor,3):.3f}_0.3K.txt'
            fs.append(plotname)

        savename = 'AFM/ex+ani' + easy_plane_string + f'/IP/plots/custom/6micro_voltage_dependence_diffusion_lengths_{int(thickness/5)}layer.png'

        # plot_diffusion_length(fs, ts, savename, 'IP', thickness)
        # plot_logarithms(fs, ts, 'IP')
        # plot_cut_offs(fs, ts, 'IP')
        plot_data_with_fitted_functions(fs, ts, 'IP', thickness)
    
    across_voltages()
    
    ## Both systems across voltages ###
    
    def both_systems_across_voltages():
        thickness = 10
        
        
        start = 1
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
        for i in range(start,55):
            voltage = (i+1)*0.025
            ts.append(round(voltage*factor,3))
            plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/{length}x50x{thickness}/tAvg_damping0.0004_V-{round(voltage*factor,3):.3f}_0.3K.txt'
            fs1.append(plotname)
            
        easy_plane = 1
        
        fs2 = []
        factor = thickness/20
        easy_plane_string = ''
        if easy_plane:
            easy_plane_string = '+hard_axis+Hfield'
        for i in range(start,55):
            voltage = (i+1)*0.025
            plotname = 'AFM/ex+ani' + easy_plane_string + f'/IP/cache/t_avg/{length}x50x{thickness}/tAvg_damping0.0004_V-{round(voltage*factor,3):.3f}_0.3K.txt'
            fs2.append(plotname)
            
        savename_both = f'AFM/custom/plots/voltage_dependence_diffusion_lengths_{int(thickness/5)}layer.png'
        
        plot_diffusion_length_both_systems(fs1, fs2, ts, savename_both, 'IP', thickness)
    
    # both_systems_across_voltages()
    
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