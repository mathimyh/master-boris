import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
import os
from itertools import chain
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 26})

path = 'C:/Users/mathimyh/master/master-boris/'

def plot_something(filename):

    f = open(filename, 'r')

    lines = f.readlines()

    x = []

    for i in range(len(lines)):
        x.append(i)

    plt.plot(x, lines)
    plt.show()

def SA_plotting(filename, plotname, title):

    f = open(filename, 'r')

    lines = f.readlines()

    data = lines[70].strip()

    data = data.split('\t')
    data = data[1:]


    # Calculate only one of the three dimensions
    #
    ys = []
    i = 0
    while i < len(data):
        ys.append(float(data[i]))
        i += 3

    xs = []

    for i in range(len(ys)):
        xs.append(i)

    plt.plot(xs, ys)
    plt.show()
    # plt.xlabel("Distance from injector (nm)")
    # plt.ylabel("Spin accumulation in x-direction")
    # plt.title(title)

    # plotname = "plots/" + plotname
    # plt.savefig(plotname, dpi=500)

    # Total length of the vector

    # ys = []

    # for i in range(0, len(data), 3):
    #     ys.append(np.sqrt(float(data[i])**2 + float(data[i+1])**2 + float(data[i+2])**2))

    # xs = []

    # for i in range(len(ys)):
    #     xs.append(i)

    # plt.plot(xs, ys)
    # plt.xlabel("Distance from injector (nm)")
    # plt.ylabel("Spin accumulation")
    # plt.title(title)

    # plotname = "plots/" + plotname
    # plt.savefig(plotname, dpi=500)

def plot_plateau(meshdims, V, damping, x_vals, MEC, ani, T, type):

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/plots/plateau/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    print(folder_name)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


    filename = type + '/' + modules_folder + ani + '/cache/plateau/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/plateau_V'  + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K.txt'
    f = open(filename, 'r')

    lines = f.readlines()
    lines = lines[11:]

    indexer = 0

    fig, ax = plt.subplots(2, 3)

    fig.set_figheight(10)
    fig.set_figwidth(14)

    direction = 1
    if ani == 'OOP':
        direction = 3

    for i in range(direction, len(lines[0].split('\t'))-3, 3):
        ts = []
        vals = []

        for line in lines:
            vec = line.split('\t')
            vec = vec[3:]
            ts.append(float(vec[0])*1e12)
            vals.append(float(vec[i]))

        text = str(x_vals[indexer]) + 'nm'
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


    plotname = type + '/' + modules_folder + ani + '/plots/plateau/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/plateau_V'  + str(V) + '_damping' + str(damping) + '_' + str(T) + 'K.png'
    fig.suptitle('<mxdmdt> over time')
    fig.tight_layout()
    fig.savefig(plotname, dpi=600)
    # plt.show()

def plot_tAvg_SA(meshdims, cellsize, t, V, damping, MEC, ani, T, type, x_start, x_stop):

    plt.figure(figsize=(10, 7))

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += '+mec'
    modules_folder += '/'

    folder_name = type + '/' + modules_folder + ani + '/plots/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = type + '/' + modules_folder + ani + '/cache/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.txt'
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

    plt.plot(xs, ys, linewidth=2)
    plt.xlabel(r'$x$ $(\mu m)$')
    plt.ylabel(r'$\mu_x$')
    plt.tight_layout()

    plotname = type + '/' + modules_folder + ani + '/plots/' + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(T) + 'K.png'
    plt.savefig(plotname, dpi=600)
    plt.show()

def plot_tAvg_SA_2D(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani):

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

def plot_tAvg_SA_z(meshdims, cellsize, t, V, damping, data, x_start, x_stop, MEC, ani):

    mec_folder = ''
    if MEC:
        mec_folder = 'MEC/'

    folder_name = ani + '/plots/' + mec_folder + 't_avg/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    filename = ani + '/cache/' + mec_folder + 't_avg/'  + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '_zdir.txt'
    filename_2D = ani + '/cache/' + mec_folder + 't_avg/'  + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/2D_tAvg_damping' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '.txt'

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

        xs = np.linspace(0, len(ys), len(ys))

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
        tAvg_data = np.zeros_like(raw_data[0]) # Haven't tried this before, but should work, right?

        for k, matrix in enumerate(raw_data):
            for i, row in enumerate(matrix):
                for j, col in enumerate(row):
                    tAvg_data[i][j] += col
                    if k == len(raw_data)-1:
                        tAvg_data[i][j] /= len(raw_data)

        ys = np.array([tAvg_data[i][int(len(tAvg_data[0])/2)] for i in range(len(tAvg_data))])
        xs = np.linspace(0, meshdims[2], len(ys))

    else:
        print("No simulations have been done with these params")
        exit()

    plt.plot(xs, np.flip(ys))
    plt.xlabel("z (nm)")
    plt.ylabel(r'$\mu_x$')
    plt.tight_layout()

    plotname = ani + '/plots/' + mec_folder + 't_avg/' +  str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) + '/tAvg_zdir_damping_' + str(damping) + '_V' + str(V) + '_' + str(data[1:-1]) + '_t' + str(t) + 'ps.png'
    plt.savefig(plotname, dpi=500)
    plt.show()

def plot_tAvg_comparison(plots, legends, savename, ani):

    indexer = 0
    plt.figure(figsize=(13,8))


    N = len(plots)
    colors = plt.cm.viridis(np.linspace(0,1,N+1))
    # colors = ['tab:blue', 'tab:red']


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
        if k > 0:
            background = 4*sum([y for y in ys[3*int(len(ys)/4):]]) / len(ys)
            ys = [y - background for y in ys]

        leg = round(ys[0] * 1e-11, 1)

        ys = [(y )/ys[0] for y in ys]
        xs = np.linspace(2.02, 4, len(ys))
        # xs = np.linspace(2.02, 4, len(ys))

        

        plt.plot(xs, ys, label=legends[k], linewidth=3, color=colors[k])
        # plt.plot(xs, ys, label=str(leg), linewidth=3, color=colors[k])

        indexer += 1

    plt.xlabel('x (μm)')
    plt.ylabel(r'$\mu$ (normalized)')

    # plt.legend(title= r'$\mu$ value underneath injector ($10^{11}$)')
    plt.legend(title='Temperature (K)')

    plt.savefig(savename, dpi=600)
    plt.show()

def plot_magnon_dispersion(meshdims, cellsize, t, V, damping, MEC, ani, T, type, dir, axis, transport = False):

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += 'mec'
    modules_folder += '/'

    time_step = 0.1e-12


    if transport:
        output_file = path + type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_axis' + axis + 'V' + str(V) + '_damping' + str(damping) + '_T' + str(T) +  '_dispersion.txt'
        savename = type + '/' + modules_folder + ani + '/plots/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_axis' + axis + 'V' + str(V) + '_damping' + str(damping) + '_T' + str(T) +  '_dispersion.png'
    else:
        output_file = path + type + '/' + modules_folder + ani + '/cache/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_axis' + axis + 'groundstate' + '_damping' + str(damping) + '_T' + str(T) +  '_dispersion.txt'
        savename = type + '/' + modules_folder + ani + '/plots/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/' + 'dir' + dir + '_axis' + axis + 'groundstate' + '_damping' + str(damping) + '_T' + str(T) +  '_dispersion.png'

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

    ax1.imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, 400))

    label = 'q' + dir

    ax1.set_xlabel(label)
    ax1.set_ylabel('f (THz)')
    # ax1.set_ylim(0, 0.1)

    plt.tight_layout()

    folder_name = type + '/' + modules_folder + ani + '/plots/dispersions/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)


    plt.savefig(savename, dpi=600)

    plt.show()

def plot_phonon_dispersion(meshdims, damping, MEC, ani, dir,time_step):

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
    titles = [r'$\mu_0$ = 5.0 $\cdot 10^{11}$', r'$\mu_0$ = 14.8 $\cdot 10^{11}$', r'$\mu_0$ = 29.6 $\cdot 10^{11}$']

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

        clim_max = 1000
        # if i == 0:
        #     clim_max = 1000

        label = 'qy'

        ax.imshow(result, origin='lower', interpolation='bilinear', extent = [-k_max, k_max,f_min, f_max], aspect ="auto", clim=(0, clim_max))

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

def plot_critical_T(meshdims, damping, MEC, ani, type):

    plt.figure(figsize=(10,6))

    modules_folder = 'ex+ani'
    if MEC:
        modules_folder += 'mec'
    modules_folder += '/'

    output_file = path + type + '/' + modules_folder + ani + '/cache/critical_T/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/critical_T.txt'

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
        

    max_y = max(ys)
    ys = [y/max_y for y in ys]

    plt.plot(xs, ys, linewidth=3)

    folder_name = type + '/' + modules_folder + ani + '/plots/critical_T/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    savename = type + '/' + modules_folder + ani + '/plots/critical_T/' + str(meshdims[0]) + 'x' + str(meshdims[1]) + 'x' + str(meshdims[2]) +  '/critical_T.png'
    plt.xlabel(r'Temperature ($K$)')
    plt.ylabel(r'${m}_{x}$')
    plt.tight_layout()
    plt.savefig(savename, dpi=600)
    plt.show()

def plot_diffusion_length(plots, ts, savename, ani):
    
    plt.figure(figsize=(10,7))


    N = len(ts)
    colors = plt.cm.viridis(np.linspace(0,1,N+1))
    # colors = ['tab:blue', 'tab:red']


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
        #     ys = ys[int(len(ys)/2):int(3*len(ys)/4)]
        #     xs = np.linspace(2.02, 3, len(ys))
        # else:
        #     ys = ys[int(len(ys)/2):]
        #     xs = np.linspace(2.02, 2.5, len(ys))

        if k == 0:
            ys = ys[int(len(ys)/2):]

        ys = [(((y )/ys[0])) for y in ys]
        xs = np.linspace(0, 1.98, len(ys))

        def f(x, a, b):
            return a + np.exp(b* x**2)

        params, params_cov = curve_fit(f, xs, ys)

        # plt.plot(ts[k], params[1], marker='v', markersize = 10, color='black')
        plt.plot(f(xs, params[0], params[1]), xs, color=colors[k])
        plt.plot(xs, ys, color=colors[k])

    plt.xlabel('Thickness (layers)')
    plt.ylabel(r'$L_d$ (μm)')

    plt.tight_layout()
    plt.savefig(savename, dpi=600)

    plt.show()

def main():

    # f1 = 'IP/cache/t_avg/4000x50x5/tAvg_damping0.0004_V-0.0115_mxdmdt.txt'
    # f2 = 'IP/cache/t_avg/4000x50x10/tAvg_damping0.0004_V-0.045_mxdmdt.txt'
    # f3 = 'IP/cache/t_avg/4000x50x15/tAvg_damping0.0004_V-0.12_mxdmdt.txt'
    # f4 = 'IP/cache/t_avg/4000x50x20/tAvg_damping0.0004_V-0.27_mxdmdt.txt'
    # f5 = 'IP/cache/t_avg/4000x50x25/tAvg_damping0.0004_V-0.4_mxdmdt.txt'
    # f6 = 'IP/cache/t_avg/4000x50x30/tAvg_damping0.0004_V-0.53_mxdmdt.txt'
    # f7 = 'IP/cache/t_avg/4000x50x35/tAvg_damping0.0004_V-0.63_mxdmdt.txt'
    # f8 = 'IP/cache/t_avg/4000x50x40/tAvg_damping0.0004_V-0.73_mxdmdt.txt'
    # f9 = 'IP/cache/t_avg/4000x50x45/tAvg_damping0.0004_V-0.8_mxdmdt.txt'

    # l1 = '1'
    # l2 = '2'
    # l3 = '3'
    # l4 = '4'
    # l5 = '5'
    # l6 = '6'
    # l7 = '7'
    # l8 = '8'
    # l9 = '9'

    # # # # # # title = 'Normalized spin accumulation with/without MEC'

    # f1 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x25/tAvg_damping0.0004_V-0.4_0.3K.txt'
    # f2 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x25/tAvg_damping0.0004_V-0.5_5K.txt'
    # f3 = 'AFM/ex+ani/IP/cache/t_avg/4000x50x5/tAvg_damping0.0004_V-0.04_0.3K.txt'

    # l1 = '0.3'
    # l2 = '5'
    # l3 = '5'

    # savename = 'AFM/ex+ani/IP/plots/custom/5layer_temp_comparison.png'

    # plot_tAvg_comparison([f1,f2], [l1,l2], savename, 'IP')

    # plot_dispersion([4000, 50, 5], 4e-4, 1, 'OOP', 'y')

    # # # FOR DISPERSIONS DOWN HERE

    f1 = 'AFM/ex+ani/IP/cache/dispersions/4000x50x5/diry_axisxV-0.0115_damping0.0004_T0.3_dispersion.txt'
    f2 = 'AFM/ex+ani/IP/cache/dispersions/4000x50x5/diry_axisxV-0.02_damping0.0004_T0.3_dispersion.txt'
    f3 = 'AFM/ex+ani/IP/cache/dispersions/4000x50x5/diry_axisxV-0.04_damping0.0004_T0.3_dispersion.txt'
    # f4 = 'IP/cache/dispersions/1000x50x150/steady/diry_axisx_dispersion.txt'
    # f5 = 'IP/cache/dispersions/1000x50x190/steady/diry_axisx_dispersion.txt'


    savename = 'AFM/ex+ani/IP/plots/custom/1layer_dispersion_torque_comparison_0.3K.png'

    plot_dispersions((f1,f2,f3), savename)

    # ### DIFFUSION LENGTHS DOWN HERE

    # f1 = 'IP/cache/t_avg/4000x50x5/tAvg_damping0.0004_V-0.0115_mxdmdt.txt'
    # f2 = 'IP/cache/t_avg/4000x50x10/tAvg_damping0.0004_V-0.045_mxdmdt.txt'
    # f3 = 'IP/cache/t_avg/4000x50x15/tAvg_damping0.0004_V-0.12_mxdmdt.txt'
    # f4 = 'IP/cache/t_avg/4000x50x20/tAvg_damping0.0004_V-0.25_mxdmdt.txt'
    # f5 = 'IP/cache/t_avg/4000x50x25/tAvg_damping0.0004_V-0.4_mxdmdt.txt'
    # f6 = 'IP/cache/t_avg/4000x50x30/tAvg_damping0.0004_V-0.53_mxdmdt.txt'
    # f7 = 'IP/cache/t_avg/4000x50x35/tAvg_damping0.0004_V-0.63_mxdmdt.txt'
    # f8 = 'IP/cache/t_avg/4000x50x40/tAvg_damping0.0004_V-0.73_mxdmdt.txt'
    # f9 = 'IP/cache/t_avg/4000x50x45/tAvg_damping0.0004_V-0.8_mxdmdt.txt'

    # ts = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    # savename = 'IP/plots/t_avg/finished/diffusion_length.png'

    # plot_diffusion_length([f1,f2,f3,f4,f5,f6,f7,f8,f9], ts, savename, 'IP')


if __name__ == '__main__':
    main()