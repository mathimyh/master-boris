import os
from pathlib import Path

path = str(Path(__file__).resolve().parent) + '/'

def make_folder(filepath):
    folder_name = '/'.join(filepath.split('/')[:-1])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

def check_file_exists(filepath):
    try:
        open(filepath, 'r')
    except FileNotFoundError:
        print('No simulations have been done with these parameters.')

class Params:
    
    def __init__(self, meshdims, cellsize, t, 
                    V, damping, MEC, ani, T, type,
                    hard_axis, Hfield):
        
        self.meshdims = meshdims
        self.cellsize = cellsize
        self.t = t
        self.V = round(V, 3)
        self.damping = damping
        self.MEC = MEC
        self.ani = ani
        self.T = T
        self.type = type
        self.hard_axis = hard_axis
        self.Hfield = Hfield

        self.modules_folder = 'ex+ani'
        if MEC:
            self.modules_folder += '+mec'
        if hard_axis:
            self.modules_folder += '+hard_axis'
        self.H_string = ''
        if Hfield != 0.:
            self.modules_folder += '+Hfield'
            self.H_string += '_H' + str(round(Hfield,1))   
        self.modules_folder += '/'
    
class TimeAvgSA(Params):

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield, x_start, x_stop, direction = 'x'):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield)
        self.direction = direction
        self.x_start = x_start
        self.x_stop = x_stop
        
        self.dir_str = ''
        if direction != 'x':
            self.dir_str = direction + '_'

    def simname(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/sims/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2])  + f"/V{self.V:.3f}" + '_damping' + str(self.damping) + '_' + str(self.T) + 'K_steady_state.bsm'

    def cachename(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/cache/' + 't_avg/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/' + self.dir_str + 'tAvg_damping' + str(self.damping) + f"_V{self.V:.3f}" + '_' + str(self.T) + 'K.txt'
    
    def plotname(self):
        return self.type + '/' + self.modules_folder + self.ani + '/plots/' + 't_avg/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/' + self.dir_str + 'tAvg_damping' + str(self.damping) + f"_V{self.V:.3f}" + '_' + str(self.T) + 'K.png'

    def set_hard_axis(self, axis_bool, hfield):
        self.hard_axis = axis_bool
        self.Hfield = hfield
        self.modules_folder = 'ex+ani'
        if self.MEC:
            self.modules_folder += '+mec'
        if self.hard_axis:
            self.modules_folder += '+hard_axis'
        self.H_string = ''
        if self.Hfield != 0.:
            self.modules_folder += '+Hfield'
            self.H_string += '_H' + str(round(self.Hfield,1))   
        self.modules_folder += '/'

class Steadystate(Params):

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield, x_vals=False):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield)
        self.x_vals = x_vals

    def simname(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/sims/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + f"/V{self.V:.3f}" + '_damping' + str(self.damping) + '_' + str(self.T) + 'K_steady_state.bsm'

    def cachename(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/cache/plateau/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/plateau_' + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_' + str(self.T) + 'K.txt'
    
    def plotname(self):
        return self.type + '/' + self.modules_folder + self.ani + '/plots/plateau/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/plateau' + str(self.damping) + f"_V{self.V:.3f}" + '_' + str(self.T) + 'K.png'
    
class MagnonDispersion(Params):

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, 
                 Hfield, component, axis, steadystate=False, triple=False):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield)
        self.component = component
        self.axis = axis
        self.steadystate = steadystate
        self.triple = triple

    def simname(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/sims/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + f"/V{self.V:.3f}" + '_damping' + str(self.damping) + '_' + str(self.T) + 'K_steady_state.bsm'

    def cachename(self):
        if not self.steadystate:
            if self.hard_axis:
                y = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/ycomponent_axis' + self.axis + 'groundstate' + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + '_dispersion.txt'
                z = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/zcomponent_axis' + self.axis + 'groundstate' + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string +  '_dispersion.txt'
                return [y,z]
            else:
                return [path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/dir' + self.component + '_axis' + self.axis + 'groundstate' + '_damping' + str(self.damping) + '_T' + str(self.T) +  '_dispersion.txt']
        else:
            if self.triple:
                # Get values from the edges and middle of the systems(in y-direction)
                n1, n2, n3 = self.cellsize, self.meshdims[1]/2, self.meshdims[1] - self.cellsize
                if self.hard_axis:
                    y1 = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/ycomponent_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + f'y={n1}_dispersion.txt'
                    z1 = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/zcomponent_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string +  f'y={n1}_dispersion.txt'
                    y2 = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/ycomponent_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + f'y={n2}_dispersion.txt'
                    z2 = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/zcomponent_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string +  f'y={n2}_dispersion.txt'
                    y3 = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/ycomponent_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + f'y={n3}_dispersion.txt'
                    z3 = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/zcomponent_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string +  f'y={n3}_dispersion.txt'
                    return [y1,z1,y2,z2,y3,z3]
                else:
                    one = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + f'_y={n1}_dispersion.txt'
                    two = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + f'_y={n2}_dispersion.txt'
                    three = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + f'_y={n3}_dispersion.txt'
                    return [one,two,three]
            else:
                if self.hard_axis:
                    y = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/ycomponent_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + '_dispersion.txt'
                    z = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/zcomponent_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string +  '_dispersion.txt'
                    return [y,z]
                else:
                    one = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + '_dispersion.txt'
                    return [one]

    def plotname(self):
        if not self.steadystate:
            return path + self.type + '/' + self.modules_folder + self.ani + '/plots/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + 'groundstate' + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + '_dispersion.png'
        else:
            return path + self.type + '/' + self.modules_folder + self.ani + '/plots/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + f"_V{self.V:.3f}" + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + '_dispersion.png'

class CriticalT(Params):

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield, max_T):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, Hfield,hard_axis)
        self.max_T = max_T

    def cachename(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/cache/critical_T/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/critical_T.txt'
    
    def plotname(self):
        return self.type + '/' + self.modules_folder + self.ani + '/plots/critical_T/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/critical_T.png'
    
class CurrentDensity(Params):
    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield)

    def cachename(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/cache/current_density/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/Jc_V'  + str(self.V) + '_damping' + str(self.damping) + '_' + str(self.T) + 'K.txt'