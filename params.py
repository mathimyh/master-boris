path = 'C:/Users/mathimyh/master/master-boris/'
import os

def make_folder(filepath):
    folder_name = '/'.join(filepath.split('/')[:-1])
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

class Params:
    
    def __init__(self, meshdims, cellsize, t, 
                    V, damping, MEC, ani, T, type,
                    hard_axis, Hfield):
        
        self.meshdims = meshdims
        self.cellsize = cellsize
        self.t = t
        self.V = V
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

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield, x_start, x_stop):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield)
        self.x_start = x_start
        self.x_stop = x_stop

    def simname(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/sims/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/V' + str(self.V) + '_damping' + str(self.damping) + '_' + str(self.T) + 'K_steady_state.bsm'

    def cachename(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/cache/' + 't_avg/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/tAvg_damping' + str(self.damping) + '_V' + str(self.V) + '_' + str(self.T) + 'K.txt'
    
    def plotname(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/plots/' + 't_avg/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/tAvg_damping' + str(self.damping) + '_V' + str(self.V) + '_' + str(self.T) + 'K.png'
    
class Steadystate(Params):

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield, x_vals=False):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield)
        self.x_vals = x_vals

    def simname(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/sims/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/V' + str(self.V) + '_damping' + str(self.damping) + '_' + str(self.T) + 'K_steady_state.bsm'

    def cachename(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/cache/plateau/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/plateau_V'  + str(self.V) + '_damping' + str(self.damping) + '_' + str(self.T) + 'K.txt'
    
    def plotname(self):
        return self.type + '/' + self.modules_folder + self.ani + '/plots/' + 't_avg/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/tAvg_damping' + str(self.damping) + '_V' + str(self.V) + '_' + str(self.T) + 'K.png'
    
class MagnonDispersion(Params):

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, 
                 Hfield, component, axis, steadystate=False, triple=False):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield)
        self.component = component
        self.axis = axis
        self.steadystate = steadystate
        self.triple = triple

    def simname(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/sims/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) + '/V' + str(self.V) + '_damping' + str(self.damping) + '_' + str(self.T) + 'K_steady_state.bsm'

    def cachename(self):
        if not self.steadystate:
            if self.hard_axis:
                y = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/ycomponent_axis' + self.axis + 'groundstate' + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + '_dispersion.txt'
                z = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/zcomponent_axis' + self.axis + 'groundstate' + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string +  '_dispersion.txt'
                return [y,z]
            else:
                return [path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/dir' + self.component + '_axis' + self.axis + 'groundstate' + '_damping' + str(self.damping) + '_T' + str(self.T) +  '_dispersion.txt']
        else:
            # Fix this for hard_axis
            if self.triple:
                one = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + 'V' + str(self.V) + '_damping' + str(self.damping) + '_T' + str(self.T) + '_y=25_dispersion.txt'
                two = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + 'V' + str(self.V) + '_damping' + str(self.damping) + '_T' + str(self.T) + '_y=5_dispersion.txt'
                three = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + 'V' + str(self.V) + '_damping' + str(self.damping) + '_T' + str(self.T) + '_y=45_dispersion.txt'
                return [one,two,three]
            else:
                if self.hard_axis:
                    y = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/ycomponent_axis' + self.axis + 'V' + str(self.V) + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + '_dispersion.txt'
                    z = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/zcomponent_axis' + self.axis + 'V' + str(self.V) + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string +  '_dispersion.txt'
                    return [y,z]
                else:
                    one = path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + 'V' + str(self.V) + '_damping' + str(self.damping) + '_T' + str(self.T) + '_dispersion.txt'
                    return [one]

    def plotname(self):
        if not self.steadystate:
            return path + self.type + '/' + self.modules_folder + self.ani + '/plots/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + 'groundstate' + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + '_dispersion.png'
        else:
            return path + self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/' + 'dir' + self.component + '_axis' + self.axis + 'V' + str(self.V) + '_damping' + str(self.damping) + '_T' + str(self.T) + self.H_string + '_dispersion.png'

class MagnonDispersionSinc(Params):

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield, component, axis):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield)
        self.component = component
        self.axis = axis

    def cachename(self):
        y = self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/sinc_ycomponent_axis' + self.axis + '_damping' + str(self.damping) + '_dispersion.txt'
        z = self.type + '/' + self.modules_folder + self.ani + '/cache/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/sinc_zcomponent_axis' + self.axis + '_damping' + str(self.damping) + '_dispersion.txt'
        return [y,z]

    def plotname(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/plots/dispersions/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/sinc_axis' + self.axis + '_damping' + str(self.damping) +  '_dispersion.png'

class CriticalT(Params):

    def __init__(self, meshdims, cellsize, t, V, damping, MEC, ani, T, type, hard_axis, Hfield, max_T):
        super().__init__(meshdims, cellsize, t, V, damping, MEC, ani, T, type, Hfield,hard_axis)
        self.max_T = max_T

    def cachename(self):
        return path + self.type + '/' + self.modules_folder + self.ani + '/cache/critical_T/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/critical_T.txt'
    
    def plotname(self):
        return self.type + '/' + self.modules_folder + self.ani + '/plots/critical_T/' + str(self.meshdims[0]) + 'x' + str(self.meshdims[1]) + 'x' + str(self.meshdims[2]) +  '/critical_T.png'