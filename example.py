import sys
import os
sys.path.insert(0, 'C:/users/mathimyh/documents/boris data/borispythonscripts/')

from NetSocks import NSClient   # type: ignore
import numpy as np

ns = NSClient(); ns.configure(True, False)
ns.cuda(1); ns.reset(); ns.clearelectrodes()

# 5nm = -0.03
# 100nm =-2.9 

# Dimensions (nm)
Lx = 1500
Ly = 50
Lz = 5
cellsize = 5
meshdims = (Lx, Ly, Lz)

# Parameters
t = 1000 # ps
V = -0.03 # mV
damping = 4e-4
T = 0.3 # K

# Set up the antiferromagnet
AFM = ns.AntiFerromagnet(np.array(meshdims)*1e-9, [cellsize*1e-9])
AFM.modules(['exchange', 'aniuni', 'SOTfield', 'transport'])
temp = str(T) + 'K'
ns.temperature(temp)
ns.setode('sLLG', 'RK4') # Stochastic LLG for temperature effects
ns.setdt(1e-15)

# Set parameters   
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
AFM.param.ea1 = (1,0,0)

# Add increased damping at edges along x-axis to prevent reflections. For smaller meshes this needs to be shortened
damping_x = 50
AFM.param.damping_AFM.setparamvar('abl_tanh', [damping_x/meshdims[0], damping_x/meshdims[0], 0, 0, 0, 0, 1, 10, damping_x])

# Set spesific params for torque
AFM.param.SHA = 1
AFM.param.flST = 1

# Current along y-direction
ns.addelectrode(np.array([(meshdims[0]/2 - 100), 0, (meshdims[2]-cellsize), (meshdims[0]/2 + 100), 0, meshdims[2]])* 1e-9)
ns.addelectrode(np.array([(meshdims[0]/2 - 100), meshdims[1], (meshdims[2]-cellsize), (meshdims[0]/2 + 100), meshdims[1], meshdims[2]]) * 1e-9)
ns.designateground('1')

# Add step function so that torque only acts on region in the injector
if meshdims[2]==5:
    width = 20
elif meshdims[2] == 100:
    width = 40
func = '(step(x-' + str(meshdims[0]/2 - width/2) + 'e-9)-step(x-' + str(meshdims[0]/2 + width/2) + 'e-9)) * (step(z-' + str(meshdims[2]-cellsize) + 'e-9)-step(z-' + str(meshdims[2]) + 'e-9))'
AFM.param.SHA.setparamvar('equation', func)
AFM.param.flST.setparamvar('equation',func)

# ns.Relax('time', 1000e-12)
ns.V([0.001*V, 'time', t*1e-12])