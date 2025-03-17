from NetSocks import NSClient # type:ignore

# Dimensions
Lx = 4000
Ly = 50
Lz = 5
cellsize = 5
meshdims = (Lx, Ly, Lz)

# Parameters
t = 100
V = -0.021
damping = 4e-4
MEC = 0
ani = 'IP'
type = 'FM'
T = 0.3

ns = NSClient(); ns.configure(True, False)
ns.cuda(1)
ns.reset()
ns.iterupdate(200)

modules = ['exchange', 'aniuni']

# Set up the antiferromagnet
FM = ns.Ferromagnet(np.array(meshdims)*1e-9, [cellsize*1e-9])
FM.modules(modules)
temp = str(T) + 'K'
ns.temperature(temp)
ns.setode('sLLG', 'RK4') # Stochastic LLG for temperature effects
ns.setdt(1e-15) 

# Set parameters. Should possibly just save these in the database really    
FM.param.grel = 1
FM.param.damping =  damping
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
if ani == 'OOP':
    FM.param.ea1 = (0,0,1) # Set it z-direction
    FM.setangle(0,90) # Just move m to z-direction, not necessary to wait every time
elif ani == 'IP':
    FM.param.ea1 = (1,0,0)
else:
    print('Choose anisotropy direction')
    exit()

# Add increased damping at edges along x-axis to prevent reflections
FM.param.damping.setparamvar('abl_tanh', [300/meshdims[0], 300/meshdims[0], 0, 0, 0, 0, 1, 10, 300]) 

# Add the magnetoelastic parameters if necessary
if MEC:
    FM.surfacefix('-z')
    ns.seteldt(1e-15)
    FM.mcellsize([cellsize*1e-9]) 
    FM.param.cC = (36e10, 17e10, 8.86e10) # N/m^2       A. Yu. Lebedev et al (1989)
    FM.param.density = 5250 #kg/m^3       found this on google
    FM.param.MEc = (-3.44e6, 7.5e6) #J/m^3   G. Wedler et al (1999) 
    FM.param.mdamping = 1e15 # Should probably be lower than this 

# FM.pbc(['x', 10])

# ns.random()

# Relax for 1000 ps to thermally excite magnons
ns.Relax(['time', 1000e-12])