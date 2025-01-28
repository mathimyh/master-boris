import plotting
import transport
import dispersion_relations

def main():
    
    # Dimensions
    Lx = 4000
    Ly = 500
    Lz = 25
    cellsize = 5
    meshdims = (Lx, Ly, Lz)

    # Parameters
    t = 1000
    V = -0.4
    data = '<mxdmdt>'
    damping = 4e-4
    MEC = 0
    ani = 'IP'
    type = 'AFM'
    T = 0.3
    x_vals = [2020, 2300, 2600, 3000, 3500, 4000]

    params = [meshdims, cellsize, t, V, damping, MEC, ani, T, type]

    transport.save_steadystate(*params, x_vals)
    # transport.time_avg_SA(*params, 2020, 4000)
    # plotting.plot_plateau(meshdims, V, damping, x_vals, MEC, ani, type)


if __name__ == '__main__':
    main()