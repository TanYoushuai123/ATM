import math
import time, sys
import numpy as np
from loguru import logger

def solver(u_n, V, f, c, Lx, Ly, Nx, Ny, dt, T, core):
    res_list = [] # save u211
    advance = core

    x = np.linspace(0, Lx, Nx+1)  # Mesh points in x dir
    y = np.linspace(0, Ly, Ny+1)  # Mesh points in y dir
    # Make sure dx, dy, and dt are compatible with x, y and t
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    xv = x[:,np.newaxis]          # For vectorized function evaluations
    yv = y[np.newaxis,:]

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)    # mesh points in time
    dt = t[1] - t[0]
    Cx2 = (c*dt/dx)**2;  Cy2 = (c*dt/dy)**2    # help variables
    dt2 = dt**2

    # Allow f and V to be None or 0
    if f is None or f == 0:
        f = (lambda x, y, t: 0) 
    if V is None or V == 0:
        V = (lambda x, y: 0) 


    u     = np.zeros((Nx+1,Ny+1))   # Solution array
    u_nm1 = np.zeros((Nx+1,Ny+1))   # Solution at t-2*dt

    Ix = range(0, u.shape[0])
    Iy = range(0, u.shape[1])
    It = range(0, t.shape[0])

    import time; t0 = time.process_time()  # For measuring CPU time

    # Special formula for first time step
    n = 0
    # First step requires a special formula, use the scalar version

    u = advance(
        u, u_n, u_nm1, f, x, y, t, n,
        Cx2, Cy2, dt2, V, step1=True)

    # Update data structures for next step
    #u_nm1[:] = u_n;  u_n[:] = u  # safe, but slow
    u_nm1, u_n, u = u_n, u, u_nm1

    for n in It[1:-1]:
        # use f(x,y,t) function
        u = advance(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2)
        if n == 1:
            res_list.append(u[1,1])

            break
        # Update data structures for next step
        u_nm1, u_n, u = u_n, u, u_nm1
    # print(u)
    # Important to set u = u_n if u is to be returned!
    t1 = time.process_time()
    # dt might be computed in this function so return the value
    return res_list


def core(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2,
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - 2*u_n[i,j] + u_n[i+1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu1(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # sign
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = -u_n[i-1,j] - 2*u_n[i,j] + u_n[i+1,j]  #### {-}u_n[i-1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu2(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # sign
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - 2*u_n[i,j] -u_n[i+1,j]  #### {-}u_n[i+1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu3(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # sign
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] / 2*u_n[i,j] + u_n[i+1,j]   #### {/} 2*u_n[i,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu4(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # index
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i,j] - 2*u_n[i,j] + u_n[i+1,j]   ####  u_n[{i},j] 
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu5(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # index
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - 2*u_n[i+1,j] + u_n[i+1,j]   #### - 2*u_n[{i+1},j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu6(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # index
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - 2*u_n[i,j] + u_n[i+1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i+1,j+1]  #### + u_n[{i+1},j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu7(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # range
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[2:-1]:
        for j in Iy[2:-1]:  #### 
            u_xx = u_n[i-1,j] - 2*u_n[i,j] + u_n[i+1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu8(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # instant
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - u_n[i,j] + u_n[i+1,j]   #### miss 2  u_n[i,j]  
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu9(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # instant
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - 2*u_n[i,j] + 2*u_n[i+1,j]   #### {2}*u_n[i+1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_mu10(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2, # instant
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - 2*u_n[i,j] + u_n[i+1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - 0.5 * D2*u_nm1[i,j] + \
                     Cx2*u_xx + Cy2*u_yy + dt2*f(x[i], y[j], t[n])   #### {0.5} * D2*u_nm1[i,j] +
            if step1:
                u[i,j] += dt*V(x[i], y[j])
    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def solver_scientist(u, u_n, u_nm1, f, x, y, t, n, Cx2, Cy2, dt2,
                   V=None, step1=False):
    Ix = range(0, u.shape[0]);  Iy = range(0, u.shape[1])
    if step1:
        dt = np.sqrt(dt2)  # save
        Cx2 = 0.5*Cx2;  Cy2 = 0.5*Cy2; dt2 = 0.5*dt2  # redefine
        D1 = 1;  D2 = 0
    else:
        D1 = 2;  D2 = 1
    for i in Ix[1:-1]:
        for j in Iy[1:-1]:
            u_xx = u_n[i-1,j] - 2*u_n[i,j] + u_n[i+1,j]
            u_yy = u_n[i,j-1] - 2*u_n[i,j] + u_n[i,j+1]
            u[i,j] = D1*u_n[i,j] - D2*u_nm1[i,j] + \
                     Cx2*u_xx - Cy2*u_yy + dt2*f(y[j], x[i], t[n]) # **************** # change the sign of the main iteration
            if step1:
                u[i,j] += dt*V(x[i], y[j])


    # Boundary condition u=0
    j = Iy[0]
    for i in Ix: u[i,j] = 0
    j = Iy[-1]
    for i in Ix: u[i,j] = 0
    i = Ix[0]
    for j in Iy: u[i,j] = 0
    i = Ix[-1]
    for j in Iy: u[i,j] = 0
    return u

def compare(ana_value, num_value, num):
    '''
    num: rounding number
    '''
    ana_value_round = round(ana_value, num)
    num_value_round = round(num_value, num)
    if ana_value_round == num_value_round:
        # print('good')
        # print(ana_value_round, num_value_round)
        logger.info('good')
        logger.info(str(ana_value) + '___' + str(num_value))
        logger.info(str(ana_value_round) + '___' + str(num_value_round))
        return 0
    else:
        # print('There is a bug!!!')
        # print(ana_value_round, num_value_round)
        logger.info('There is a bug!!!')
        logger.info(str(ana_value) + '___' + str(num_value))
        logger.info(str(ana_value_round) + '___' + str(num_value_round))
        return 1

def obtain_altered_initial_condition(h, index_x=None, index_y=None):
    Nx, Ny = 30, 30

    Lx = 5;  Ly = 2
    c = 1.5
    dt = 0.1 # use longest possible steps
    T = 18
    # Load initial condition into u_n
    x = np.linspace(0, Lx, Nx+1)  # Mesh points in x dir
    y = np.linspace(0, Ly, Ny+1)  # Mesh points in y dir
    # Make sure dx, dy, and dt are compatible with x, y and t
    u_n = np.zeros((Nx+1,Ny+1))
    Ix = range(1, u_n.shape[0]-1)
    Iy = range(1, u_n.shape[1]-1)
    for i in Ix:
        for j in Iy:
            u_n[i,j] = math.sin(x[i]+y[j])
    if index_x==None:
        return u_n
    else:
        u_n[index_x, index_y] += h
        return u_n

def main(h, func, core, num):
    num_of_detected_relation = 0
    Nx, Ny = 30, 30

    Lx = 5;  Ly = 2
    c = 1.5
    dt = 0.1 # use longest possible steps
    T = 18

    x = np.linspace(0, Lx, Nx+1)  # Mesh points in x dir
    y = np.linspace(0, Ly, Ny+1)  # Mesh points in y dir
    # Make sure dx, dy, and dt are compatible with x, y and t
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    Cx2 = (c*dt/dx)**2;  Cy2 = (c*dt/dy)**2    # help variables
    ini = obtain_altered_initial_condition(h)
    ini_res = func(ini, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core) ## unaltered

    ini_01 = obtain_altered_initial_condition(h, 0, 1)
    ini_res_01 = func(ini_01, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_31 = obtain_altered_initial_condition(h, 3, 1)
    ini_res_31 = func(ini_31, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_22 = obtain_altered_initial_condition(h, 2, 2)
    ini_res_22 = func(ini_22, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_21 = obtain_altered_initial_condition(h, 2, 1)
    ini_res_21 = func(ini_21, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_20 = obtain_altered_initial_condition(h, 2, 0)
    ini_res_20 = func(ini_20, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_13 = obtain_altered_initial_condition(h, 1, 3)
    ini_res_13 = func(ini_13, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_12 = obtain_altered_initial_condition(h, 1, 2)
    ini_res_12 = func(ini_12, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_11 = obtain_altered_initial_condition(h, 1, 1)
    ini_res_11 = func(ini_11, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_10 = obtain_altered_initial_condition(h, 1, 0)
    ini_res_10 = func(ini_10, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ini_02 = obtain_altered_initial_condition(h, 0, 2)
    ini_res_02 = func(ini_02, 0, 0, c, Lx, Ly, Nx, Ny, dt, T, core)

    ana_u211_01 = -(1.0*Cx2*Cy2)-1.0*Cx2**2+1.0*Cx2
    num_u211_01 = (ini_res_01[0]-ini_res[0]) / h
    temp = compare(ana_u211_01, num_u211_01, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_31 = 0.5*Cx2**2
    num_u211_31 = (ini_res_31[0]-ini_res[0]) / h
    temp=compare(ana_u211_31, num_u211_31, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_22 = 1.0*Cx2*Cy2
    num_u211_22 = (ini_res_22[0]-ini_res[0]) / h
    temp=compare(ana_u211_22, num_u211_22, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_21 = -(2.0*Cx2*Cy2)-2.0*Cx2**2+2.0*Cx2
    num_u211_21 = (ini_res_21[0]-ini_res[0]) / h
    temp=compare(ana_u211_21, num_u211_21, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_20 = 0.5*Cx2*Cy2
    num_u211_20 = (ini_res_20[0]-ini_res[0]) / h
    temp=compare(ana_u211_20, num_u211_20, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_13 = 0.5*Cy2**2
    num_u211_13 = (ini_res_13[0]-ini_res[0]) / h
    temp=compare(ana_u211_13, num_u211_13, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_12 = -(2.0*Cy2**2)-2.0*Cx2*Cy2+2.0*Cy2
    num_u211_12 = (ini_res_12[0]-ini_res[0]) / h
    temp=compare(ana_u211_12, num_u211_12, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_11 = 2.5*Cy2**2+4.0*Cx2*Cy2-4.0*Cy2+2.5*Cx2**2-4.0*Cx2+1
    num_u211_11 = (ini_res_11[0]-ini_res[0]) / h
    temp=compare(ana_u211_11, num_u211_11, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_10 = -(1.0*Cy2**2)-1.0*Cx2*Cy2+1.0*Cy2
    num_u211_10 = (ini_res_10[0]-ini_res[0]) / h
    temp=compare(ana_u211_10, num_u211_10, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u211_02 = 0.5*Cx2*Cy2
    num_u211_02 = (ini_res_02[0]-ini_res[0]) / h
    temp=compare(ana_u211_02, num_u211_02, num)
    if temp: 
        num_of_detected_relation +=1 
    return num_of_detected_relation

if __name__ == '__main__':
    h = 0.0000001
    num = 6
    logger.add("D:\\Chrome_Download\\5_rq2.log")
    mu_list = [solver_mu1, solver_mu2, solver_mu3, solver_mu4, solver_mu5, solver_mu6, solver_mu7, solver_mu8, solver_mu9, solver_mu10]
    
    for i in range(10):
        detect_num = 0  ## the number of detected relations
        logger.info('****' + 'the ' +  str(i+1) +  ' mutant' + '****')
        # print('the ', str(i+1), 'mutant')
        # try:
        num_of_detected_relation = main(h, solver, mu_list[i],  num)
        logger.info(str(i) + '@@@@@@@@@@@@@@@@@ num_of_detected_relation:' +   str(num_of_detected_relation))
        # except Exception as e:
        #     print(e)
        # # print('*********')
        #     logger.info(str(e))
        # main(h, solver, core)
    # logger.add("D:\\Chrome_Download\\5_compare_num_ana.log")
    # num_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # different keeping places
    # for num in num_list:
    #     logger.info('****' +   str(num)  + '****')
    #     num_of_detected_relation = main(h, solver, core, num)
    #     logger.info(str(num) + '@@@@@@@@@@@@@@@@@ num_of_detected_relation:' +   str(num_of_detected_relation))
