import math
import sys, time
import numpy as np
from loguru import logger

def f(x, t):
    return 0

def solver(u_n, a, f, L, dt, F, T):
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_scientist(u_n, a, f, L, dt, F, T):
    """
    please simulate possible bug in this function and do not change the initial value
    ##### The altered code could not run, and note by # **************** # 
    """
    res_list = []
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)


    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx+1] = 0 # **************** # change the index of the main iteration
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    return res_list  # u_n holds latest u

def solver_mu1(u_n, a, f, L, dt, F, T): #sign
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = -u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])    #### {-}u_n[i]

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu2(u_n, a, f, L, dt, F, T): # sign
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] / 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])   #### {/} 2*u_n[i]

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu3(u_n, a, f, L, dt, F, T): # sign
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] + 2*u_n[i] * u_n[i+1]) + \
                   dt*f(x[i], t[n])  #### {+} 2*u_n[i]

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu4(u_n, a, f, L, dt, F, T): # index
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i-1] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])  #### u_n[{i-1}]

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu5(u_n, a, f, L, dt, F, T): # index
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i-1] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])   #### u[{i-1}]

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu6(u_n, a, f, L, dt, F, T): # index
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i+1] + u_n[i+1]) + \
                   dt*f(x[i], t[n])  #### - 2*u_n[{i+1}]

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu7(u_n, a, f, L, dt, F, T): # range
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx-3):   ####  range(1, {Nx-3})
            u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu8(u_n, a, f, L, dt, F, T): # instant
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + 0.5 * F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])  #### {0.5} * F*(u_n[i-1] 

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu9(u_n, a, f, L, dt, F, T): # instant
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + 2 * u_n[i+1]) + \
                   dt*f(x[i], t[n])  #### {2} * u_n[i+1])

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u

def solver_mu10(u_n, a, f, L, dt, F, T): # instant
    res_list = [] # take the corners of the second step
    import time;  t0 = time.process_time()  # For measuring the CPU time

    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u   = np.zeros(Nx+1)

    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = 2 * u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + \
                   dt*f(x[i], t[n])  #### {2} * u_n[i] 

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if n==1:
            res_list.append(u[1])
            res_list.append(u[-2])

        # Switch variables before next step
        #u_n[:] = u  # safe, but slow
        u_n, u = u, u_n

    t1 = time.process_time()
    
    return res_list  # u_n holds latest u


def obtain_altered_initial_condition(h, index=None):
    a = 3.5
    L = 1.5
    Nx = 30
    F = 0.5
    # Compute dt from Nx and F
    dx = L/Nx;  dt = F/a*dx**2

    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) 
    u_n = np.zeros(Nx+1)
    for i in range(0, Nx+1):
        u_n[i] = math.sin(x[i])
    u_n[0] = 0
    u_n[Nx] = 0
    if index==None:
        return u_n
    else:
        u_n[index] += h
        return u_n

def compare(ana_value, num_value, num):
    '''
    num: rounding number
    '''
    ana_value_round = round(ana_value, num)
    num_value_round = round(num_value, num)
    if ana_value_round == num_value_round:
        # print('good')        # print(ana_value_round, num_value_round)
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

def main(h, func, num):
    num_of_detected_relation = 0
    a = 3.5
    L = 1.5
    Nx = 30
    F = 0.5
    T=2
    # Compute dt from Nx and F
    dx = L/Nx;  dt = F/a*dx**2

    dx = np.sqrt(a*dt/F)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) 
    k = F
    n=Nx
    ini = obtain_altered_initial_condition(h) ## unaltered
    res_list = func(ini, a, f, L, dt, F, T)
    ini_0 = obtain_altered_initial_condition(h, 0) ## altered 1 initial condition
    res_list_0 = func(ini_0, a, f, L, dt, F, T)
    ini_1 = obtain_altered_initial_condition(h, 1) ## altered 1 initial condition
    res_list_1 = func(ini_1, a, f, L, dt, F, T)
    ini_2 = obtain_altered_initial_condition(h, 2) ## altered 2 initial condition
    res_list_2 = func(ini_2, a, f, L, dt, F, T)
    ini_3 = obtain_altered_initial_condition(h, 3) ## altered 3 initial condition
    res_list_3 = func(ini_3, a, f, L, dt, F, T)
    ini_n_3 = obtain_altered_initial_condition(h, n-3) ## altered n-3 initial condition
    res_list_n_3 = func(ini_n_3, a, f, L, dt, F, T)
    ini_n_2 = obtain_altered_initial_condition(h, n-2) ## altered n-2 initial condition
    res_list_n_2 = func(ini_n_2, a, f, L, dt, F, T)
    ini_n_1 = obtain_altered_initial_condition(h, n-1) ## altered n-1 initial condition
    res_list_n_1 = func(ini_n_1, a, f, L, dt, F, T)
    ini_n = obtain_altered_initial_condition(n) ## altered n initial condition
    res_list_n = func(ini_n, a, f, L, dt, F, T)
    
    true_deri_u21_u00 = k-2*k**2
    num_deri_u21_u00 = (res_list_0[0] - res_list[0])/h
    temp=compare(true_deri_u21_u00, num_deri_u21_u00, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u21_u01 = 5*k**2-4*k+1
    num_deri_u21_u01 = (res_list_1[0] - res_list[0])/h
    temp=compare(true_deri_u21_u01, num_deri_u21_u01, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u21_u02 = 2*k-4*k**2
    num_deri_u21_u02 = (res_list_2[0] - res_list[0])/h
    temp=compare(true_deri_u21_u02, num_deri_u21_u02, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u21_u03 = k**2
    num_deri_u21_u03 = (res_list_3[0] - res_list[0])/h
    temp=compare(true_deri_u21_u03, num_deri_u21_u03, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u0n_3 = k**2
    num_deri_u2n_1_u0n_3 = (res_list_n_3[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u0n_3, num_deri_u2n_1_u0n_3, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u0n_2 = 2*k-4*k**2
    num_deri_u2n_1_u0n_2 = (res_list_n_2[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u0n_2, num_deri_u2n_1_u0n_2, num)
    # print(true_deri_u2n_1_u0n_2, num_deri_u2n_1_u0n_2)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u0n_1 = 5*k**2-4*k+1
    num_deri_u2n_1_u0n_1 = (res_list_n_1[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u0n_1, num_deri_u2n_1_u0n_1, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u0n = k-2*k**2
    num_deri_u2n_1_u0n = (res_list_n[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u0n, num_deri_u2n_1_u0n, num)
    if temp: 
        num_of_detected_relation +=1 
    return num_of_detected_relation

if __name__ == '__main__':
    h = 0.0000001
    num = 6
    logger.add("D:\\Chrome_Download\\6_rq2.log")
    mu_list = [solver_mu1, solver_mu2, solver_mu3, solver_mu4, solver_mu5, solver_mu6, solver_mu7, solver_mu8, solver_mu9, solver_mu10]
    
    for i in range(10):
        detect_num = 0  ## the number of detected relations
        logger.info('****' + 'the ' +  str(i+1) +  ' mutant' + '****')
        # print('the ', str(i+1), 'mutant')
        try:
            num_of_detected_relation = main(h, mu_list[i], num)
            logger.info(str(i) + '@@@@@@@@@@@@@@@@@ num_of_detected_relation:' +   str(num_of_detected_relation))
        except Exception as e:
            print(e)
        # print('*********')
            logger.info(str(e))
        
        # main(h, solver)
    # logger.add("D:\\Chrome_Download\\6_compare_num_ana.log")
    # num_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # different keeping places
    # for num in num_list:
    #     logger.info('****' +   str(num)  + '****')
    #     num_of_detected_relation = main(h, solver, num)
    #     logger.info(str(num) + '@@@@@@@@@@@@@@@@@ num_of_detected_relation:' +   str(num_of_detected_relation))
