import numpy as np
import math
from loguru import logger

def solver(I, V, f, c, L, dt, C, T):
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i] + 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_scientist(I, V, f, c, L, dt, C, T):
    
    """
    please simulate possible bug in this function and do not change the initial value
    Solve u_tt=c^2*u_xx + f on (0,L)x(0,T].
    ##### DEFINE detect the bug successfully. And then we manually find the bug, and note by # **************** # .
    """
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):                                                 # change the index of main iteration              
            u[i]=- u_nm1[i] + 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i]) + dt**2*f(x[i], t[n])# **************** # .

        # Insert boundary conditions
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu1(I, V, f, c, L, dt, C, T): # sign
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=+ u_nm1[i] + 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions  #### {+} u_nm1[i] 
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu2(I, V, f, c, L, dt, C, T): # sign
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i] - 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions     #### {-} 2*u_n[i]
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu3(I, V, f, c, L, dt, C, T): # sign
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] * dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0         ####  {*} dt*V(x[i])
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i] + 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu4(I, V, f, c, L, dt, C, T): # index
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i-1] + 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions   #### - u_nm1[{i-1}] 
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu5(I, V, f, c, L, dt, C, T): # index
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i] + 2*u_n[i-1] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions   #### 2*u_n[{i-1}]
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu6(I, V, f, c, L, dt, C, T): # index
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i+1] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0    #### = u_n[{i+1}]
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i] + 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu7(I, V, f, c, L, dt, C, T): # range
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx-1):   #### range(1, {Nx-1})
            u[i]=- u_nm1[i] + 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu8(I, V, f, c, L, dt, C, T): # instant
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i] + 1*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions   #### + {1}*u_n[i]
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu9(I, V, f, c, L, dt, C, T): # instant
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 0.5*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i] + 2*u_n[i] + 2*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions  #### {2}*C2*(u_n[i-1]
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

def solver_mu10(I, V, f, c, L, dt, C, T): # instant
    
    res_list = [] # the second and the third time step
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1) # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1) # Mesh points in space

    # Make sure dx and dt are compatible with x and t
    dx = x[1] - x[0]
    dt = t[1] - t[0]
    if f is None or f == 0 :
        f = lambda x, t: 0
    if V is None or V == 0:
        V = lambda x: 0
    C2 = C**2 * 0.6 # Help variable in the scheme
    u = np.zeros(Nx+1) # Solution array at new time level
    u_n = I # Solution at 1 time level back
    u_nm1 = np.zeros(Nx+1) # Solution at 2 time levels back
    import time; t0 = time.process_time() # Measure CPU time
    # Load initial condition into u_n

    # Special formula for first time step 
    n = 0
    for i in range(1, Nx):
        u[i] = u_n[i] + dt*V(x[i]) + 1*C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + 0.5*dt**2*f(x[i], t[n])
    u[0] = 0; u[Nx]=0    ####  {1}*C2*(u_n[i-1]
    # Switch variables before next step
    u_nm1[:] = u_n; u_n[:] = u
    for n in range(1, Nt):
        if n == 2:
            res_list.append(u[1])
            res_list.append(u[-2])
        if n == 3:
            res_list.append(u[1])
            res_list.append(u[-2])
        # Update all inner points at time t[n+1]
        for i in range(1, Nx):
            u[i]=- u_nm1[i] + 2*u_n[i] + C2*(u_n[i-1] - 2*u_n[i] + u_n[i+1]) + dt**2*f(x[i], t[n])
        # Insert boundary conditions
        u[0] = 0; u[Nx]=0
        # Switch variables before next step
        u_nm1[:] = u_n; u_n[:] = u
    cpu_time = time.process_time() - t0
    return res_list

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
        return 0
    else:
        # print('There is a bug!!!')
        # print(ana_value_round, num_value_round)
        logger.info('There is a bug!!!')
        logger.info(str(ana_value) + '___' + str(num_value))
        return 1

def obtain_altered_initial_condition(h, index=None):
    L = 0.75 # sapce
    x0 = 0.8*L
    a = 0.005
    freq = 440
    wavelength = 2*L
    c = freq*wavelength
    omega = 2*math.pi*freq
    num_periods = 1
    T = 2*math.pi/omega*num_periods # total time; it is small 
    C = 1
    dt = L/50./c

    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x_mesh = np.linspace(0, L, Nx+1)
    I = np.zeros(Nx+1)
    for i in range(0,Nx+1):
        I[i] = math.sin(x_mesh[i]*math.pi*10)
    I[0] = 0
    I[Nx] = 0
    if index == None:
        return I
    else:
        I[index] += h
        return I

def main(h, func, num):
    num_of_detected_relation = 0
    L = 0.75 # sapce
    x0 = 0.8*L
    a = 0.005
    freq = 440
    wavelength = 2*L
    c = freq*wavelength
    omega = 2*math.pi*freq
    num_periods = 1
    T = 2*math.pi/omega*num_periods # total time; it is small 
    C = 1
    dt = L/50./c

    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x_mesh = np.linspace(0, L, Nx+1)
    g = C**2 * 0.6
    n = Nx
    ini = obtain_altered_initial_condition(h)
    ini_res = func(ini, 0, 0, c, L, dt, C, T) ## unaltered

    ini_0 = obtain_altered_initial_condition(h, 0)
    ini_res_0 = func(ini_0, 0, 0, c, L, dt, C, T)

    ini_1 = obtain_altered_initial_condition(h, 1)
    ini_res_1 = func(ini_1, 0, 0, c, L, dt, C, T)

    ini_2 = obtain_altered_initial_condition(h, 2)
    ini_res_2 = func(ini_2, 0, 0, c, L, dt, C, T)

    ini_3 = obtain_altered_initial_condition(h, 3)
    ini_res_3 = func(ini_3, 0, 0, c, L, dt, C, T)

    ini_n = obtain_altered_initial_condition(h, n)
    ini_res_n = func(ini_n, 0, 0, c, L, dt, C, T)

    ini_n_1 = obtain_altered_initial_condition(h, n-1)
    ini_res_n_1 = func(ini_n_1, 0, 0, c, L, dt, C, T)

    ini_n_2 = obtain_altered_initial_condition(h, n-2)
    ini_res_n_2 = func(ini_n_2, 0, 0, c, L, dt, C, T)

    ini_n_3 = obtain_altered_initial_condition(h, n-3)
    ini_res_n_3 = func(ini_n_3, 0, 0, c, L, dt, C, T)
    
    ana_u21_u00 = 1.0*g-1.0*g**2
    num_u21_u00 = (ini_res_0[0] - ini_res[0])/h
    temp = compare(ana_u21_u00, num_u21_u00, num)
    if temp: 
        num_of_detected_relation +=1 

    ana_u21_u01 = 2.5*g**2-4.0*g+1
    num_u21_u01 = (ini_res_1[0] - ini_res[0])/h
    temp = compare(ana_u21_u01, num_u21_u01, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u21_u02 = 2.0*g-2.0*g**2
    num_u21_u02 = (ini_res_2[0] - ini_res[0])/h
    temp = compare(ana_u21_u02, num_u21_u02, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u21_u03 = 0.5*g**2
    num_u21_u03 = (ini_res_3[0] - ini_res[0])/h
    temp = compare(ana_u21_u03, num_u21_u03, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u2n_1_u0n = 1.0*g-1.0*g**2
    num_u2n_1_u0n = (ini_res_n[1] - ini_res[1])/h
    temp = compare(ana_u2n_1_u0n, num_u2n_1_u0n, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u2n_1_u0n_1 = 2.5*g**2-4.0*g+1
    num_u2n_1_u0n_1 = (ini_res_n_1[1] - ini_res[1])/h
    temp = compare(ana_u2n_1_u0n_1, num_u2n_1_u0n_1, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u2n_1_u0n_2 = 2.0*g-2.0*g**2
    num_u2n_1_u0n_2 = (ini_res_n_2[1] - ini_res[1])/h
    temp = compare(ana_u2n_1_u0n_2, num_u2n_1_u0n_2, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u2n_1_u0n_3 = 0.5*g**2
    num_u2n_1_u0n_3 = (ini_res_n_3[1] - ini_res[1])/h
    temp = compare(ana_u2n_1_u0n_3, num_u2n_1_u0n_3, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u31_u00 = 2.5*g**3-4*g**2+1.5*g
    num_u31_u00 = (ini_res_0[2] - ini_res[2])/h
    temp = compare(ana_u31_u00, num_u31_u00, num)
    if temp: 
        num_of_detected_relation +=1 
    ana_u31_u01 = -(7.0*g**3)+15.0*g**2-9.0*g+1
    num_u31_u01 = (ini_res_1[2] - ini_res[2])/h
    temp = compare(ana_u31_u01, num_u31_u01, num)
    if temp: 
        num_of_detected_relation +=1 
    return num_of_detected_relation  ### 

if __name__ == "__main__":
    h = 0.0000001
    num = 6
    logger.add("D:\\Chrome_Download\\2_rq2.log")
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
            logger.info(e)
        
    # logger.add("D:\\Chrome_Download\\2_compare_num_ana.log")
    # num_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # different keeping places
    # for num in num_list:
    #     logger.info('****' +   str(num)  + '****')
    #     num_of_detected_relation = main(h, solver, num)
    #     logger.info(str(num) + '@@@@@@@@@@@@@@@@@ num_of_detected_relation:' +   str(num_of_detected_relation))
