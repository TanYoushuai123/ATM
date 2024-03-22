import numpy as np
import copy
import math
import time
from loguru import logger

def solver_motivation(vis, nxc, nt, u, u00, u0n):
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) - (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) - (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) - (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) - (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list


def solver(vis, nxc, nt, u, u00, u0n):
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list

def solver_scientist(vis, nxc, nt, u, u00, u0n):
    """
    please simulate possible bug in this function and do not change the initial value
    ##### DEFINE detect the bug successfully. And then we manually find the bug, and note by # **************** # .
    """
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    for it in range(nt + 1):
        # print(u)
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])  
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(ip[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    # **************** # change the index of the main iteration
    return res_list

def solver_mu1(vis, nxc, nt, u, u00, u0n): # sign
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))

    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = -un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list   #### {-}un[i]

def solver_mu2(vis, nxc, nt, u, u00, u0n): # sign
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i] + (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list  #### {+} (un[i] * dt *

def solver_mu3(vis, nxc, nt, u, u00, u0n): # sign
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) / (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list     #### {/} (vis*dt*(un[int(ip[i])]-

def solver_mu4(vis, nxc, nt, u, u00, u0n): # index
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i-1] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list   #### = un[{i-1}]

def solver_mu5(vis, nxc, nt, u, u00, u0n): # index
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i+1] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list #### = un[{i+1}]

def solver_mu6(vis, nxc, nt, u, u00, u0n): # index
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(im[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list  #### + (vis*dt*(un[int({im[i]})]

def solver_mu7(vis, nxc, nt, u, u00, u0n): # range
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx-6):   ####  range({nx-6}):
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list

def solver_mu8(vis, nxc, nt, u, u00, u0n): # instant
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i] - (2*un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list  #### {2}*un[i] * dt *

def solver_mu9(vis, nxc, nt, u, u00, u0n): # instant
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = 0.5*un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list  #### {0.5}*un[i]

def solver_mu10(vis, nxc, nt, u, u00, u0n): # instant
    res_list = []
    nx = nxc + 1
    t = np.linspace(0, 0.5, nt+1)
    Lx = 2 * math.pi
    dx = Lx / nxc
    un = np.zeros(nx) 
    ip = np.zeros(nx) 
    im = np.zeros(nx) 
    dt = t[2] - t[1]
    # dx = x[1] - x[0]
    for i in range(nx):
        # boundary condition
        ip[i] = i + 1
        im[i] = i - 1
    ip[nx-1] = 1
    im[0] = nx - 2

    ## The first step
    un = copy.deepcopy(u) 
    for i in range(nx):
        if i == 0:
            u[i] = un[i] - (un[i] * dt * (un[i] - u00)/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+u00)/(dx*dx))
        elif i == nx-1:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(u0n-2*un[i]+un[int(im[i])])/(dx*dx))
        else:
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
  
    for it in range(nt + 1):
        if it == 1:
            res_list.append(u[0]) # corner values of the second step
            res_list.append(u[-1])     
        un = copy.deepcopy(u) 
        for i in range(nx):
            u[i] = un[i] - (un[i] * dt * (un[i] - un[int(im[i])])/dx) + (0.3*vis*dt*(un[int(ip[i])]-2*un[i]+un[int(im[i])])/(dx*dx))
    return res_list   #### {0.3}*vis*dt*(un[

def obtain_altered_initial_condition(h, index=None):
    nxc = 30 
    nt = 30
    t = np.linspace(0, 0.5, nt+1)
    dt = t[2] - t[1]
    Lx = 2 * math.pi
    dx = Lx / nxc
    ini = [math.sin(0 + n * dx)  for n in range(math.ceil(Lx / dx) + 1)]
    # ini = ini[1:-1]
    if index==None:
        return ini, dt, dx
    else:
        ini[index] += h
        return ini

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
        logger.info(str(ana_value) +'____' + str(num_value))
        return 0
    else:
        # print('There is a bug!!!')
        # print(ana_value_round, num_value_round)
        logger.info('There is a bug!!!')
        logger.info(str(ana_value) + '____' +str(num_value))
        return 1


def main(h, func, num):
    '''
    We take u21 u2n_1 as function f
    As for u21, the independent variables are: u00, u01, u02, u03, u0n_3, u0n_2, u0n_1
    With respect to u2n_1, the independent variables are: 
    '''
    num_of_detected_relation = 0
    vis = 0.3
    nxc = 30 # number of discretize points in x
    nt = 30   # number of discretize points in t
    ini, dt, dx = obtain_altered_initial_condition(h) ## unaltered initial condition
    u00, u01, u02, u03, u0n_3, u0n_2, u0n_1, u0n = ini[-2], ini[0], ini[1], ini[2], ini[-3], ini[-2], ini[-1], ini[1]
    res_list = func(vis,nxc,nt,ini, u00, u0n)
    ini_0, dt, dx = obtain_altered_initial_condition(h) ## altered 0 initial condition
    res_list_0 = func(vis,nxc,nt,ini_0, h+u00, u0n)
    ini_n, dt, dx = obtain_altered_initial_condition(h) ## altered n initial condition
    res_list_n = func(vis,nxc,nt,ini_n, u00, u0n+h)
    ini_1 = obtain_altered_initial_condition(h, 0) ## altered 1 initial condition
    res_list_1 = func(vis,nxc,nt,ini_1, u00, u0n)
    ini_2 = obtain_altered_initial_condition(h, 1) ## altered 2 initial condition
    res_list_2 = func(vis,nxc,nt,ini_2, u00, u0n)
    ini_3 = obtain_altered_initial_condition(h, 2) ## altered 3 initial condition
    res_list_3 = func(vis,nxc,nt,ini_3, u00, u0n)
    ini_n_3 = obtain_altered_initial_condition(h, -3) ## altered n-3 initial condition
    res_list_n_3 = func(vis,nxc,nt,ini_n_3, u00, u0n)
    ini_n_2 = obtain_altered_initial_condition(h, -2) ## altered n-2 initial condition
    res_list_n_2 = func(vis,nxc,nt,ini_n_2, u00, u0n)
    ini_n_1 = obtain_altered_initial_condition(h, -1) ## altered n-1 initial condition
    res_list_n_1 = func(vis,nxc,nt,ini_n_1, u00, u0n)

    ## for u21
    x = dt/dx
    y = vis*dt/(dx*dx)
    true_deri_u21_u00 = u0n_3*x*y**2-2*u0n_2*x*y**2+u0n_1*x*y**2-2*u02*x*y**2+4*u01*x*y**2-2*u00*x*y**2-2*y**2+u0n_2*u0n_3*x**2*y+u01*u0n_3*x**2*y-u0n_2**2*x**2*y-2*u01*u0n_2*x**2*y+u01*u0n_1*x**2*y-2*u01*u02*x**2*y+6*u01**2*x**2*y-4*u00*u01*x**2*y+u0n_2*x*y-4*u01*x*y+y+u01*u0n_2*u0n_3*x**3-u01*u0n_2**2*x**3+2*u01**3*x**3-2*u00*u01**2*x**3+u01*u0n_2*x**2-2*u01**2*x**2+u01*x
    num_deri_u21_u00 = (res_list_0[0] - res_list[0])/h
    temp=compare(true_deri_u21_u00, num_deri_u21_u00, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u21_u01 = -(2*u0n_3*x*y**2)+4*u0n_2*x*y**2-2*u0n_1*x*y**2+4*u02*x*y**2-8*u01*x*y**2+4*u00*x*y**2+5*y**2-2*u0n_2*u0n_3*x**2*y-2*u01*u0n_3*x**2*y+u00*u0n_3*x**2*y+2*u0n_2**2*x**2*y+4*u01*u0n_2*x**2*y-2*u00*u0n_2*x**2*y-2*u01*u0n_1*x**2*y+u00*u0n_1*x**2*y+4*u01*u02*x**2*y-2*u00*u02*x**2*y-12*u01**2*x**2*y+12*u00*u01*x**2*y-2*u00**2*x**2*y+u0n_3*x*y-4*u0n_2*x*y+u0n_1*x*y-u02*x*y+12*u01*x*y-4*u00*x*y-4*y-2*u01*u0n_2*u0n_3*x**3+u00*u0n_2*u0n_3*x**3+2*u01*u0n_2**2*x**3-u00*u0n_2**2*x**3-4*u01**3*x**3+6*u00*u01**2*x**3-2*u00**2*u01*x**3+u0n_2*u0n_3*x**2-u0n_2**2*x**2-2*u01*u0n_2*x**2+u00*u0n_2*x**2+6*u01**2*x**2-4*u00*u01*x**2+u0n_2*x-4*u01*x+u00*x+1
    num_deri_u21_u01 = (res_list_1[0] - res_list[0])/h
    temp=compare(true_deri_u21_u01, num_deri_u21_u01, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u21_u02 = u0n_3*x*y**2-2*u0n_2*x*y**2+u0n_1*x*y**2-2*u02*x*y**2+4*u01*x*y**2-2*u00*x*y**2-4*y**2+u0n_2*u0n_3*x**2*y-u0n_2**2*x**2*y+2*u01**2*x**2*y-2*u00*u01*x**2*y+u0n_2*x*y-2*u02*x*y-u01*x*y+2*y
    num_deri_u21_u02 = (res_list_2[0] - res_list[0])/h
    temp=compare(true_deri_u21_u02, num_deri_u21_u02, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u21_u03 = y**2
    num_deri_u21_u03 = (res_list_3[0] - res_list[0])/h
    temp=compare(true_deri_u21_u03, num_deri_u21_u03, num)
    # print(true_deri_u21_u03, num_deri_u21_u03)
    if temp: 
        num_of_detected_relation +=1 

    true_deri_u21_u0n_3 = u02*x*y**2-2*u01*x*y**2+u00*x*y**2+y**2+u02*u0n_2*x**2*y-2*u01*u0n_2*x**2*y+u00*u0n_2*x**2*y-u01**2*x**2*y+u00*u01*x**2*y+u0n_2*x*y+u01*x*y-u01**2*u0n_2*(x**3)+u00*u01*u0n_2*(x**3)+u01*u0n_2*x**2
    num_deri_u21_u0n_3 = (res_list_n_3[0] - res_list[0])/h
    temp=compare(true_deri_u21_u0n_3, num_deri_u21_u0n_3, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u21_u0n_2 = -(2*u02*x*y**2)+4*u01*x*y**2-2*u00*x*y**2-2*y**2+u02*u0n_3*x**2*y-2*u01*u0n_3*x**2*y+u00*u0n_3*x**2*y-2*u02*u0n_2*x**2*y+4*u01*u0n_2*x**2*y-2*u00*u0n_2*x**2*y+2*u01**2*x**2*y-2*u00*u01*x**2*y+u0n_3*x*y-2*u0n_2*x*y+u02*x*y-4*u01*x*y+u00*x*y+y-u01**2*u0n_3*(x**3)+u00*u01*u0n_3*(x**3)+2*u01**2*u0n_2*(x**3)-2*u00*u01*u0n_2*(x**3)+u01*u0n_3*x**2-2*u01*u0n_2*x**2-u01**2*x**2+u00*u01*x**2+u01*x
    num_deri_u21_u0n_2 = (res_list_n_2[0] - res_list[0])/h
    temp=compare(true_deri_u21_u0n_2, num_deri_u21_u0n_2, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u21_u0n_1 = u02*x*y**2-2*u01*x*y**2+u00*x*y**2+y**2-u01**2*x**2*y+u00*u01*x**2*y+u01*x*y
    num_deri_u21_u0n_1 = (res_list_n_1[0] - res_list[0])/h
    temp=compare(true_deri_u21_u0n_1, num_deri_u21_u0n_1, num)
    if temp: 
        num_of_detected_relation +=1 
    ## for u2n_1
    true_deri_u2n_1_u01 = y**2+u02*x*y
    num_deri_u2n_1_u01 = (res_list_1[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u01, num_deri_u2n_1_u01, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u02 = -(2*y**2)-2*u02*x*y+u01*x*y+y
    num_deri_u2n_1_u02 = (res_list_2[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u02, num_deri_u2n_1_u02, num)
    if temp: 
        num_of_detected_relation +=1 

    true_deri_u2n_1_u03 = y**2
    num_deri_u2n_1_u03 = (res_list_3[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u03, num_deri_u2n_1_u03, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u0n_3 = u0n_2*x*y**2-2*u0n_1*x*y**2+u0n*x*y**2+y**2+u0n_2**2*x**2*y-u0n_1*u0n_2*x**2*y+u0n*u0n_2*x**2*y-u0n_1**2*x**2*y+u0n_2*x*y+u0n_1*x*y+u0n_1*u0n_2**2*x**3-u0n_1**2*u0n_2*x**3+u0n_1*u0n_2*x**2
    num_deri_u2n_1_u0n_3 = (res_list_n_3[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u0n_3, num_deri_u2n_1_u0n_3, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u0n_2 = u0n_3*x*y**2-6*u0n_2*x*y**2+9*u0n_1*x*y**2-4*u0n*x*y**2-4*y**2+2*u0n_2*u0n_3*x**2*y-u0n_1*u0n_3*x**2*y+u0n*u0n_3*x**2*y-3*u0n_2**2*x**2*y-4*u0n_1*u0n_2*x**2*y-2*u0n*u0n_2*x**2*y+9*u0n_1**2*x**2*y-2*u0n*u0n_1*x**2*y+u0n_3*x*y-8*u0n_1*x*y+u0n*x*y+2*y+2*u0n_1*u0n_2*u0n_3*x**3-u0n_1**2*u0n_3*x**3-3*u0n_1*u0n_2**2*x**3+2*u0n_1**3*x**3+u0n_1*u0n_3*x**2-3*u0n_1**2*x**2+2*u0n_1*x
    num_deri_u2n_1_u0n_2 = (res_list_n_2[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u0n_2, num_deri_u2n_1_u0n_2, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u0n_1 = -(2*u0n_3*x*y**2)+9*u0n_2*x*y**2-12*u0n_1*x*y**2+5*u0n*x*y**2+5*y**2-u0n_2*u0n_3*x**2*y-2*u0n_1*u0n_3*x**2*y-2*u0n_2**2*x**2*y+18*u0n_1*u0n_2*x**2*y-2*u0n*u0n_2*x**2*y-15*u0n_1**2*x**2*y+4*u0n*u0n_1*x**2*y+u0n_3*x*y-8*u0n_2*x*y+14*u0n_1*x*y-2*u0n*x*y-4*y+u0n_2**2*u0n_3*x**3-2*u0n_1*u0n_2*u0n_3*x**3-u0n_2**3*x**3+6*u0n_1**2*u0n_2*x**3-4*u0n_1**3*x**3+u0n_2*u0n_3*x**2-6*u0n_1*u0n_2*x**2+6*u0n_1**2*x**2+2*u0n_2*x-4*u0n_1*x+1
    num_deri_u2n_1_u0n_1 = (res_list_n_1[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u0n_1, num_deri_u2n_1_u0n_1, num)
    if temp: 
        num_of_detected_relation +=1 
    true_deri_u2n_1_u0n = u0n_3*x*y**2-4*u0n_2*x*y**2+5*u0n_1*x*y**2-2*u0n*x*y**2-2*y**2+u0n_2*u0n_3*x**2*y-u0n_2**2*x**2*y-2*u0n_1*u0n_2*x**2*y+2*u0n_1**2*x**2*y+u0n_2*x*y-2*u0n_1*x*y+y
    num_deri_u2n_1_u0n = (res_list_n[1] - res_list[1])/h
    temp=compare(true_deri_u2n_1_u0n, num_deri_u2n_1_u0n, num)
    if temp: 
        num_of_detected_relation +=1 
    return num_of_detected_relation

if __name__ == "__main__":
    h = 0.0000001
    num = 6
    logger.add("D:\\Chrome_Download\\7_rq2.log")
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
        
    #     # main(h, solver)
    # logger.add("D:\\Chrome_Download\\7_compare_num_ana.log")
    # num_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # different keeping places
    # for num in num_list:
    #     logger.info('****' +   str(num)  + '****')
    #     num_of_detected_relation = main(h, solver, num)
    #     logger.info(str(num) + '@@@@@@@@@@@@@@@@@ num_of_detected_relation:' +   str(num_of_detected_relation))
