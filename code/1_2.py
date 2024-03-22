import numpy as np
import math
import matplotlib.pyplot as plt
from loguru import logger

def solver_scientist(I, w, dt, T):
    """
    please simulate possible bug in this function and do not change the initial value
    Solve u'' + w**2*u = 0 for t in (0,T], u(0)=I and u'(0)=0,
    by a central finite difference method with time step dt.
    ##### DEFINE detect the bug successfully. And then we manually find the bug, and note by # **************** # .
    """
    res_list = []
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n] + u[n-1] - dt*2*w**2*u[n]  # **************** # change the sign of main iteration
        break
    res_list.append(u[2])
    return u

def solver(I, w, dt, T):
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n-1] - dt**2*w**2*u[n]
        break
    res_list.append(u[2])
    return res_list

def solver_mu1(I, w, dt, T): # sign
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n-1] + dt**2*w**2*u[n]   #### {+} dt**2*w**2*u[n]
        break
    res_list.append(u[2])
    return res_list

def solver_mu2(I, w, dt, T): # sign
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt): 
        u[n+1] = -2*u[n] - u[n-1] - dt**2*w**2*u[n] #### {-2}*u[n]
        break
    res_list.append(u[2])
    return res_list

def solver_mu3(I, w, dt, T): # sign
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n-1] / dt**2*w**2*u[n]  #### {/} dt**2*w**2*u[n]
        break
    res_list.append(u[2])
    return res_list

def solver_mu4(I, w, dt, T): # index
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n-1] - u[n-1] - dt**2*w**2*u[n]   #### 2*u[{n-1}]
        break
    res_list.append(u[2])
    return res_list

def solver_mu5(I, w, dt, T): # index
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n] - dt**2*w**2*u[n]  #### - u[{n}] 
        break
    res_list.append(u[2])
    return res_list

def solver_mu6(I, w, dt, T): # index
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n-1] - dt**2*w**2*u[n+1]   #### dt**2*w**2*u[{n+1}] 
        break
    res_list.append(u[2])
    return res_list

def solver_mu7(I, w, dt, T): # range
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(2, Nt):     #### range({2}, Nt):
        u[n+1] = 2*u[n] - u[n-1] - dt**2*w**2*u[n]
        break
    res_list.append(u[2])
    return res_list

def solver_mu8(I, w, dt, T): # instant
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 1.5*u[n] - u[n-1] - dt**2*w**2*u[n]   #### {1.5}*u[n]
        break
    res_list.append(u[2])
    return res_list

def solver_mu9(I, w, dt, T): # instant
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - 2*u[n-1] - dt**2*w**2*u[n]   #### - {2}*u[n-1]
        break
    res_list.append(u[2])
    return res_list

def solver_mu10(I, w, dt, T): # instant
    res_list = [] # we choose the first step and the second step
    dt = float(dt)
    Nt = int(round(T/dt))
    u = np.zeros(Nt+1)
    t = np.linspace(0, Nt*dt, Nt+1)  # Nt+1 elements in total
    u[0] = I
    u[1] = u[0] - 0.5*dt**2*w**2*u[0]
    res_list.append(u[1])
    for n in range(1, Nt):
        u[n+1] = 2*u[n] - u[n-1] - 0.5*dt**2*w**2*u[n]  #### - {0.5}*dt**2
        break
    res_list.append(u[2])
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

def main(h, func, num):
    num_of_detected_relation = 0
    # num = 1  # num=9 10 11, 12,13,14,15,solver都会出错。
    I = 1
    w = 2*math.pi
    # w = 0.0001
    dt = 0.05
    num_periods = 5
    P = 2*math.pi/w # one period
    T = P*num_periods
    res_list = func(I, w, dt, T)
    res_list_h = func(I+h, w, dt, T)
    # ana_value_u0 = 1 - 0.5*dt**2*w**2
    # num_value_u0 = (res_list_h[0] - res_list[0]) / h
    ana_value_u1 = 1 - 2 * dt**2*w**2 + 0.5*dt**4*w**4
    num_value_u1 = (res_list_h[1] - res_list[1]) / h
    # print(type(ana_value_u0), type(num_value_u0))
    # compare(ana_value_u0, num_value_u0, num)
    temp = compare(ana_value_u1, num_value_u1, num)
    if temp:
        num_of_detected_relation += 1
    return num_of_detected_relation

if __name__ == "__main__":
    h = 0.0000001
    num = 6
    logger.add("D:\\Chrome_Download\\1_rq2.log")
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
        
    
    
    # logger.add("D:\\Chrome_Download\\1_compare_num_ana.log")
    # num_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # different keeping places
    # for num in num_list:
    #     logger.info('****' +   str(num)  + '****')
    #     num_of_detected_relation = main(h, solver, num)
    #     logger.info(str(num) + '@@@@@@@@@@@@@@@@@ num_of_detected_relation:' +   str(num_of_detected_relation))
