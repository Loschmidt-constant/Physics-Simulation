# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 16:50:08 2019

@author: A.Nishimura
"""

from matplotlib import pyplot as plt
import math

def draw(dt):
    #初期値
    k = 10
    a = 0.1
    b = 0
    
    #初期条件
    t = 0
    theta = a
    phi = b
    
    theta_list = []
    phi_list = []
    t_list = []
    
    while t < 5:
        t += dt
        theta += phi*dt
        phi += -k*math.sin(theta)*dt
        theta_list.append(theta)
        t_list.append(t)
    plt.plot(t_list, theta_list)
    plt.xlabel('t')
    plt.ylabel('theta')

if __name__ == "__main__": 
    dt_list=[0.5,0.25,0.10,0.01]
    for dt in dt_list:
        draw(dt)       
    plt.legend(['dt = 0.5','dt = 0.25','dt = 0.10','dt = 0.01'])
    plt.show()