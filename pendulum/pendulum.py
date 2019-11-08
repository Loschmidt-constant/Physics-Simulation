# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 14:03:08 2019

@author: A.Nishimura
"""

from matplotlib import pyplot as plt
import math

k = 10
a = 0.1
b = 0

dt = 0.001
t = 0
theta = a
phi = b

theta_list = []
phi_list = []
t_list = []

while t < 10:
    t += dt
    theta += phi*dt
    phi += -k*math.sin(theta)*dt
    theta_list.append(theta)
    phi_list.append(phi)
    t_list.append(t)

plt.legend([])
plt.plot(t_list, theta_list)
plt.xlabel('t'), plt.ylabel('theta')
plt.show()
