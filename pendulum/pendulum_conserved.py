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

dt = 0.01
t = 0
theta = a
phi = b
J = 0

theta_list = []
conserved_list = []
t_list = []

while t < 10:
    t += dt
    theta += phi*dt
    phi += -k*math.sin(theta)*dt
    J = 0.5*(phi**2) - k*math.cos(theta)
   #theta_list.append(theta)
    conserved_list.append(J)
    t_list.append(t)

print(conserved_list)

plt.plot(t_list, conserved_list, label="conserved quantity")
plt.xlabel('t'), plt.ylabel('J')
plt.show()