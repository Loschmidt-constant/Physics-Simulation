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

dt = 0.05
t = 0
theta = a
phi = b

theta_list = []
theta2_list = []
t_list = []

while t < 10:
    t += dt
    theta += phi*dt
    phi += -k*math.sin(theta)*dt
    theta_list.append(theta)
    theta2 = a*math.cos(math.sqrt(k)*t)+(b/math.sqrt(k))*math.sin(math.sqrt(k)*t)
    theta2_list.append(theta2)
    t_list.append(t)

plt.plot(t_list, theta_list,marker="o",label="Numerical solution")
plt.plot(t_list, theta2_list,label="Analytical solution")
plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=8)
plt.xlabel('t'), plt.ylabel('theta')
plt.show()