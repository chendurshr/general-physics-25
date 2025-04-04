# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 04:17:21 2025

@author: chend
"""

import numerical_method
#import matplotlib.pyplot as plt
import numpy as np

dt = 0.05
t_euler_list = []
T_euler_list = []
total_energy_list = []
ext_temp = 10
t_max = 4
radiator_setting_list = np.array([[20, 22], [16, 18]])

for i in radiator_setting_list:
    low, high = i
    title = "Room temperature fluctuation for different radiator thresholds"

    t_euler, T_euler, total_energy = numerical_method.euler_method(
        ext_temp, low, high, t_max, dt, title, plot_cond=False)
    t_euler_list.append(t_euler)
    T_euler_list.append(T_euler)
    total_energy_list.append(total_energy)


numerical_method.multi_plot(
    t_euler_list, T_euler_list, ext_temp, radiator_setting_list, title, total_energy_list)
