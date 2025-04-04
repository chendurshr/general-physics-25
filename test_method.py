# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 06:46:31 2025

@author: chend
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
# T_init = 0.0       # Initial room temperature (°C)
Q_max = 1500.0      # Heater power (W)
t_max = 3        # Simulation time in seconds (1 hour)
dt = .1          # Time step in seconds

# glass parameters
k_glass = .8       # W/m·K (thermal conductivity of glass)
d = 0.05          # m (thickness of glass wall)
A = 2.0             # m² (surface area of glass)
V = 20.0            # m³ (room volume)
c_air = 1005        # J/kg·K (specific heat capacity of air)
# k = (k_glass * A) / (m_air * c_air * d)  # Effective cooling coefficient (1/s)

# Define ODE function with heater logic


def air_mass(temperature):
    temperature_kelvin = temperature + 273.15
    air_density = 1.225 * (288.15 / temperature_kelvin)
    air_mass = air_density * V
    return air_mass


def heating_cooling_eqn(ext_temp, radiator_setting, int_temp, heater_on):
    T_low, T_high = radiator_setting
    if int_temp < T_low:
        heater_on = True  # Turn heater ON

    elif int_temp > T_high:
        heater_on = False  # Turn heater OFF

    m_air = air_mass(int_temp)
    k = (k_glass * A) / (m_air * c_air * d)
    Q_in = Q_max if heater_on else 0  # Heater power

    return (Q_in / (m_air * c_air)) - k * (int_temp - ext_temp), heater_on


def single_plot(time_arr, temp_arr, ext_temp, radiator_setting, date):
    # Plot results
    plt.plot(time_arr/3600, temp_arr, label="Temperature", color="b")
    plt.axhline(radiator_setting[0], linestyle="--",
                color="g", label=r"$T_{low}$ (Heater ON)")
    plt.axhline(radiator_setting[1], linestyle="--",
                color="r", label=r"$T_{high}$ (Heater OFF)")
    plt.axhline(ext_temp, linestyle=":", color="black",
                label="Outside Temperature")
    plt.xlabel("Time (hours)")
    plt.ylabel("Temperature (°C)")
    plt.title("Room Ambient Temperature on " + date)
    plt.legend()
    plt.ylim(ext_temp - 1, radiator_setting[1] + 1)
    plt.grid()
    plt.show()


def multi_plot(time_arr_list, temp_arr_list, ext_temp, radiator_setting_list, title, total_energy_list):
    # Plot results
    plt.figure(figsize=(12, 10))
    colors_list = [["maroon", "navy", "orange"], ["red", "cyan", "green"]]
    max_temp = np.max((radiator_setting_list))
    for colors, time_arr, temp_arr, radiator_setting, total_energy in zip(colors_list, time_arr_list, temp_arr_list, radiator_setting_list, total_energy_list):
        plt.plot(time_arr/3600, temp_arr,
                 label=f"Total energy:{total_energy:.0f} Wh", color=colors[2])
        plt.axhline(radiator_setting[0], linestyle="--",
                    color=colors[1], label=r"$T_{low}$")
        plt.axhline(radiator_setting[1], linestyle="--",
                    color=colors[0], label=r"$T_{high}$")
    plt.axhline(ext_temp, linestyle="--", color="black",
                label="Outside Temperature")
    plt.xlabel("Time (hours)")
    plt.ylabel("Temperature (°C)")
    plt.title(title)
    plt.legend()
    plt.ylim(ext_temp - 1, max_temp + 1)
    plt.grid()
    plt.show()
# Solve using Euler's method


def euler_method(ext_temp, rad_t_low, rad_t_high, duration, time_step, date, plot_cond=False):
    duration_seconds = duration * 3600
    time_array = np.arange(0, duration_seconds, time_step)
    temperature_array = np.zeros_like(time_array)
    temperature_array[0] = ext_temp
    heater_on = False  # Heater initially OFF
    total_energy = 0
    for i in range(1, len(time_array)):
        dTdt, heater_on = heating_cooling_eqn(
            ext_temp, (rad_t_low, rad_t_high), temperature_array[i-1], heater_on)
        if heater_on:
            total_energy += time_step * Q_max
        temperature_array[i] = temperature_array[i-1] + time_step * dTdt

    total_energy /= 3600 * duration
    if plot_cond:
        single_plot(time_array, temperature_array, ext_temp,
                    (rad_t_low, rad_t_high), date)
    return time_array, temperature_array, total_energy


#t_euler, T_euler, total_energy = euler_method(T_outside, radiator_setting, t_max-2, dt, plot_cond=True)

# print(total_energy)
# =============================================================================
#
# for temp in test_array:
#     t_euler, T_euler, total_energy = euler_method(
#         temp, 17, 20, t_max-2, dt, plot_cond=True)
#     t_euler_list.append(t_euler)
#     T_euler_list.append(T_euler)
#     total_energy_list.append(total_energy)
#
# =============================================================================
