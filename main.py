# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 12:43:46 2025

@author: chend

1. Talk about hysterisis
"""

import numpy as np
import scipy.integrate as spi
import os
import scipy.special as sp
import numerical_method
import matplotlib.pyplot as plt

num_rooms = 7000


def std_calc(lower, mid, upper):
    lower_std = (mid - lower) / 4
    upper_std = (upper - mid) / 4
    mean_std = np.mean((lower_std, upper_std))
    duration = upper - lower - 3 * mean_std
    return mean_std, duration


std_calc = np.vectorize(std_calc)

# Define Gaussian function


def gaussian(max_val, mean, std, x):
    return max_val * np.exp(-((x - mean) ** 2) / (2 * std ** 2))


def gaussian_area_analytical(peak_val, sigma):
    return peak_val * sigma * np.sqrt(2 * np.pi)


def generalized_gaussian(max_val, mean, std, x, p):
    return max_val * np.exp(-((x - mean) / std)**p)


def generalized_gaussian_integral(A, sigma, p):
    return A * sigma * (2 / p) * sp.gamma(1 / p)


generalized_gaussian_integral = np.vectorize(generalized_gaussian_integral)
gaussian_area_analytical = np.vectorize(gaussian_area_analytical)


# Numerical integration over a wide range
#result, error = spi.quad(gaussian, -np.inf, np.inf, args=(A, mu, sigma))

# Analytical solution
#analytical_result = A * sigma * np.sqrt(2 * np.pi)
irradiance_data_path = r"C:\Users\chend\Desktop\Projects\Yr3\Sem2\genphys\irradiance data"
irr_data_path_list = os.listdir(irradiance_data_path)
irr_data_path_list = [os.path.join(irradiance_data_path, path)
                      for path in irr_data_path_list]


temp_data_path = r"C:\Users\chend\Desktop\Projects\Yr3\Sem2\genphys\temperature data"
temp_data_path_list = os.listdir(temp_data_path)
temp_data_path_list = [os.path.join(temp_data_path, path)
                       for path in temp_data_path_list]
mean_val_list = []
peak_at_mean_list = []
lower_val_list = []
upper_val_list = []
std_list = []
total_irradiance_normal_gauss_list = []
total_irradiance_general_gauss_list = []
temp_list = []
duration_list = []

for data_path in irr_data_path_list:
    data = np.genfromtxt(data_path, delimiter=",")
    mean_val_list.append(data[:, 0])
    peak_at_mean_list.append(data[:, 1])
    lower_val_list.append(data[:, 2])
    upper_val_list.append(data[:, 3])

    std_val, duration_val = std_calc(data[:, 2], data[:, 0], data[:, 3])
    std_list.append(std_val)
    duration_list.append(duration_val)

    total_irradiance_normal_gauss_list.append(gaussian_area_analytical(
        data[:, 1], std_val))
    total_irradiance_general_gauss_list.append(
        generalized_gaussian_integral(data[:, 1], std_val, 0.475))


for data_path in temp_data_path_list:
    data = np.genfromtxt(data_path, delimiter=",")
    temp_list.append(data)

print("data acq flag")

t_euler_list = []
T_euler_list = []
total_energy_list = []
total_energy_month = []
month_list = ["September", "October", "November",
              "December", "January", "February"]
dt = 0.1

for month in range(len(month_list)):
    print(month)
    total_energy_month = []
    for i in range(0, len(temp_list[month])):
        date = month_list[month] + " " + str(i+1)
        t_max = duration_list[month][i]
        t_euler, T_euler, total_energy = numerical_method.improved_euler_method(
            temp_list[month][i], 20, 24, t_max, dt, date, plot_cond=False)
        t_euler_list.append(t_euler)
        T_euler_list.append(T_euler)
        total_energy_month.append(total_energy)
    total_energy_list.append(total_energy_month)
print("total energy calc flag")

fig2, axs2 = plt.subplots(2, 3, figsize=(18, 10))
axs2 = axs2.flatten()
min_sol_irrad_list = []
for month in range(6):
    solar_watt = total_irradiance_general_gauss_list[month] / \
        duration_list[month]
    ax = axs2[month]
    sol_irrad_min = np.argmin(solar_watt)
    min_sol_irrad_list.append(solar_watt[sol_irrad_min])
    ax.bar(
        range(1, len(temp_list[month]) + 1),
        solar_watt * 1000,
        width=1,
        edgecolor="k"
    )
    ax.bar(sol_irrad_min + 1, solar_watt[sol_irrad_min]*1000, width=1,
           edgecolor="k", color="midnightblue", label="Minimum Daily Solar power")
    ax.set_title(f"{month_list[month]}", fontsize=14)
    ax.set_xlabel("Day of the month")
    ax.set_ylabel("Solar Power (W/m²)")
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    ax.legend()
# fig2.suptitle("Average Solar Power per Day (First 6 Months)",
    #     fontsize=18, y=1.02)
plt.tight_layout()

plt.show()

fig, axs = plt.subplots(2, 3, figsize=(18, 10))
axs = axs.flatten()  # To index subplots in a 1D loop
max_num_panels_list = []
for month in range(len(month_list)):
    solar_watt = total_irradiance_general_gauss_list[month] / \
        duration_list[month]
    solar_output = solar_watt * \
        numerical_method.solar_panel_efficiency(temp_list[month])
    num_panels = total_energy_list[month]/(solar_output*1000)
    num_panels *= num_rooms
    day_max = np.argmax(num_panels)
    max_num_panels_list.append(num_panels[day_max])
    ax = axs[month]
    ax.bar(range(1, len(temp_list[month]) + 1), num_panels,
           width=1, edgecolor="k", label='Number of Panels')
    ax.bar(day_max + 1, num_panels[day_max], width=1,
           edgecolor="k", color="tomato", label=fr"Max Panel Requirement($m^2$) ({num_panels[day_max]:.1f})")

    ax.set_title(f"{month_list[month]}", fontsize=14)
    ax.set_xlabel("Day of the month")
    ax.set_ylabel("Number of panels")
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    ax.legend()

#plt.suptitle("Number of Solar Panels Needed Each Day", fontsize=18, y=1.02)
plt.tight_layout()
plt.show()
