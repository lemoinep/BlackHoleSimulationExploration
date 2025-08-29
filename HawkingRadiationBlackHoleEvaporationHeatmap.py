
# Author(s): Dr. Patrick Lemoine

# Hawking EvapHeat: Visualizing Black Hole Radiation and Evaporation
# Visual information : This heatmap represents how close a given temperature is to the expected Hawking temperature for black holes of different masses. 
# Dark (low difference) regions indicate agreement  with Hawking radiation theory, visually emphasizing key physics 
# such as the final explosive evaporation phase at low masses and high temperatures.

import numpy as np
import matplotlib.pyplot as plt

# Constants
hbar = 1.0545718e-34  # Reduced Planck constant (J.s)
c = 3e8               # Speed of light (m/s)
G = 6.67430e-11       # Gravitational constant (m^3/kg/s^2)
k_B = 1.380649e-23    # Boltzmann constant (J/K)
M_solar = 1.98847e30  # Solar mass (kg)
CMB_temp = 2.7        # Cosmic Microwave Background temperature (K)

def hawking_temperature(M):
    return (hbar * c**3) / (8 * np.pi * G * M * k_B)

def evaporation_time(M):
    return (5120 * np.pi * G**2 * M**3) / (hbar * c**4)

# Create logarithmic mass and temperature grids
mass_vals = np.logspace(11, 15, 400)  # Mass range (kg)
temp_vals = np.logspace(1, 20, 400)   # Temperature range (K)

M_grid, T_grid = np.meshgrid(mass_vals, temp_vals)

# Compute Hawking temperature for each mass grid point
T_hawking = hawking_temperature(M_grid)

# Compute the absolute difference between grid temperatures and Hawking temp
difference = np.abs(T_grid - T_hawking)

# Plot heatmap of log10 of absolute temperature difference
plt.figure(figsize=(12, 8))
heatmap = plt.contourf(np.log10(mass_vals / M_solar), np.log10(temp_vals),
                       np.log10(difference), levels=100, cmap='inferno')

plt.colorbar(heatmap, label='log10 Absolute Temperature Difference')
plt.xlabel('log10 Black Hole Mass (Solar Masses)')
plt.ylabel('log10 Temperature (Kelvin)')

plt.title('Heatmap Visualization of Hawking Temperature Fit\n'
          'Low difference values highlight expected Hawking temperature for given mass',
          color='red', fontsize=14)

# Annotations for interpretation
plt.text(-12, 5, 'Lower mass → Higher temperature → Faster evaporation,\nFinal Explosive Phase region',
         color='red', fontsize=12)
plt.text(-10, -1, 'Temperature below CMB (~2.7 K)\nMinimal evaporation region', color='blue', fontsize=12)

plt.show()
