# -*- coding: utf-8 -*-
"""
Script to run the rotated maps plotting routine.

Version:  v1.2 (2024-10-22)
Created on Wed Jun  5 09:23:30 2024

@author: jgazer (Jonathan Gasser, SwRI)
"""

import numpy as np
import matplotlib.pyplot as plt
from data_io import load_data
from map_rotated import plot_rotated_map

# --- Set Rotation Angles (in deg) ---
alfa = 35  # angle between old pole and new pole
phi = 45  # longitude where the new pole is
theta = 20 # end rotation towards positive azimuth

# --- Load Datafile ---
filename = 'example_data2.csv'
lons, lats, data = load_data(filename)

# lons = np.arange(0,301,31)
# lats = np.arange(-90,91,13)
# data = np.random.rand(12, 30)

# --- Create Mollweide Plot ---
fig, ax = plot_rotated_map(lons, lats, data, alf=alfa, phi=phi, th=theta,
                           title= 'Example Mollweide Map', grid_color= 'k')

## --- Save Plot as Postscript or PNG file:
# plt.savefig('figure.ps', format='ps')
plt.savefig('mollweide_example.png')

plt.show()
input('Press Enter to close the image and continue...')
