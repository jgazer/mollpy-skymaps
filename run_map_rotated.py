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
alf = 35  # angle between old pole and new pole
phi = 45  # longitude where the new pole is
th = 20 # end rotation towards positive azimuth

# --- Load Datafile ---
lons, lats, data = load_data('example_data2.csv')

# --- Create Mollweide Plot ---
fig, ax = plot_rotated_map(lons, lats, data, alf=alf, phi=phi, th=th,
                           title= 'Example Mollweide Map', grid_color= 'k')

## --- Save Plot as Postscript or PNG file:
# plt.savefig('figure.ps', format='ps')
# plt.savefig('mollweide_example.png')

plt.show()
