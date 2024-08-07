##To plot probability density function 

import sys
import os
import argparse

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
librarydir = os.path.dirname(parentdir) + '/0_libraries'
sys.path.append(librarydir)

import llib_customfunction
import llib_figures as figs
import llib_cellcycles as cc
import matplotlib.pyplot as plt
import numpy as np
import statistics as stat
import seaborn as sns
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from math import nan
from math import sqrt

#plt.style.use('nishitstyle')    

###########FUNCTIONS################

def print_to_file(filename, column_vector, no_of_columns):
    column_length_vec = [len(column_vector[i]) for i in range(no_of_columns)]

    if len(set(column_length_vec)) == 1:
        column_length = column_length_vec[0]    
        file = open(filename, "w+")
        for j in range(column_length):
            for i in range(no_of_columns):
                file.write('%s ' % column_vector[i][j])  
            file.write('\n')
        file.close()    
    else :
        print('All columns must have the same length.')

def filter_data(scatter, filter_percent):
    x, y = scatter
    # Compute percentiles
    x_10, x_90 = np.percentile(x, [filter_percent, 100-filter_percent])
    y_10, y_90 = np.percentile(y, [filter_percent, 100-filter_percent])
    
    # Filter data
    filtered_data = [(xi, yi) for xi, yi in zip(x, y) if x_10 < xi < x_90 and y_10 < yi < y_90]
    
    # Unzip the filtered data to get new x and y lists
    x_new, y_new = zip(*filtered_data)
    
    return [x_new, y_new]  

def find_nearest_index(times, target_time):
    """Finds the index of the nearest time in the list to the target time."""
    times = np.array(times)
    index = np.abs(times - target_time).argmin()
    return index
def cell_cycle_time_series(filename):
    with open(filename, 'r+') as file: 
        line_index = 0
        cells = []
        birth_times = []
        cell_time = []
        cell_obs = []    
        for line in file: 
            data_point = line.split()
            if line_index == 0:
                prev_cell_id = data_point[0]
                cell_time.append(float(data_point[1]))
                cell_obs.append(float(data_point[2]))
                birth_times.append(float(data_point[1]))  # Add birth time
            else:
                cell_id = data_point[0]
                if  cell_id == prev_cell_id:
                    cell_time.append(float(data_point[1]))
                    cell_obs.append(float(data_point[2]))
                else:
                    cells.append([cell_time, cell_obs])
                    cell_time = []
                    cell_obs = []    
                    cell_time.append(float(data_point[1]))
                    cell_obs.append(float(data_point[2]))
                    birth_times.append(float(data_point[1]))  # Add birth time
                    prev_cell_id = cell_id
            line_index += 1
    cells.append([cell_time, cell_obs])            
    return cells, birth_times

def density_at_birth_times(density_data, birth_times):
    densities_at_birth = []
    for i in range(len(density_data)):
        cell_time = density_data[i][0]
        cell_density = density_data[i][1]
        birth_time = birth_times[i]
        
        # Find the density at birth time
        # Assuming birth time matches the first time point for each cell
        density_at_birth = cell_density[0] if cell_time[0] == birth_time else None
        densities_at_birth.append(density_at_birth)
    
    return densities_at_birth

# Function to get density data during mitosis
def cells_obs_in_mitosis(cells_obs_exp, mitotic_times):
    density_during_mitosis = []
    
    for i in range(len(cells_obs_exp)):
        density_for_cell = []
        
        for k in range(len(cells_obs_exp[i][0])):
            if isinstance(cells_obs_exp[i][0][k], list):
                if any(t > mitotic_times[i] for t in cells_obs_exp[i][0][k]):
                    density_for_cell.append(cells_obs_exp[i][1][k])
            else:
                if cells_obs_exp[i][0][k] > mitotic_times[i]:
                    density_for_cell.append(cells_obs_exp[i][1][k])
        
        density_during_mitosis.append(density_for_cell)
    
    return density_during_mitosis

def convert_to_floats_ignore_none(cells_obs_exp, cells_vol_exp, tau):
    mit_times, no_cells = cc.mitotic_times(cells_vol_exp, tau), len(cells_obs_exp)
    observations = [cells_obs_exp[n][1][mit_times[0][n]] if mit_times[0][n] < len(cells_obs_exp[n][1]) else None for n in range(no_cells)]
    
    # Filter out None values and convert the rest to floats
    float_observations = [float(obs) for obs in observations if obs is not None]
    
    return float_observations
## density at EOS

def cells_obs_at_EOS(cells_obs_exp):
    no_cells = len(cells_obs_exp)
    time_window = [0, 5]
    signa_delta_minus, signa_delta_plus = 1, 0
    # Assuming cc.observable_slope_change_times_ratio_method returns appropriate structures
    area_slope_change_times, area_slope_change_frame, delta_a_minus_at_slope_change, delta_a_plus_at_slope_change, der_ratio_at_slope_change = cc.observable_slope_change_times_ratio_method(cells_area, tau_slope_change, time_window, signa_delta_minus, signa_delta_plus)
    
    # Assuming area_slope_change_frame and area_slope_change_times are correctly returned
    spreadingphase_times = (area_slope_change_frame, area_slope_change_times)
    
    # Extract observable values at spreadingphase_times for each cell
    return [cells_obs_exp[n][1][spreadingphase_times[0][n]] for n in range(no_cells)]

# Load the data
#massfile = 'mass_tracks.dat'
#densityfile = 'density_tracks.dat'




##########DEFINESHELLPARS##########

# DEFINE PARSER
parser = argparse.ArgumentParser(description="Analysis of coupling mass growth and density: you need to input 4 files and 2 parameters")

# DEFINE ARGUMENTS
parser.add_argument("massfile", type=str, help="The mass input file (e.g., mass_tracks.dat)")
parser.add_argument("volfile", type=str, help="The volume input file (e.g., vol_tracks.dat)")
parser.add_argument("areafile", type=str, help="The area input file (e.g., area_tracks.dat)")
parser.add_argument("densityfile", type=str, help="The density input file (e.g., density_tracks.dat)")
parser.add_argument("--tau", type=float, default=1, help="Derivative time step (h)")

# UNPACK ARGUMENTS
args = parser.parse_args()      

###############LOADDATA####################

massfile, volfile, areafile, densityfile = args.massfile, args.volfile, args.areafile, args.densityfile

# BASIC TIMESERIES
cells_mass = cc.cell_cycle_time_series(massfile)
cells_vol = cc.cell_cycle_time_series(volfile)
cells_area = cc.cell_cycle_time_series(areafile)
cells_density = cc.cell_cycle_time_series(densityfile)
cells_id = cc.cell_index_original_file(massfile)
no_cells = len(cells_mass)
#print(len(cells_mass), len(cells_vol), len(cells_area), len(cells_density))

# BASIC PARAMETERS FOR DERIVATIVES
dt = 0.25
phase_step = 0.01
tau_gem = 2
tau_slope_change = 0.75
tau = 0.25 ## tau has got to take the form k*0.25 because of the discrete nature of the experiment

# RELEVANT CELL CYCLE TIMES
division_times = cc.division_times(cells_mass)
mitotic_times = cc.mitotic_times(cells_vol, tau)
time_window, signa_delta_minus, signa_delta_plus = [0, 5], 1, 0
area_slope_change_times, area_slope_change_frame, delta_a_minus_at_slope_change, delta_a_plus_at_slope_change, der_ratio_at_slope_change = cc.observable_slope_change_times_ratio_method(cells_area, tau_slope_change,  
                                                                                                                                time_window, signa_delta_minus, signa_delta_plus)
spreadingphase_times = [area_slope_change_frame, area_slope_change_times]
print (spreadingphase_times)
######################OBSERVABLES#######################
cells_mass, birth_times = cell_cycle_time_series(massfile)
cells_density, _ = cell_cycle_time_series(densityfile)  # We don't need birth times from here

# Find densities at birth times
densities_at_birth = density_at_birth_times(cells_density, birth_times)
## BIRTH
density_stat_around_birth = cc.value_in_time_range(cells_density, no_cells*[0], no_cells*[1])
cv_density_around_birth = stat.stdev(density_stat_around_birth)/stat.mean(density_stat_around_birth)

#print("Densities at birth times:", densities_at_birth)
# DENSITY TRACK AFTER SPREADING AND BEFORE MITOSIS
#cells_density_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_density, spreadingphase_times[1], mitotic_times[1]))
density_at_EOS= cells_obs_at_EOS(cells_density)
## EOS
density_stat_around_eos = cc.value_in_time_range(cells_density, np.subtract(area_slope_change_times, 0.5), np.add(area_slope_change_times, 0.5))
cv_density_around_eos = stat.stdev(density_stat_around_eos)/stat.mean(density_stat_around_eos)


#Mitosis
## MITOSIS
density_stat_around_mitosis = cc.value_in_time_range(cells_density, np.subtract(mitotic_times[1], 1), mitotic_times[1])
cv_density_around_mitosis = stat.stdev(density_stat_around_mitosis)/stat.mean(density_stat_around_mitosis)

density_during_mitosis = convert_to_floats_ignore_none(cells_density, cells_vol, tau)

#density_at_EOS = cells_obs_at_spreading_phase(cells_density, spreadingphase_times)

# DISTRIBUTIONS
#hist_density_as_bm = cc.all_values(cells_density_as_bm)
hist_density_birth = densities_at_birth
hist_density_birth_new = density_stat_around_birth
hist_density_mitosis = density_during_mitosis
hist_density_EOS= density_at_EOS
hist_density_mitosis_new = density_stat_around_mitosis
hist_density_EOS_new = density_stat_around_eos
# STATISTICS
#mean_density_as_bm = stat.mean(hist_density_as_bm)
#stdev_density_as_bm = stat.stdev(hist_density_as_bm)
#cv_density_as_bm = stdev_density_as_bm/mean_density_as_bm
mean_density_birth = stat.mean(hist_density_birth)
stdev_density_birth = stat.stdev(hist_density_birth)
cv_density_birth = stdev_density_birth/mean_density_birth
mean_density_birth_new= stat.mean(density_stat_around_birth)

mean_density_mitosis = stat.mean(hist_density_mitosis)
stdev_density_mitosis = stat.stdev(hist_density_mitosis)
cv_density_mitosis = stdev_density_mitosis/mean_density_mitosis
mean_density_mitosis_new= stat.mean(density_stat_around_mitosis)


mean_density_EOS = stat.mean(hist_density_EOS)
stdev_density_EOS = stat.stdev(hist_density_EOS)
cv_density_EOS = stdev_density_EOS/ mean_density_EOS
mean_density_EOS_new = stat.mean(density_stat_around_eos)
######################PRINTING#######################

print('density CV mitosis = %f' % cv_density_mitosis)
print('density CV birth= %f' % cv_density_birth)
print('density CV EOS = %f' % cv_density_EOS)
print('density CV birth new= %f' % cv_density_around_birth)
print('density CV EOS new = %f' % cv_density_around_eos)
print('density CV mitosis new = %f' % cv_density_around_mitosis)
#print_to_file('HeLa_pooled_mass_grate_as_bm_allvalues_derstep%.2f.dat' % tau, [hist_grate_mass_as_bm], 1)
#print_to_file('HeLa_pooled_vol_grate_as_bm_allvalues_derstep%.2f.dat' % tau, [hist_grate_vol_as_bm], 1)

######################PLOTTING#######################

# COLORS

cool_purple = '#810f7c'
purple_mitosis = '#3c073a'
purple_EOS = '#c617be'
# PLOT HISTOGRAM
width = 2
height = width/1.618
#fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99) 
#for axis in ['top','bottom','left','right']:
    #ax.spines[axis].set_linewidth(0.8)    
#ax.hist(x=hist_density_as_bm, bins=6, density=True, color=cool_purple)
#ax.axvline(x=mean_density_as_bm,  linewidth=0.8, linestyle='--', color='black', 
           #label=r'avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_as_bm)  
#ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=6)
#ax.set_ylabel(r'pdf', fontsize=6)
#ax.tick_params(axis='x', length=4, width=0.6, labelsize=6)
#ax.tick_params(axis='y', length=4, width=0.6, labelsize=6)
#ax.legend(fontsize=6, loc=(0.6, 0.8))
#fig.savefig('Hela_density_as_bm_hist.pdf', bbox_inches='tight')  
#plt.close(fig) 

# PLOT DENSITY
#width = 2
#height = width/1.618
#fig, ax = plt.subplots(figsize=(width, height), dpi=100)
#fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99) 
#for axis in ['top','bottom','left','right']:
    #ax.spines[axis].set_linewidth(0.8)    
#sns.kdeplot(hist_density_as_bm, color=cool_purple, fill=True)
#ax.axvline(x=mean_density_as_bm,  linewidth=0.8, linestyle='--', color='black', 
           #label=r'avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_as_bm)  
#ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=6)
#ax.set_ylabel(r'pdf', fontsize=6)
#ax.tick_params(axis='x', length=4, width=0.6, labelsize=6)
#ax.tick_params(axis='y', length=4, width=0.6, labelsize=6)
#ax.legend(fontsize=6, loc=(0.6, 0.8))
#fig.savefig('Hela_density_as_bm_density.pdf', bbox_inches='tight') 
#plt.close(fig)

width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99) 
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.8)    
sns.kdeplot(hist_density_birth, color=cool_purple, fill=True)
ax.axvline(x=mean_density_birth,  linewidth=0.8, linestyle='--', color='black', 
           label=r'avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_birth)  
ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=6)
ax.set_ylabel(r'pdf', fontsize=6)
ax.tick_params(axis='x', length=4, width=0.6, labelsize=6)
ax.tick_params(axis='y', length=4, width=0.6, labelsize=6)
ax.legend(fontsize=6, loc=(0.6, 0.8))
fig.savefig('Hela_density_as_birth.pdf', bbox_inches='tight') 
plt.close(fig)



#without legend
cool_blue_current = '#8cc5e3'
cool_orange = '#FF9F00'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
width = 6.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)

# Plot KDE for birth density
sns.kdeplot(hist_density_birth, color=cool_purple, fill=True, label='Density at birth')
ax.axvline(x=mean_density_birth, linewidth=1, linestyle='--', color='black',
           label=r'Birth avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_birth)

# Plot KDE for mitosis density
sns.kdeplot(hist_density_mitosis, color=purple_mitosis, fill=True, label='Density at mitosis')
ax.axvline(x=mean_density_mitosis, linewidth=1, linestyle='--', color='black',
           label=r'Mitosis avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_mitosis)


# Plot kde for EOS density
sns.kdeplot(hist_density_EOS, color=purple_EOS, fill=True, label='Density at EOS')
ax.axvline(x=mean_density_EOS, linewidth=1, linestyle='--', color='black',
           label=r'EOS avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_EOS)

ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=20)
ax.set_ylabel(r'pdf', fontsize=20)
ax.tick_params(axis='x', length=4, width=0.6, labelsize=16)
ax.tick_params(axis='y', length=4, width=0.6, labelsize=16)
ax.legend(fontsize=10, loc='upper right')

fig.savefig('Hela_density_birth_mitosis_EOS.pdf', bbox_inches='tight')
plt.close(fig)

scaling_factor=1.3
width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)

# Plot KDE for birth density
sns.kdeplot(hist_density_birth, color=cool_purple, fill=True, label='Density at birth')
ax.axvline(x=mean_density_birth, linewidth=0.8, linestyle='--', color='black',
           label=r'Birth avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_birth)

# Plot KDE for birth new density
sns.kdeplot(hist_density_birth_new, color=purple_mitosis, fill=True, label='Density at birth new')
ax.axvline(x=mean_density_birth_new, linewidth=0.8, linestyle='--', color='black',
           label=r'birth avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_birth_new)
ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=10)
ax.set_ylabel(r'pdf', fontsize=10)
ax.tick_params(axis='x', length=4, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=4, width=0.6, labelsize=8)
ax.legend(fontsize=2, loc='right')
fig.savefig('Hela_density_birth_birthnew.pdf', bbox_inches='tight')
plt.close(fig)


scaling_factor=1.3
width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)

# Plot KDE for mitosis density
sns.kdeplot(hist_density_mitosis, color=cool_purple, fill=True, label='Density at mitosis')
ax.axvline(x=mean_density_mitosis, linewidth=0.8, linestyle='--', color='black',
           label=r'mitosis avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_mitosis)

# Plot KDE for mitosis new density
sns.kdeplot(hist_density_mitosis_new, color=purple_mitosis, fill=True, label='Density at mitosis new')
ax.axvline(x=mean_density_mitosis_new, linewidth=0.8, linestyle='--', color='black',
           label=r'Mitosis avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_mitosis_new)
ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=10)
ax.set_ylabel(r'pdf', fontsize=10)
ax.tick_params(axis='x', length=4, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=4, width=0.6, labelsize=8)
ax.legend(fontsize=2, loc='right')
fig.savefig('Hela_density_mitosis_mitosisnew.pdf', bbox_inches='tight')
plt.close(fig)


scaling_factor=1.3
width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)




# Plot kde for EOS density
sns.kdeplot(hist_density_EOS, color=purple_EOS, fill=True, label='Density at EOS')
ax.axvline(x=mean_density_EOS, linewidth=0.8, linestyle='--', color='black',
           label=r'EOS avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_EOS)
# Plot KDE for mitosis density
sns.kdeplot(hist_density_EOS_new, color=purple_mitosis, fill=True, label='Density at EOS new')
ax.axvline(x=mean_density_EOS_new, linewidth=0.8, linestyle='--', color='black',
           label=r'EOS new avg = %.2f pg $\mathrm{\mu m}^{-3}$' % mean_density_EOS_new)

ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=10)
ax.set_ylabel(r'pdf', fontsize=10)
ax.tick_params(axis='x', length=4, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=4, width=0.6, labelsize=8)
ax.legend(fontsize=2, loc='right')

fig.savefig('Hela_density_EOS_EOSnew.pdf', bbox_inches='tight')
plt.close(fig)
