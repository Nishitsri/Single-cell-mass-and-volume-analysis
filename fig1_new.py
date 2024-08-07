import sys
import os
import argparse

#currentdir = os.path.dirname(os.path.realpath(__file__))
#parentdir = os.path.dirname(currentdir)
#librarydir = os.path.dirname(parentdir) + '/0_libraries'
#sys.path.append(librarydir)

import llib_customfunction
import llib_figures as figs
import llib_cellcycles as cc
import matplotlib.pyplot as plt
import numpy as np
import statistics as stat
from scipy.stats import gaussian_kde
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from math import nan
from math import sqrt

#plt.style.use('mystyle')

###########FUNCTIONS################

def mean_stdev_cv_vs_time(cells_obs, dt):
    time_obs, obs_stat_at_fixed_time = cc.obs_stat_at_fixed_time(cells_obs, dt)
    mean_obs_at_fixed_time = [ stat.mean(obs_stat_at_fixed_time[k]) for k in range(len(time_obs))]
    stdev_obs_at_fixed_time = [ stat.stdev(obs_stat_at_fixed_time[k]) if len(obs_stat_at_fixed_time[k]) > 1 else None for k in range(len(time_obs)) ]
    cv_obs_at_fixed_time = [ stdev_obs_at_fixed_time[k]/mean_obs_at_fixed_time[k] if len(obs_stat_at_fixed_time[k]) > 1 else None for k in range(len(time_obs)) ]

    return mean_obs_at_fixed_time, stdev_obs_at_fixed_time, cv_obs_at_fixed_time, time_obs

def mean_stdev_cv_vs_phase(cells_obs, dphase):
    phase_obs, obs_stat_at_fixed_phase = cc.obs_stat_at_fixed_phase(cells_obs, dphase)
    mean_obs_at_fixed_phase = [ stat.mean(obs_stat_at_fixed_phase[k]) for k in range(len(phase_obs))]
    stdev_obs_at_fixed_phase = [ stat.stdev(obs_stat_at_fixed_phase[k]) if len(obs_stat_at_fixed_phase[k]) > 1 else None for k in range(len(phase_obs)) ]
    cv_obs_at_fixed_phase = [ stdev_obs_at_fixed_phase[k]/mean_obs_at_fixed_phase[k] if len(obs_stat_at_fixed_phase[k]) > 1 else None for k in range(len(phase_obs)) ]

    return mean_obs_at_fixed_phase, stdev_obs_at_fixed_phase, cv_obs_at_fixed_phase, phase_obs    

##########DEFINESHELLPARS##########

# DEFINE PARSER
parser = argparse.ArgumentParser(description="Analysis of coupling mass growth and density: you need to input 4 files and 2 parameters")

# DEFINE ARGUMENTS
parser.add_argument("massfile", type=str, help="The mass input file (e.g., mass_tracks.dat)")
parser.add_argument("volfile", type=str, help="The volume input file (e.g., vol_tracks.dat)")
parser.add_argument("areafile", type=str, help="The area input file (e.g., area_tracks.dat)")
parser.add_argument("densityfile", type=str, help="The density input file (e.g., density_tracks.dat)")

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
tau = 1 ## tau has got to take the form k*0.25 because of the discrete nature of the experiment

# RELEVANT CELL CYCLE TIMES
division_times = cc.division_times(cells_mass)
mitotic_times = cc.mitotic_times(cells_vol, tau)
time_window, signa_delta_minus, signa_delta_plus = [0, 5], 1, 0
area_slope_change_times, area_slope_change_frame, delta_a_minus_at_slope_change, delta_a_plus_at_slope_change, der_ratio_at_slope_change = cc.observable_slope_change_times_ratio_method(cells_area, tau_slope_change,	
																																time_window, signa_delta_minus, signa_delta_plus)
spreadingphase_times = [area_slope_change_frame, area_slope_change_times]

######################OBSERVABLES#######################

# TRACKS

## COMPUTE NORMALIZED TRACKS
cells_mass_normalized = cc.normalize_cell_cycle_time_series(cells_mass)
cells_vol_normalized = cc.normalize_cell_cycle_time_series(cells_vol)
cells_density_normalized = cc.normalize_cell_cycle_time_series(cells_density)
cells_area_normalized = cc.normalize_cell_cycle_time_series(cells_area)

## COMPUTE TRACKS IN CELL CYCLE TIME
cells_mass_cc_time = [ [ list(np.divide(cells_mass[n][0], cells_mass[n][0][len(cells_mass[n][0])-1])), cells_mass[n][1]] for n in range(no_cells) ]
cells_vol_cc_time = [ [ list(np.divide(cells_vol[n][0], cells_vol[n][0][len(cells_vol[n][0])-1])), cells_vol[n][1]] for n in range(no_cells) ]
cells_density_cc_time = [ [ list(np.divide(cells_density[n][0], cells_density[n][0][len(cells_density[n][0])-1])), cells_density[n][1]] for n in range(no_cells) ]

## COMPUTE NORMALIZED TRACKS IN CELL CYCLE TIME
cells_mass_normalized_cc_time = [ [ list(np.divide(cells_mass_normalized[n][0], cells_mass_normalized[n][0][len(cells_mass_normalized[n][0])-1])), cells_mass_normalized[n][1]] for n in range(no_cells) ]
cells_vol_normalized_cc_time = [ [ list(np.divide(cells_vol_normalized[n][0], cells_vol_normalized[n][0][len(cells_vol_normalized[n][0])-1])), cells_vol_normalized[n][1]] for n in range(no_cells) ]
cells_density_normalized_cc_time = [ [ list(np.divide(cells_density_normalized[n][0], cells_density_normalized[n][0][len(cells_density_normalized[n][0])-1])), cells_density_normalized[n][1]] for n in range(no_cells) ]
cells_area_cc_time = [ [ list(np.divide(cells_area[n][0], cells_area[n][0][len(cells_area[n][0])-1])), cells_area[n][1]] for n in range(no_cells) ]


######################STATISTICS#######################

dt_window = 1
dphase_window = 0.04

# DIMENSIONFUL

## FIXED TIME
mean_mass_at_fixed_time, stdev_mass_at_fixed_time, cv_mass_at_fixed_time, time_mass = mean_stdev_cv_vs_time(cells_mass, dt_window)
mean_volume_at_fixed_time, stdev_volume_at_fixed_time, cv_volume_at_fixed_time, time_volume = mean_stdev_cv_vs_time(cells_vol, dt_window)
mean_density_at_fixed_time, stdev_density_at_fixed_time, cv_density_at_fixed_time, time_density = mean_stdev_cv_vs_time(cells_density, dt_window)
mean_area_at_fixed_time, stdev_area_at_fixed_time, cv_area_at_fixed_time, time_area = mean_stdev_cv_vs_time(cells_area, dt_window)

## CELL CYCLE TIME
mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase, cv_mass_at_fixed_phase, phase_mass = mean_stdev_cv_vs_phase(cells_mass, dphase_window)
mean_volume_at_fixed_phase, stdev_volume_at_fixed_phase, cv_volume_at_fixed_phase, phase_volume = mean_stdev_cv_vs_phase(cells_vol, dphase_window)
mean_density_at_fixed_phase, stdev_density_at_fixed_phase, cv_density_at_fixed_phase, phase_density = mean_stdev_cv_vs_phase(cells_density, dphase_window)
mean_area_at_fixed_phase, stdev_area_at_fixed_phase, cv_area_at_fixed_phase, phase_area = mean_stdev_cv_vs_phase(cells_area, dphase_window)

# DIMENSIONLESS

## NORMALIZED
mean_mass_normalized_at_fixed_time, stdev_mass_normalized_at_fixed_time, cv_mass_normalized_at_fixed_time, time_mass_normalized = mean_stdev_cv_vs_time(cells_mass_normalized, dt_window)
mean_volume_normalized_at_fixed_time, stdev_volume_normalized_at_fixed_time, cv_volume_normalized_at_fixed_time, time_volume_normalized = mean_stdev_cv_vs_time(cells_vol_normalized, dt_window)
mean_density_normalized_at_fixed_time, stdev_density_normalized_at_fixed_time, cv_density_normalized_at_fixed_time, time_density_normalized = mean_stdev_cv_vs_time(cells_density_normalized, dt_window)
mean_area_normalized_at_fixed_time, stdev_area_normalized_at_fixed_time, cv_area_normalized_at_fixed_time, time_area_normalized = mean_stdev_cv_vs_time(cells_area_normalized, dt_window)

## CELL CYCLE TIME
mean_mass_normalized_at_fixed_phase, stdev_mass_normalized_at_fixed_phase, cv_mass_normalized_at_fixed_phase, phase_mass_normalized = mean_stdev_cv_vs_phase(cells_mass_normalized, dphase_window)
mean_volume_normalized_at_fixed_phase, stdev_volume_normalized_at_fixed_phase, cv_volume_normalized_at_fixed_phase, phase_volume_normalized = mean_stdev_cv_vs_phase(cells_vol_normalized, dphase_window)
mean_density_normalized_at_fixed_phase, stdev_density_normalized_at_fixed_phase, cv_density_normalized_at_fixed_phase, phase_density_normalized = mean_stdev_cv_vs_phase(cells_density_normalized, dphase_window)
mean_area_normalized_at_fixed_phase, stdev_area_normalized_at_fixed_phase, cv_area_normalized_at_fixed_phase, phase_area_normalized = mean_stdev_cv_vs_phase(cells_area_normalized, dphase_window)

######################PRINTING#######################

######################PLOTTING#######################

# COLORS
cool_red = '#e41a1c'
#cool_red_outer = '#a90908'
cool_blue = '#8cc5e3'
cool_purple = '#984ea3'
cool_green = '#4daf4a'
cool_red_old = '#a50f15'
cool_blue_old = '#253494'
cool_purple_old = '#810f7c'
cool_green_old = '#1b9e77'
#
#
cool_orange = '#FF9F00'
cool_purple_new = '#4a2377'
# AVERAGES VS CELL CYCLE TIME

#width = 1.5
#height = width/1.618
#plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['ps.fonttype'] = 42
#plt.rcParams['font.family'] = ['sans-serif']
#plt.rcParams['font.sans-serif'] = ['Arial']
#fig, ax = plt.subplots(figsize=(width, height), dpi=100)
#fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
#for axis in ['top','bottom','left','right']:
    #ax.spines[axis].set_linewidth(1.0)
    #ax.spines[['right', 'top']].set_visible(False)
#ax.plot(phase_mass, mean_mass_at_fixed_phase, 'o', color=cool_red, markeredgecolor='black', markeredgewidth='0.3', markersize=4, alpha=1)
#ax.fill_between(phase_mass, np.subtract(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), 
                #np.add(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), color='red', alpha=0.2) 
#ax.set_xlabel(r'cell cycle time $t/\tau_{div}$ ([0, 1])', fontsize=10)
#ax.set_xlabel(r'cell cycle phase', fontsize=10)
#ax.set_ylabel('mean mass\n' r'$\langle M \rangle$ (pg)', fontsize=10)
#ax.set_yticks([200, 300, 400, 500])
#ax.set_yticks([200, 300, 400])
#ax.tick_params(axis='x', length=2, width=0.6, labelsize=6)
#ax.tick_params(axis='y', length=2, width=0.6, labelsize=6)
#fig.savefig('Hela_pooled_mean_mass_vs_cell_cycle_phase.svg', bbox_inches='tight')  
#plt.close(fig)












scaling_factor=1.3
width = 2
height = width/1.618
#plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['Arial']
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_mass, mean_mass_at_fixed_phase, 'o', color=cool_orange, markeredgecolor= 'black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_mass, np.subtract(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), 
                np.add(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), color='orange', alpha=0.2) 
#ax.set_xlabel(r'cell cycle time $t/\tau_{div}$ ([0, 1])', fontsize=10)
ax.set_xlabel(r'cell cycle phase', fontsize=12)
ax.set_ylabel('mean mass\n' r'$\langle M \rangle$ (pg)', fontsize=12)
#ax.set_yticks([200, 300, 400, 500])
ax.set_yticks([200, 300, 400])
ax.tick_params(axis='x', length=2, width=0.6, labelsize=10)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=10)
fig.savefig('Hela_pooled_mean_mass_vs_cell_cycle_phase_coolorangesmallerdots_blackborder.pdf', bbox_inches='tight')  
plt.close(fig)


## used cool_blue for making volume plots

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
ax.plot(phase_volume, mean_volume_at_fixed_phase, 'o', color=cool_blue, markeredgecolor='black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_volume, np.subtract(mean_volume_at_fixed_phase, stdev_volume_at_fixed_phase), 
                np.add(mean_volume_at_fixed_phase, stdev_volume_at_fixed_phase), color='blue', alpha=0.2) 
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean volume\n'r'$\langle V \rangle$ ($\mathrm{\mu m}^{3}$)', fontsize=10)
ax.set_yticks([1500, 2500, 3500])
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_volume_vs_cell_cycle_phase_coolbluewnew_blackborder.pdf', bbox_inches='tight')  
plt.close(fig)




## used cool_blue old for making volume plots

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
ax.plot(phase_volume, mean_volume_at_fixed_phase, 'o', color=cool_blue_old, markeredgecolor='black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_volume, np.subtract(mean_volume_at_fixed_phase, stdev_volume_at_fixed_phase), 
                np.add(mean_volume_at_fixed_phase, stdev_volume_at_fixed_phase), color='blue', alpha=0.2) 
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean volume\n'r'$\langle V \rangle$ ($\mathrm{\mu m}^{3}$)', fontsize=10)
ax.set_yticks([1500, 2500, 3500])
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_volume_vs_cell_cycle_phase_coolbluewold_blackborder.pdf', bbox_inches='tight')  
plt.close(fig)

##used cool_green for plotting area

width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_area, mean_area_at_fixed_phase, 'o', color=cool_green, markeredgecolor='black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_area, np.subtract(mean_area_at_fixed_phase, stdev_area_at_fixed_phase), 
                np.add(mean_area_at_fixed_phase, stdev_area_at_fixed_phase), color='green', alpha=0.2) 
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean \n spreading \n area'r'$\langle A \rangle$ ($\mathrm{\mu m}^{2}$)', fontsize=10)
#ax.set_ylabel(r'$\langle A \rangle$ ($\mathrm{\mu m}^{2}$)''\nmean spreading area', fontsize=8)
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_area_vs_cell_cycle_phase_coolgreennew_blackborder.pdf', bbox_inches='tight')  
plt.close(fig) 



##used cool_green old for plotting area

width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_area, mean_area_at_fixed_phase, 'o', color=cool_green_old, markeredgecolor='black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_area, np.subtract(mean_area_at_fixed_phase, stdev_area_at_fixed_phase), 
                np.add(mean_area_at_fixed_phase, stdev_area_at_fixed_phase), color='green', alpha=0.2) 
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean \n spreading \n area'r'$\langle A \rangle$ ($\mathrm{\mu m}^{2}$)', fontsize=10)
#ax.set_ylabel(r'$\langle A \rangle$ ($\mathrm{\mu m}^{2}$)''\nmean spreading area', fontsize=8)
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_area_vs_cell_cycle_phase_coolgreenold_blackborder.pdf', bbox_inches='tight')  
plt.close(fig) 



width = 2
height = width/1.618
#plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = ['sans-serif']
plt.rcParams['font.sans-serif'] = ['Arial']
fig, ax1 = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(1.0)
    ax1.spines[['top']].set_visible(False)
ax1.plot(phase_density, mean_density_at_fixed_phase, 'o', color=cool_purple_new, markeredgecolor='black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax1.fill_between(phase_density, np.subtract(mean_density_at_fixed_phase, stdev_density_at_fixed_phase), 
                np.add(mean_density_at_fixed_phase, stdev_density_at_fixed_phase), color='purple', alpha=0.2) 
ax1.set_xlabel(r'cell cycle phase', fontsize=8)
ax1.set_ylabel('mean density\n' r'$\langle \rho \rangle$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=8)
ax1.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax1.tick_params(axis='y', length=2, width=0.6, labelsize=8)
ax2 = ax1.twinx() 
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(1.0)
    ax2.spines[['top']].set_visible(False)
ax2.plot(phase_mass, mean_mass_at_fixed_phase, 'o', color=cool_orange, markeredgecolor='black', markeredgewidth='0.2', markersize=3, alpha=1)
ax2.fill_between(phase_mass, np.subtract(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), 
                np.add(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), color='orange', alpha=0.2) 
#ax.set_xlabel(r'cell cycle time $t/\tau_{div}$ ([0, 1])', fontsize=10)
ax2.set_xlabel(r'cell cycle phase', fontsize=10)
ax2.set_ylabel('mean mass\n' r'$\langle M \rangle$ (pg)', fontsize=10)
#ax.set_yticks([200, 300, 400, 500])
ax2.set_yticks([200, 300, 400])
ax2.tick_params(axis='x', length=2, width=0.6, labelsize=6)
ax2.tick_params(axis='y', length=2, width=0.6, labelsize=6)
#fig.savefig('Hela_pooled_mean_mass_mean_density_vs_cell_cycle_phase_purple_orange.pdf', bbox_inches='tight')  
plt.close(fig)


## used cool_purple_old for making density plots
width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
#plt.rcParams['font.family'] = ['sans-serif']
#plt.rcParams['font.sans-serif'] = ['Arial']
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_density, mean_density_at_fixed_phase, 'o', color=cool_purple_old, markeredgecolor='black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_density, np.subtract(mean_density_at_fixed_phase, stdev_density_at_fixed_phase), 
                np.add(mean_density_at_fixed_phase, stdev_density_at_fixed_phase), color='purple', alpha=0.2) 
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean density\n' r'$\langle \rho \rangle$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=10)
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_density_vs_cell_cycle_phase_coolpurpleold_blackborder.pdf', bbox_inches='tight')  
plt.close(fig)

## used cool_purple for making density plots
width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
#plt.rcParams['font.family'] = ['sans-serif']
#plt.rcParams['font.sans-serif'] = ['Arial']
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_density, mean_density_at_fixed_phase, 'o', color=cool_purple, markeredgecolor='black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_density, np.subtract(mean_density_at_fixed_phase, stdev_density_at_fixed_phase), 
                np.add(mean_density_at_fixed_phase, stdev_density_at_fixed_phase), color='purple', alpha=0.2) 
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean density\n' r'$\langle \rho \rangle$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=10)
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_density_vs_cell_cycle_phase_coolpurple_blackborder.pdf', bbox_inches='tight')  
plt.close(fig)


##used cool_orange for mass


width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
#plt.rcParams['font.family'] = ['sans-serif']
#plt.rcParams['font.sans-serif'] = ['Arial']
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_mass, mean_mass_at_fixed_phase, 'o', color=cool_orange, markeredgecolor= 'black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_mass, np.subtract(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), 
                np.add(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), color='orange', alpha=0.2) 
#ax.set_xlabel(r'cell cycle time $t/\tau_{div}$ ([0, 1])', fontsize=10)
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean mass\n' r'$\langle M \rangle$ (pg)', fontsize=10)
#ax.set_yticks([200, 300, 400, 500])
ax.set_yticks([200, 300, 400])
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_mass_vs_cell_cycle_phase_coolorangesmallerdots_blackborder_editable.pdf', bbox_inches='tight')  
plt.close(fig)

##old red for mass


width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
#plt.rcParams['font.family'] = ['sans-serif']
#plt.rcParams['font.sans-serif'] = ['Arial']
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_mass, mean_mass_at_fixed_phase, 'o', color=cool_red, markeredgecolor= 'black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_mass, np.subtract(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), 
                np.add(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), color='red', alpha=0.2) 
#ax.set_xlabel(r'cell cycle time $t/\tau_{div}$ ([0, 1])', fontsize=10)
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean mass\n' r'$\langle M \rangle$ (pg)', fontsize=10)
#ax.set_yticks([200, 300, 400, 500])
ax.set_yticks([200, 300, 400])
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_mass_vs_cell_cycle_phase_coolred_blackborder_editable.pdf', bbox_inches='tight')  
plt.close(fig)

width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
#plt.rcParams['font.family'] = ['sans-serif']
#plt.rcParams['font.sans-serif'] = ['Arial']
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_mass, mean_mass_at_fixed_phase, 'o', color=cool_red_old, markeredgecolor= 'black', markeredgewidth='0.2', markersize=5*scaling_factor, alpha=1)
ax.fill_between(phase_mass, np.subtract(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), 
                np.add(mean_mass_at_fixed_phase, stdev_mass_at_fixed_phase), color='red', alpha=0.2) 
#ax.set_xlabel(r'cell cycle time $t/\tau_{div}$ ([0, 1])', fontsize=10)
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean mass\n' r'$\langle M \rangle$ (pg)', fontsize=10)
#ax.set_yticks([200, 300, 400, 500])
ax.set_yticks([200, 300, 400])
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_mean_mass_vs_cell_cycle_phase_coolredold_blackborder_editable.pdf', bbox_inches='tight')  
plt.close(fig)

