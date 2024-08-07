
# Do an EXPONENTIAL TEST on the data
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
    sterr_obs_at_fixed_time = [ stdev_obs_at_fixed_time[k]/sqrt(len(obs_stat_at_fixed_time[k])) if len(obs_stat_at_fixed_time[k]) > 1 else None for k in range(len(time_obs)) ]
    cv_obs_at_fixed_time = [ stdev_obs_at_fixed_time[k]/mean_obs_at_fixed_time[k] if len(obs_stat_at_fixed_time[k]) > 1 else None for k in range(len(time_obs)) ]

    return mean_obs_at_fixed_time, sterr_obs_at_fixed_time, stdev_obs_at_fixed_time, cv_obs_at_fixed_time, time_obs

def mean_stdev_cv_vs_phase(cells_obs, dphase):
    phase_obs, obs_stat_at_fixed_phase = cc.obs_stat_at_fixed_phase(cells_obs, dphase)
    mean_obs_at_fixed_phase = [ stat.mean(obs_stat_at_fixed_phase[k]) for k in range(len(phase_obs))]
    stdev_obs_at_fixed_phase = [ stat.stdev(obs_stat_at_fixed_phase[k]) if len(obs_stat_at_fixed_phase[k]) > 1 else None for k in range(len(phase_obs)) ]
    sterr_obs_at_fixed_phase = [ stdev_obs_at_fixed_phase[k]/sqrt(len(obs_stat_at_fixed_phase[k])) if len(obs_stat_at_fixed_phase[k]) > 1 else None for k in range(len(phase_obs)) ]
    cv_obs_at_fixed_phase = [ stdev_obs_at_fixed_phase[k]/mean_obs_at_fixed_phase[k] if len(obs_stat_at_fixed_phase[k]) > 1 else None for k in range(len(phase_obs)) ]

    return mean_obs_at_fixed_phase, sterr_obs_at_fixed_phase, stdev_obs_at_fixed_phase, cv_obs_at_fixed_phase, phase_obs   

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
tau = 0.25 ## tau has got to take the form k*0.25 because of the discrete nature of the experiment

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

## COMPUTE TRACKS IN CELL CYCLE TIME
cells_mass_cc_time = [ [ list(np.divide(cells_mass[n][0], cells_mass[n][0][len(cells_mass[n][0])-1])), cells_mass[n][1]] for n in range(no_cells) ]
cells_vol_cc_time = [ [ list(np.divide(cells_vol[n][0], cells_vol[n][0][len(cells_vol[n][0])-1])), cells_vol[n][1]] for n in range(no_cells) ]
cells_density_cc_time = [ [ list(np.divide(cells_density[n][0], cells_density[n][0][len(cells_density[n][0])-1])), cells_density[n][1]] for n in range(no_cells) ]

## COMPUTE NORMALIZED TRACKS IN CELL CYCLE TIME
cells_mass_normalized_cc_time = [ [ list(np.divide(cells_mass_normalized[n][0], cells_mass_normalized[n][0][len(cells_mass_normalized[n][0])-1])), cells_mass_normalized[n][1]] for n in range(no_cells) ]
cells_vol_normalized_cc_time = [ [ list(np.divide(cells_vol_normalized[n][0], cells_vol_normalized[n][0][len(cells_vol_normalized[n][0])-1])), cells_vol_normalized[n][1]] for n in range(no_cells) ]
cells_density_normalized_cc_time = [ [ list(np.divide(cells_density_normalized[n][0], cells_density_normalized[n][0][len(cells_density_normalized[n][0])-1])), cells_density_normalized[n][1]] for n in range(no_cells) ]

# DERIVATIVE TIME STEP
tau = 1

# COMPUTE GROWTH RATE
cells_grate_mass = cc.growth_rate_linear_fit(cells_mass, tau)
cells_grate_vol = cc.growth_rate_linear_fit(cells_vol, tau)

######################STATISTICS#######################

# MASS, VOLUME AND DENSITY
dt_window = 1
dphase_window = 0.03

## CELL CYCLE TIME
mean_mass_at_fixed_phase, sterr_mass_at_fixed_time, stdev_mass_at_fixed_phase, cv_mass_at_fixed_phase, phase_mass = mean_stdev_cv_vs_phase(cells_mass, dphase_window)
mean_volume_at_fixed_phase, sterr_volume_at_fixed_time, stdev_volume_at_fixed_phase, cv_volume_at_fixed_phase, phase_volume = mean_stdev_cv_vs_phase(cells_vol, dphase_window)
mean_density_at_fixed_phase, sterr_density_at_fixed_time, stdev_density_at_fixed_phase, cv_density_at_fixed_phase, phase_density = mean_stdev_cv_vs_phase(cells_density, dphase_window)

# DIMENSIONLESS

## CELL CYCLE TIME
mean_mass_normalized_at_fixed_phase, sterr_mass_normalized_at_fixed_phase, stdev_mass_normalized_at_fixed_phase, cv_mass_normalized_at_fixed_phase, phase_mass_normalized = mean_stdev_cv_vs_phase(cells_mass_normalized, dphase_window)
mean_volume_normalized_at_fixed_phase, sterr_volume_normalized_at_fixed_phase, sstdev_volume_normalized_at_fixed_phase, cv_volume_normalized_at_fixed_phase, phase_volume_normalized = mean_stdev_cv_vs_phase(cells_vol_normalized, dphase_window)
mean_density_normalized_at_fixed_phase, sterr_density_normalized_at_fixed_phase, sstdev_density_normalized_at_fixed_phase, cv_density_normalized_at_fixed_phase, phase_density_normalized = mean_stdev_cv_vs_phase(cells_density_normalized, dphase_window)

# GROWTH RATES
dt_window = 1
dphase_window = 0.01

## MASS GROWTH RATE
mean_grate_mass_at_fixed_time, sterr_grate_mass_at_fixed_time, stdev_grate_mass_at_fixed_time, cv_grate_mass_at_fixed_time, time_grate_mass = mean_stdev_cv_vs_time(cells_grate_mass, dt_window)
mean_grate_mass_at_fixed_phase, sterr_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase, cv_grate_mass_at_fixed_phase, phase_grate_mass = mean_stdev_cv_vs_phase(cells_grate_mass, dphase_window)

## VOL GROWTH RATE
mean_grate_vol_at_fixed_time, sterr_grate_vol_at_fixed_time, stdev_grate_vol_at_fixed_time, cv_grate_vol_at_fixed_time, time_grate_vol = mean_stdev_cv_vs_time(cells_grate_vol, dt_window)
mean_grate_vol_at_fixed_phase, sterr_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase, cv_grate_vol_at_fixed_phase, phase_grate_vol = mean_stdev_cv_vs_phase(cells_grate_vol, dphase_window)

## CONDITIONED MEAN BINNING WINDOW SIZE
window_size = 800

## EXPONENTIAL TEST
cells_mass_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_mass, spreadingphase_times[1], mitotic_times[1]))
cells_vol_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_vol, spreadingphase_times[1], mitotic_times[1]))

cells_der_mass_as_bm = cc.der_linear_fit(cells_mass_as_bm, tau)
cells_der_vol_as_bm = cc.der_linear_fit(cells_vol_as_bm, tau)

dermass_vs_mass, corr_dermass_vs_mass = cc.scatter_plot_vectors(cells_mass_as_bm, cells_der_mass_as_bm)
binned_dermass_vs_mass = cc.binningdata_bins_mean(dermass_vs_mass[0], dermass_vs_mass[1], window_size)

dervol_vs_vol, corr_dervol_vs_vol = cc.scatter_plot_vectors(cells_vol_as_bm, cells_der_vol_as_bm)
binned_dervol_vs_vol = cc.binningdata_bins_mean(dervol_vs_vol[0], dervol_vs_vol[1], window_size)

# GROWTH RATES CORRELATION
cells_grate_mass_as_bm = cc.growth_rate_linear_fit(cells_mass_as_bm, tau)
cells_grate_vol_as_bm = cc.growth_rate_linear_fit(cells_vol_as_bm, tau)

gratemass_vs_gratevol, corr_gratemass_vs_gratevol = cc.scatter_plot_vectors(cells_grate_mass_as_bm, cells_grate_vol_as_bm)
binned_gratemass_X_vs_gratevol_Y = cc.binningdata_bins_mean(gratemass_vs_gratevol[0], gratemass_vs_gratevol[1], window_size)
binned_gratevol_Y_vs_gratemass_X = cc.binningdata_bins_mean(gratemass_vs_gratevol[1], gratemass_vs_gratevol[0], window_size)

decile_filter = 5
filter_window_size = 500
filter_gratemass_vs_gratevol = filter_data(gratemass_vs_gratevol, decile_filter)
filter_binned_gratemass_X_vs_gratevol_Y = cc.binningdata_bins_mean(filter_gratemass_vs_gratevol[0], filter_gratemass_vs_gratevol[1], filter_window_size)
filter_binned_gratevol_Y_vs_gratemass_X = cc.binningdata_bins_mean(filter_gratemass_vs_gratevol[1], filter_gratemass_vs_gratevol[0], filter_window_size)

######################PRINTING#######################
# COLORS

cool_red = '#a50f15'
cool_blue = '#253494'
cool_purple = '#810f7c'
cool_green = '#74c476'

# COLORS

# COLORS
cool_red = '#e41a1c'
#cool_red_outer = '#a90908'
cool_blue = '#8cc5e3'
cool_purple = '#810f7c'
cool_green = '#4daf4a'
cool_red_old = '#a50f15'
cool_blue_old = '#253494'
cool_purple_old = '#810f7c'
cool_green_old = '#1b9e77'
#
#
cool_orange = '#FF9F00'
cool_purple_new = '#4a2377'
######################PLOTTING#######################

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
ax.plot(binned_dermass_vs_mass[0], binned_dermass_vs_mass[1], 'o', color=cool_orange, markeredgecolor='black', markeredgewidth='0.3', markersize=8,  alpha=1) 
#ax.errorbar(binned_dermass_vs_mass[0], binned_dermass_vs_mass[1], yerr=binned_dermass_vs_mass[5],
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1) 
ax.set_xlabel('mass' r'(pg)', fontsize=10)
ax.set_ylabel('mean derivative\n' r'$\langle \frac{dM}{dt} \rangle$ (pg $\mathrm{h}^{-1}$)', fontsize=10)
ax.set_ylim(5, 15)
ax.tick_params(axis='x', length=2, width=0.6, labelsize=8)
ax.tick_params(axis='y', length=2, width=0.6, labelsize=8)
fig.savefig('Hela_pooled_exptest_mass.pdf', bbox_inches='tight')  
plt.close(fig) 

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
ax.plot(binned_dervol_vs_vol[0], binned_dervol_vs_vol[1], 'o', color=cool_blue, markeredgecolor='black', markeredgewidth='0.3', markersize=8,  alpha=1)     
#ax.errorbar(binned_dervol_vs_vol[0], binned_dervol_vs_vol[1], yerr=binned_dervol_vs_vol[5],
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1)    
ax.set_xlabel('volume'r'($\mathrm{\mu m}^3$)', fontsize=10)
ax.set_ylabel('mean derivative \n'r'$\langle \frac{dV}{dt} \rangle$ ($\mathrm{\mu m}^3$ $\mathrm{h}^{-1}$)', fontsize=10)
ax.set_ylim(40, 130)
ax.set_yticks([50, 75, 100, 125])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=8)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=8)
fig.savefig('Hela_pooled_exptest_vol.pdf', bbox_inches='tight')  
plt.close(fig) 