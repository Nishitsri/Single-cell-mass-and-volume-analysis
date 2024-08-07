##to measure target density in the dataset

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

#plt.style.use('mystyle')

############FUNCTION###############

def spring_function(x, k, x0):
    return np.multiply(np.subtract(x, x0), -k)

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
parser.add_argument("--tau", type=float, default=1, help="Derivative time step (h)")
parser.add_argument("--window_size", type=int, default=150, help="Number of points within a window of the binned scatter plot for correlations")

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

# DENSITY TRACKS
cells_density_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_density, spreadingphase_times[1], mitotic_times[1]))

# DERIVATIVE TIME STEP
tau = 2

# DENSITY DERIVATIVES
cells_density_der_as_bm = cc.der_linear_fit(cells_density_as_bm, tau)

# DENSITY SCATTER/BINNED PLOTS

## UNFILTERED DATA
binning_window_size = 1000 ## number of points windows a single windows size

density_scatter_as_bm, corr_density_scatter_as_bm = cc.scatter_plot_vectors(cells_density_as_bm, cells_density_der_as_bm)
slope, intercept, r_value, p_value, std_err = linregress(density_scatter_as_bm[0], density_scatter_as_bm[1])

density_scatter_as_bm_binned_mean = cc.binningdata_bins_mean(density_scatter_as_bm[0], density_scatter_as_bm[1], binning_window_size)
density_scatter_as_bm_binned_median = cc.binningdata_bins_median(density_scatter_as_bm[0], density_scatter_as_bm[1], binning_window_size)
density_scatter_as_bm_par, _ = curve_fit(spring_function, density_scatter_as_bm_binned_median[0], density_scatter_as_bm_binned_median[1], p0=[0.15, 0.1])
density_scatter_as_bm_binned_median_fit = [
density_scatter_as_bm_binned_median[0], spring_function(density_scatter_as_bm_binned_median[0], density_scatter_as_bm_par[0], density_scatter_as_bm_par[1])]

density_scatter_as_bm_par_mean, _ = curve_fit(spring_function, density_scatter_as_bm_binned_mean[0], density_scatter_as_bm_binned_mean[1], p0=[0.15, 0.1])


## FILTERED DATA
## filtering data may be useful to remove outliers - we try out the same plot filtering out a certain percentage of the distribution 
decile_filter = 1
filter_binning_window_size = 1000 ## number of points windows a single windows size

### SCATER PLOT
filter_density_scatter_as_bm = filter_data(density_scatter_as_bm, decile_filter)
filter_slope, filter_intercept, filter_r_value, filter_p_value, filter_std_err = linregress(filter_density_scatter_as_bm[0], filter_density_scatter_as_bm[1])
### BINNED PLOT
filter_density_scatter_as_bm_binned_mean = cc.binningdata_bins_mean(filter_density_scatter_as_bm[0], filter_density_scatter_as_bm[1], filter_binning_window_size)
filter_density_scatter_as_bm_binned_median = cc.binningdata_bins_median(filter_density_scatter_as_bm[0], filter_density_scatter_as_bm[1], filter_binning_window_size)
filter_density_scatter_as_bm_par, _ = curve_fit(spring_function, filter_density_scatter_as_bm_binned_median[0], filter_density_scatter_as_bm_binned_median[1], p0=[0.15, 0.1])
filter_density_scatter_as_bm_binned_median_fit = [
filter_density_scatter_as_bm_binned_median[0], spring_function(filter_density_scatter_as_bm_binned_median[0], filter_density_scatter_as_bm_par[0], filter_density_scatter_as_bm_par[1])
]
filter_density_scatter_as_bm_par_mean, _ = curve_fit(spring_function, filter_density_scatter_as_bm_binned_mean[0], filter_density_scatter_as_bm_binned_mean[1], p0=[0.15, 0.1])
###############PRINTING####################
print('target density from fit = ', filter_density_scatter_as_bm_par_mean[1])
print('slope from fit = ', filter_density_scatter_as_bm_par_mean[0])
print('target density from fit with filter = ', density_scatter_as_bm_par_mean[1])
print('slope from fit with filter = ', density_scatter_as_bm_par_mean[0])
df = pd.DataFrame({'target density from fit = ': filter_density_scatter_as_bm_par_mean[1],'slope from fit = ': filter_density_scatter_as_bm_par_mean[0],'target density from fit with filter = ': density_scatter_as_bm_par_mean[1],'slope from fit with filter = ': density_scatter_as_bm_par_mean[0], 'r2': r_value, 'p value': p_value} , index = [0])
df.to_csv('stats_helaadherent_tau%.2fh.csv' % tau)
###############PLOTTING####################

# COLORS
cool_red = '#a50f15'
cool_blue = '#253494'
cool_purple = '#810f7c'

# TARGET DENSITY DENSITY PLOT UNFILTERD
width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.7)
    ax.spines[['right', 'top']].set_visible(False)
xy = np.vstack([density_scatter_as_bm[0], density_scatter_as_bm[1]])
z = gaussian_kde(xy)(xy)    
ax.errorbar(density_scatter_as_bm_binned_mean[0], density_scatter_as_bm_binned_mean[1], xerr=density_scatter_as_bm_binned_mean[4], yerr=density_scatter_as_bm_binned_mean[5],
    marker='o', color='whitesmoke', capsize=4, linestyle='', markeredgecolor='black', markeredgewidth='0.3', markersize=6,  alpha=1, label='binned avg')
#ax.plot(density_scatter_as_bm_binned_mean[0], density_scatter_as_bm_binned_mean[1], 'o',
    #color='gainsboro', markeredgecolor='black', markeredgewidth='1', markersize=10, alpha=1, label='binned average')
sns.scatterplot(x=density_scatter_as_bm[0], y=density_scatter_as_bm[1], 
    hue=z, palette="magma", s=5, edgecolor="black", linewidth=0, legend=False, ax=ax)
ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=7)
ax.set_ylabel('mean density derivative\n' r'$\langle \frac{d\rho}{dt}\rangle$ (pg $\mathrm{\mu m}^{-3}hr^{-1})$', fontsize=7)
#ax.set_xticks([0.05, 0.1, 0.15, 0.2, 0.25])
#ax.set_yticks([-0.03, 0, 0.03, 0.06, 0.09])
#ax.set_xlim(-0.03, 0.09)
#ax.set_ylim(-0.025, 0.025)
ax.tick_params(axis='x', length=2, width=0.7, labelsize=6)
ax.tick_params(axis='y', length=2, width=0.7, labelsize=6)
#ax.legend(fontsize=7, loc='lower left')
fig.savefig('Helaadherent_target_density_plot_scatter_tau%.2f.pdf'%tau, bbox_inches='tight' )
plt.close(fig)

# TARGET DENSITY BINNED PLOT UNFILTERED
width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.7)
    ax.spines[['right', 'top']].set_visible(False)
ax.axvline(x=density_scatter_as_bm_par_mean[1], linestyle='--', linewidth=2, color='dimgray', label='target density')
ax.plot([0, 1], spring_function([0, 1], density_scatter_as_bm_par_mean[0], density_scatter_as_bm_par_mean[1]), '-', color='black', linewidth=2, label='lin fit')
ax.errorbar(density_scatter_as_bm_binned_mean[0], density_scatter_as_bm_binned_mean[1], xerr=density_scatter_as_bm_binned_mean[4], yerr=density_scatter_as_bm_binned_mean[5],
    marker='o', color=cool_purple, capsize=4, linestyle='', markeredgecolor='black', markeredgewidth='0.3', markersize=6,  alpha=1, label='binned avg')
#ax.plot(density_scatter_as_bm_binned_mean[0], density_scatter_as_bm_binned_mean[1], 'o',
    #color=cool_purple, markeredgecolor='black', markeredgewidth='1', markersize=10, alpha=1, label=f'p = {p_value:.2e}\nr = {r_value:.2f}')
ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=7)
ax.set_ylabel('mean density derivative\n' r'$\langle \frac{d\rho}{dt}\rangle$ (pg $\mathrm{\mu m}^{-3}hr^{-1})$', fontsize=7)
ax.set_xlim(0.1, 0.2)
ax.set_xticks([0.1, 0.15, 0.2])
ax.set_ylim(-0.005, 0.005)
#ax.set_yticks([-0.002, -0.001, 0, 0.001, 0.002])
ax.tick_params(axis='x', length=2, width=0.7, labelsize=6)
ax.tick_params(axis='y', length=2, width=0.7, labelsize=6)
#ax.legend(fontsize=5)
fig.savefig('Helaadherent_target_density_plot_nofilter_binned_tau%.2f.pdf'%tau, bbox_inches='tight')
plt.close(fig)

# TARGET DENSITY DENSITY PLOT FILTERED
width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.7)
    ax.spines[['right', 'top']].set_visible(False)
xy = np.vstack([filter_density_scatter_as_bm[0], filter_density_scatter_as_bm[1]])
z = gaussian_kde(xy)(xy)    
ax.errorbar(filter_density_scatter_as_bm_binned_mean[0], filter_density_scatter_as_bm_binned_mean[1], xerr=filter_density_scatter_as_bm_binned_mean[4], yerr=filter_density_scatter_as_bm_binned_mean[5],
    marker='o', color='whitesmoke', capsize=4, linestyle='', markeredgecolor='black', markeredgewidth='0.3', markersize=6,  alpha=1, label='binned avg')
#ax.plot(filter_density_scatter_as_bm_binned_mean[0], filter_density_scatter_as_bm_binned_mean[1], color='gainsboro',linestyle='None',marker='o', markeredgecolor='black', markeredgewidth='1', markersize=10, label='binned average')
sns.scatterplot(x=filter_density_scatter_as_bm[0], y=filter_density_scatter_as_bm[1],hue=z, palette="magma", s=5, edgecolor="black", linewidth=0, legend=False, ax=ax)
ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=7)
ax.set_ylabel('mean density derivative\n' r'$\langle \frac{d\rho}{dt}\rangle$ (pg $\mathrm{\mu m}^{-3}hr^{-1})$', fontsize=7)
ax.set_xticks([0.1, 0.15, 0.2])
#ax.set_yticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.set_xlim(0.1, 0.22)
ax.set_ylim(-0.025, 0.025)  #use for 1 and 2 hours
#ax.set_ylim(-0.005, 0.005)
ax.tick_params(axis='x', length=2, width=0.7, labelsize=6)
ax.tick_params(axis='y', length=2, width=0.7, labelsize=6)
#ax.legend(fontsize=5, loc='lower left')
fig.savefig('Helaadherent_target_density_plot_filter%gpercent_tau%.2fh.pdf' % (decile_filter, tau), bbox_inches='tight')
plt.close(fig)

# TARGET DENSITY BINNED PLOT FILTERED
width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.7)
    ax.spines[['right', 'top']].set_visible(False)
#ax.errorbar(scatter_as_bm_binned_mean[0], scatter_as_bm_binned_mean[1], yerr=scatter_as_bm_binned_mean[5], 
#    marker='o',color='black', capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=5,  alpha=1)
ax.axvline(x=filter_density_scatter_as_bm_par_mean[1], linestyle='--', linewidth=0.7, color='dimgray', label='target density')
ax.plot([0, 1], spring_function([0, 1], filter_density_scatter_as_bm_par_mean[0], filter_density_scatter_as_bm_par_mean[1]), '-', color='black', linewidth=0.7, label='lin fit')
ax.errorbar(filter_density_scatter_as_bm_binned_mean[0], filter_density_scatter_as_bm_binned_mean[1], xerr= filter_density_scatter_as_bm_binned_mean[4],yerr= filter_density_scatter_as_bm_binned_mean[5], color=cool_purple,linestyle='None',marker='o', markeredgecolor='black', markeredgewidth='0.3',capsize = 4, markersize=6, label='binned avg')
ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=7)
ax.set_ylabel('mean density derivative\n' r'$\langle \frac{d\rho}{dt}\rangle$ (pg $\mathrm{\mu m}^{-3}hr^{-1})$', fontsize=7)
ax.set_xticks([0.1, 0.15, 0.2])
#ax.set_yticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.set_xlim(0.10, 0.20)
ax.set_ylim(-0.002, 0.002) #use for tau 1 and 2 hours
#ax.set_ylim(-0.005, 0.005)
ax.tick_params(axis='x', length=2, width=0.7, labelsize=6)
ax.tick_params(axis='y', length=2, width=0.7, labelsize=6)
#ax.legend(fontsize=7)
fig.savefig('Helaadherent_target_density_binned_plot_filter_%g_tau%.2f.pdf' % (decile_filter,tau),  bbox_inches = 'tight' )
plt.close(fig)


# TARGET DENSITY DENSITY PLOT FILTERED
width = 2
height = width/1.5
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.7)
    ax.spines[['right', 'top']].set_visible(False)
xy = np.vstack([filter_density_scatter_as_bm[0], filter_density_scatter_as_bm[1]])
z = gaussian_kde(xy)(xy)
#ax.axvline(x=filter_density_scatter_as_bm_par_mean[1], linestyle='--', linewidth=2, color='dimgray', label='target density')
ax.plot([0, 1], spring_function([0, 1], filter_density_scatter_as_bm_par_mean[0], filter_density_scatter_as_bm_par_mean[1]), '-', color='black', linewidth=1.5, label='lin fit')
ax.errorbar(filter_density_scatter_as_bm_binned_mean[0], filter_density_scatter_as_bm_binned_mean[1], xerr=filter_density_scatter_as_bm_binned_mean[4], yerr=filter_density_scatter_as_bm_binned_mean[5],
    marker='o', color=cool_purple, capsize=4, linestyle='', markeredgecolor='black', markeredgewidth='0.3', markersize=6,  alpha=1, label='binned avg')
#ax.plot(filter_density_scatter_as_bm_binned_mean[0], filter_density_scatter_as_bm_binned_mean[1], 'o',
#    color=cool_purple, markeredgecolor='black', markeredgewidth='1', markersize=10, alpha=1, label='binned avg')
sns.scatterplot(x=filter_density_scatter_as_bm[0], y=filter_density_scatter_as_bm[1],
    hue=z, palette="magma", alpha=0.5, s=3, edgecolor="black", linewidth=0, legend=False, ax=ax)
#ax.text(0.092, -0.006, 'r = %.2f\np = %.2e' % (filter_r_value, filter_p_value), color='gray', fontsize=5)
ax.set_xlabel(r'density $\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=7)
ax.set_ylabel('mean density \n derivative\n' r'$\langle \frac{d\rho}{dt}\rangle$' '\n (pg' r'$\mathrm{\mu m}^{-3}hr^{-1})$', fontsize=7)
ax.set_xlim(0.10, 0.25)
ax.set_xticks([0.1, 0.15, 0.2])
ax.set_ylim(-0.002, 0.002)
#ax.set_yticks([-0.008, -0.004, 0, 0.004, 0.008])
#ax.set_yticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.tick_params(axis='x', length=2, width=0.7, labelsize=6)
ax.tick_params(axis='y', length=2, width=0.7, labelsize=6)
#ax.legend(fontsize=3, loc='upper right')
fig.savefig('Helaadherent_target_density_plot_filter%gpercent_tau%.2fh_allinone.pdf' % (decile_filter, tau), bbox_inches='tight')
plt.close(fig)
