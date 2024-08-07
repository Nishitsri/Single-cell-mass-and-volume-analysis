##measuring correlations between different variables
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

###########FUNCTIONS################

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

# EXPONENTIAL TEST
cells_mass_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_mass, spreadingphase_times[1], mitotic_times[1]))
cells_vol_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_vol, spreadingphase_times[1], mitotic_times[1]))

# DERIVATIVE TIME STEP
tau = 1

# BINNING WINDIW SIZE
window_size = 800
filter_window_size = 500

# GROWTH RATES
cells_grate_mass_as_bm = cc.growth_rate_linear_fit(cells_mass_as_bm, tau)
cells_grate_vol_as_bm = cc.growth_rate_linear_fit(cells_vol_as_bm, tau)

gratemass_vs_gratevol, corr_gratemass_vs_gratevol = cc.scatter_plot_vectors(cells_grate_mass_as_bm, cells_grate_vol_as_bm)
binned_gratemass_X_vs_gratevol_Y = cc.binningdata_bins_mean(gratemass_vs_gratevol[0], gratemass_vs_gratevol[1], window_size)
binned_gratevol_Y_vs_gratemass_X = cc.binningdata_bins_mean(gratemass_vs_gratevol[1], gratemass_vs_gratevol[0], window_size)

decile_filter = 5
filter_gratemass_vs_gratevol = filter_data(gratemass_vs_gratevol, decile_filter)
filter_binned_gratemass_X_vs_gratevol_Y = cc.binningdata_bins_mean(filter_gratemass_vs_gratevol[0], filter_gratemass_vs_gratevol[1], filter_window_size)
filter_binned_gratevol_Y_vs_gratemass_X = cc.binningdata_bins_mean(filter_gratemass_vs_gratevol[1], filter_gratemass_vs_gratevol[0], filter_window_size)

######################PRINTING#######################

######################PLOTTING#######################

# COLORS
cool_red = '#a50f15'
cool_blue = '#253494'
cool_purple = '#810f7c'

# GROWTH RATE CORRELATION WITH DENSITY PLOT
width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
xy = np.vstack([gratemass_vs_gratevol[0], gratemass_vs_gratevol[1]])
z = gaussian_kde(xy)(xy)    
#ax.errorbar(binned_gratemass_X_vs_gratevol_Y[0], binned_gratemass_X_vs_gratevol_Y[1], yerr=binned_gratemass_X_vs_gratevol_Y[5], 
#    marker='o',color='black', capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=5,  alpha=1) 
ax.plot(binned_gratemass_X_vs_gratevol_Y[0], binned_gratemass_X_vs_gratevol_Y[1], 'o',
    color='gainsboro', markeredgecolor='black', markeredgewidth='1', markersize=8,  alpha=1)
sns.scatterplot(x=gratemass_vs_gratevol[0], y=gratemass_vs_gratevol[1], 
    hue=z, palette="magma_r", s=3, edgecolor="black", linewidth=0, legend=False, ax=ax)
ax.set_xlabel(r'$\lambda_M$ ''\nmass growth rate', fontsize=10)
ax.set_ylabel(r'$\lambda_V$''\nvolume growth rate', fontsize=10)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
fig.savefig('Hela_growth_rate_cross_correlation_density_plot_no_filter.pdf', bbox_inches='tight')  
plt.close(fig) 

width = 2
height = 2#width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
xy = np.vstack([filter_gratemass_vs_gratevol[0], filter_gratemass_vs_gratevol[1]])
z = gaussian_kde(xy)(xy)    
#ax.errorbar(filter_binned_gratemass_X_vs_gratevol_Y[0], filter_binned_gratemass_X_vs_gratevol_Y[1], yerr=filter_binned_gratemass_X_vs_gratevol_Y[5], 
#    marker='o',color='black', capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=5,  alpha=1) 
ax.plot(filter_binned_gratemass_X_vs_gratevol_Y[0], filter_binned_gratemass_X_vs_gratevol_Y[1], 'o',
    color='gainsboro', markeredgecolor='black', markeredgewidth='1', markersize=8, alpha=1, label='binned average')
sns.scatterplot(x=filter_gratemass_vs_gratevol[0], y=filter_gratemass_vs_gratevol[1], 
    hue=z, palette="magma_r", s=5, edgecolor="black", linewidth=0, legend=False, ax=ax)
ax.set_xlabel(r'$\lambda_M$ ''\nmass growth rate', fontsize=10)
ax.set_ylabel(r'$\lambda_V$''\nvolume growth rate', fontsize=10)
ax.set_xticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.set_yticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.set_xlim(-0.03, 0.09)
ax.set_ylim(-0.03, 0.09)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
ax.legend(fontsize=7, loc='lower right')
fig.savefig('Hela_growth_rate_cross_correlation_density_plot_filter_%gpercent.pdf' % decile_filter, bbox_inches='tight')  
plt.close(fig) 

# GROWTH RATE DISTRIBUTIONS
width = 1.5
height = 0.33
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99) 
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
sns.histplot(x=filter_gratemass_vs_gratevol[0], stat="density", kde=True, color=cool_red, ax=ax)
ax.set_xlabel(r'$\lambda_M$ ''\nmass growth rate', fontsize=10)
ax.set_ylabel(r'density', fontsize=10)
ax.set_xlim(-0.03, 0.09)
ax.set_xticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
fig.savefig('Hela_growth_rate_cross_correlation_growth_rate_mass_histfilter_%gpercent.pdf' % decile_filter, bbox_inches='tight')  
plt.close(fig) 

width = 2
height = 0.33
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99) 
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
sns.histplot(x=filter_gratemass_vs_gratevol[1], stat="density", kde=True, color=cool_blue, ax=ax)
ax.set_xlabel(r'$\lambda_V$ ''\nvolume growth rate', fontsize=10)
ax.set_ylabel(r'density', fontsize=10)
ax.set_xlim(-0.03, 0.09)
ax.set_xticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
fig.savefig('Hela_growth_rate_cross_correlation_growth_rate_volume_histfilter_%gpercent.pdf' % decile_filter, bbox_inches='tight')  
plt.close(fig) 

# GROWTH RATE CORRELATION BINNED
width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot([-0.04, 0.09], [-0.04, 0.09], '--', color='black', linewidth=2, label='x=y')   
ax.plot(binned_gratemass_X_vs_gratevol_Y[0], binned_gratemass_X_vs_gratevol_Y[1], 'o',
    color='black', markeredgecolor='gray', markeredgewidth='1', markersize=10, alpha=1) 
#ax.errorbar(binned_gratemass_X_vs_gratevol_Y[0], binned_gratemass_X_vs_gratevol_Y[1], yerr=binned_gratemass_X_vs_gratevol_Y[5], 
#    marker='o',color='black', capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=5,  alpha=1) 
ax.set_xlabel(r'$\lambda_M$ ''\nmass growth rate', fontsize=10)
ax.set_ylabel(r'$\lambda_V$''\nvolume growth rate', fontsize=10)
ax.set_xticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.set_yticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.set_xlim(-0.04, 0.09)
ax.set_ylim(-0.04, 0.09)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
ax.legend(fontsize=7)
fig.savefig('Hela_growth_rate_cross_correlation_binned_plot_massX_volY.pdf', bbox_inches='tight')  
plt.close(fig) 

width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot([-0.04, 0.09], [-0.04, 0.09], '--', color='black', linewidth=2, label='x=y')   
ax.plot(binned_gratevol_Y_vs_gratemass_X[0], binned_gratevol_Y_vs_gratemass_X[1], 'o',
    color='black', markeredgecolor='gray', markeredgewidth='1', markersize=10, alpha=1) 
#ax.errorbar(binned_gratemass_X_vs_gratevol_Y[0], binned_gratemass_X_vs_gratevol_Y[1], yerr=binned_gratemass_X_vs_gratevol_Y[5], 
#    marker='o',color='black', capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=5,  alpha=1) 
ax.set_xlabel(r'$\lambda_V$ ''\nvolume growth rate', fontsize=10)
ax.set_ylabel(r'$\lambda_M$''\nmass growth rate', fontsize=10)
ax.set_xticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.set_yticks([-0.03, 0, 0.03, 0.06, 0.09])
ax.set_xlim(-0.04, 0.09)
ax.set_ylim(-0.04, 0.09)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
ax.legend(fontsize=7)
fig.savefig('Hela_growth_rate_cross_correlation_binned_plot_volX_massY.pdf', bbox_inches='tight')  
plt.close(fig) 

