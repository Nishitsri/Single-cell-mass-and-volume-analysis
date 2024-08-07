##To measure mass and volume growth derivatives

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

# DERIVATIVE TIME STEP
tau = 1

# COMPUTE DERIVATIVES
cells_der_mass = cc.der_linear_fit(cells_mass, tau)
cells_der_vol = cc.der_linear_fit(cells_vol, tau)

# COMPUTE GROWTH RATE
cells_grate_mass = cc.growth_rate_linear_fit(cells_mass, tau)
cells_grate_vol = cc.growth_rate_linear_fit(cells_vol, tau)

######################STATISTICS#######################

# AVG BINNING WINDOW SIZE
dt_window = 1
dphase_window = 0.01

# MASS GROWTH RATE
mean_grate_mass_at_fixed_time, sterr_grate_mass_at_fixed_time, stdev_grate_mass_at_fixed_time, cv_grate_mass_at_fixed_time, time_grate_mass = mean_stdev_cv_vs_time(cells_grate_mass, dt_window)
mean_grate_mass_at_fixed_phase, sterr_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase, cv_grate_mass_at_fixed_phase, phase_grate_mass = mean_stdev_cv_vs_phase(cells_grate_mass, dphase_window)

# VOL GROWTH RATE
mean_grate_vol_at_fixed_time, sterr_grate_vol_at_fixed_time, stdev_grate_vol_at_fixed_time, cv_grate_vol_at_fixed_time, time_grate_vol = mean_stdev_cv_vs_time(cells_grate_vol, dt_window)
mean_grate_vol_at_fixed_phase, sterr_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase, cv_grate_vol_at_fixed_phase, phase_grate_vol = mean_stdev_cv_vs_phase(cells_grate_vol, dphase_window)

# CONDITIONED MEAN BINNING WINDOW SIZE
window_size = 1200

# EXPONENTIAL TEST
cells_mass_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_mass, spreadingphase_times[1], mitotic_times[1]))
cells_vol_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_vol, spreadingphase_times[1], mitotic_times[1]))

cells_der_mass_as_bm = cc.der_linear_fit(cells_mass_as_bm, tau)
cells_der_vol_as_bm = cc.der_linear_fit(cells_vol_as_bm, tau)

dermass_vs_mass, corr_dermass_vs_mass = cc.scatter_plot_vectors(cells_mass_as_bm, cells_der_mass_as_bm)
binned_dermass_vs_mass = cc.binningdata_bins_mean(dermass_vs_mass[0], dermass_vs_mass[1], window_size)

dervol_vs_vol, corr_dervol_vs_vol = cc.scatter_plot_vectors(cells_vol_as_bm, cells_der_vol_as_bm)
binned_dervol_vs_vol = cc.binningdata_bins_mean(dervol_vs_vol[0], dervol_vs_vol[1], window_size)

######################PRINTING#######################

######################PLOTTING#######################

# COLORS

cool_red = '#a50f15'
cool_blue = '#253494'
cool_purple = '#810f7c'

# EXPONENTIAL TEST

width = 1.7
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
#ax.plot(binned_dermass_vs_mass[0], binned_dermass_vs_mass[1], 'o', color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=8,  alpha=1)
ax.errorbar(binned_dermass_vs_mass[0], binned_dermass_vs_mass[1], yerr=binned_dermass_vs_mass[5],   marker='o', color=cool_red, capsize=4, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=8,  alpha=1)
ax.set_xlabel('Mass' r'$\langle M \rangle$ pg', fontsize=10)
ax.set_ylabel('mean mass \n derivative \n'r'$\langle \frac{dM}{dt} \rangle$ (pg $\mathrm{h}^{-1}$)', fontsize=10)
ax.set_ylim(5, 15)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
fig.savefig('Hela_pooled_exptest_mass.pdf', bbox_inches='tight')  
plt.close(fig) 

width = 1.7
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
#ax.plot(binned_dervol_vs_vol[0], binned_dervol_vs_vol[1], 'o', color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=8,  alpha=1)
ax.errorbar(binned_dervol_vs_vol[0], binned_dervol_vs_vol[1], yerr=binned_dervol_vs_vol[5],marker='o', color=cool_blue, capsize=4, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=8,  alpha=1)
ax.set_xlabel('Volume'r'$\langle V \rangle$ ($\mathrm{\mu m}^3$)', fontsize=10)
ax.set_ylabel('mean volume \n derivative \n'r'$\langle \frac{dV}{dt} \rangle$ ($\mathrm{\mu m}^3$ $\mathrm{h}^{-1}$)', fontsize=10)
ax.set_ylim(40, 130)
ax.set_yticks([50, 80, 110])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
fig.savefig('Hela_pooled_exptest_vol.pdf', bbox_inches='tight')  
plt.close(fig) 

# GROWTH RATE VS CELL CYCLE TIME

width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_grate_mass, mean_grate_mass_at_fixed_phase, 'o', 
    color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=8,  alpha=1, label='mass growth rate')
ax.fill_between(phase_grate_mass, np.subtract(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), 
                np.add(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), color='red', alpha=0.2) 
ax.plot(phase_grate_vol, mean_grate_vol_at_fixed_phase, 'o', 
    color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=8, alpha=1, label='volume growth rate')
ax.fill_between(phase_grate_vol, np.subtract(mean_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase), 
                np.add(mean_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase), color='blue', alpha=0.2) 
#ax.errorbar(phase_grate_mass, mean_grate_mass_at_fixed_phase, yerr=sterr_grate_mass_at_fixed_phase,
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass') 
#ax.errorbar(phase_grate_vol, mean_grate_vol_at_fixed_phase, yerr=sterr_grate_vol_at_fixed_phase,
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='volume') 
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean growth \n rate 'r'$\langle \lambda \rangle$ ($\mathrm{h}^{-1}$)', fontsize=10)
#ax.set_yticks([-0.05, 0, 0.05, 0.1, 0.15, 0.2])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
ax.legend(fontsize=4)
fig.savefig('Hela_pooled_growth_rate_vs_cell_cycle_phase.pdf', bbox_inches='tight')  
plt.close(fig) 

#without legend

width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_grate_mass, mean_grate_mass_at_fixed_phase, 'o',
    color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=8,  alpha=1, label='mass growth rate')
ax.fill_between(phase_grate_mass, np.subtract(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase),
                np.add(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), color='red', alpha=0.2)
ax.plot(phase_grate_vol, mean_grate_vol_at_fixed_phase, 'o',
    color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=8, alpha=1, label='volume growth rate')
ax.fill_between(phase_grate_vol, np.subtract(mean_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase),
                np.add(mean_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase), color='blue', alpha=0.2)
#ax.errorbar(phase_grate_mass, mean_grate_mass_at_fixed_phase, yerr=sterr_grate_mass_at_fixed_phase,
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass')
#ax.errorbar(phase_grate_vol, mean_grate_vol_at_fixed_phase, yerr=sterr_grate_vol_at_fixed_phase,
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='volume')
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean growth \n rate 'r'$\langle \lambda \rangle$ ($\mathrm{h}^{-1}$)', fontsize=10)
#ax.set_yticks([-0.05, 0, 0.05, 0.1, 0.15, 0.2])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
#ax.legend(fontsize=4)
fig.savefig('Hela_pooled_growth_rate_vs_cell_cycle_phase.pdf', bbox_inches='tight')
plt.close(fig)



width = 6.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_grate_mass, mean_grate_mass_at_fixed_phase, 'o', 
    color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=18,  alpha=1, label='mass growth rate')
ax.fill_between(phase_grate_mass, np.subtract(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), 
                np.add(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), color='red', alpha=0.2) 
ax.plot(phase_grate_vol, mean_grate_vol_at_fixed_phase, 'o', 
    color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=14, alpha=1, label='volume growth rate')
ax.fill_between(phase_grate_vol, np.subtract(mean_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase), 
                np.add(mean_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase), color='blue', alpha=0.2) 
#ax.errorbar(phase_grate_mass, mean_grate_mass_at_fixed_phase, yerr=sterr_grate_mass_at_fixed_phase,
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass') 
#ax.errorbar(phase_grate_vol, mean_grate_vol_at_fixed_phase, yerr=sterr_grate_vol_at_fixed_phase,
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='volume') 
ax.set_xlabel(r'cell cycle phase', fontsize=12)
ax.set_ylabel('mean growth \n rate 'r'$\langle \lambda \rangle$ ($\mathrm{h}^{-1}$)', fontsize=12)
#ax.set_yticks([-0.05, 0, 0.05, 0.1, 0.15, 0.2])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=12)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=12)
ax.legend(fontsize=4)
fig.savefig('Hela_pooled_growth_rate_vs_cell_cycle_phase_big.pdf', bbox_inches='tight')
plt.close(fig) 

#without legend


width = 6.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_grate_mass, mean_grate_mass_at_fixed_phase, 'o',
    color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=14,  alpha=1, label='mass growth rate')
ax.fill_between(phase_grate_mass, np.subtract(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase),
                np.add(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), color='red', alpha=0.2)
ax.plot(phase_grate_vol, mean_grate_vol_at_fixed_phase, 'o',
    color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=14, alpha=1, label='volume growth rate')
ax.fill_between(phase_grate_vol, np.subtract(mean_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase),
                np.add(mean_grate_vol_at_fixed_phase, stdev_grate_vol_at_fixed_phase), color='blue', alpha=0.2)
#ax.errorbar(phase_grate_mass, mean_grate_mass_at_fixed_phase, yerr=sterr_grate_mass_at_fixed_phase,
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass')
#ax.errorbar(phase_grate_vol, mean_grate_vol_at_fixed_phase, yerr=sterr_grate_vol_at_fixed_phase,
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='volume')
ax.set_xlabel(r'cell cycle phase', fontsize=20)
ax.set_ylabel('mean growth \n rate 'r'$\langle \lambda \rangle$ ($\mathrm{h}^{-1}$)', fontsize=18)
ax.set_yticks([-0.05, 0.05, 0.15, 0.25])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=16)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=16)
#ax.legend(fontsize=4)
fig.savefig('Hela_pooled_growth_rate_vs_cell_cycle_phase_big_withoutlegend.pdf',bbox_inches='tight')
plt.close(fig)


width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_grate_mass, mean_grate_mass_at_fixed_phase, 'o', color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass') 
ax.fill_between(phase_grate_mass, np.subtract(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), 
                np.add(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), color='red', alpha=0.2)
ax.plot(phase_grate_vol[::2], mean_grate_vol_at_fixed_phase[::2], 'o', color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=8, alpha=1, label='volume')
ax.fill_between(phase_grate_vol[::2], np.subtract(mean_grate_vol_at_fixed_phase[::2], stdev_grate_vol_at_fixed_phase[::2]), 
                np.add(mean_grate_vol_at_fixed_phase[::2], stdev_grate_vol_at_fixed_phase[::2]), color='blue', alpha=0.2) 
#ax.errorbar(phase_grate_mass, mean_grate_mass_at_fixed_phase, yerr=sterr_grate_mass_at_fixed_phase,
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass') 
#ax.errorbar(phase_grate_vol, mean_grate_vol_at_fixed_phase, yerr=sterr_grate_vol_at_fixed_phase,
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='volume') 
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean growth \n rate 'r'$\langle \lambda \rangle$ ($\mathrm{h}^{-1}$)', fontsize=10)
#ax.set_yticks([-0.05, 0, 0.05, 0.1, 0.15, 0.2])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
ax.legend(fontsize=4)
fig.savefig('Hela_pooled_growth_rate_vs_cell_cycle_phase_vol_sliced.pdf', bbox_inches='tight')  
plt.close(fig) 


##without legend
width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_grate_mass, mean_grate_mass_at_fixed_phase, 'o', color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass')
ax.fill_between(phase_grate_mass, np.subtract(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase),
                np.add(mean_grate_mass_at_fixed_phase, stdev_grate_mass_at_fixed_phase), color='red', alpha=0.2)
ax.plot(phase_grate_vol[::2], mean_grate_vol_at_fixed_phase[::2], 'o', color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=8, alpha=1, label='volume')
ax.fill_between(phase_grate_vol[::2], np.subtract(mean_grate_vol_at_fixed_phase[::2], stdev_grate_vol_at_fixed_phase[::2]),
                np.add(mean_grate_vol_at_fixed_phase[::2], stdev_grate_vol_at_fixed_phase[::2]), color='blue', alpha=0.2)
#ax.errorbar(phase_grate_mass, mean_grate_mass_at_fixed_phase, yerr=sterr_grate_mass_at_fixed_phase,
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass')
#ax.errorbar(phase_grate_vol, mean_grate_vol_at_fixed_phase, yerr=sterr_grate_vol_at_fixed_phase,
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='volume')
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean growth \n rate 'r'$\langle \lambda \rangle$ ($\mathrm{h}^{-1}$)', fontsize=10)
#ax.set_yticks([-0.05, 0, 0.05, 0.1, 0.15, 0.2])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
#ax.legend(fontsize=4)
fig.savefig('Hela_pooled_growth_rate_vs_cell_cycle_phase_vol_sliced_withoutlegend.pdf', bbox_inches='tight')
plt.close(fig)




width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_grate_mass[::2], mean_grate_mass_at_fixed_phase[::2], 'o', color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=10, alpha=1, label='mass') 
ax.fill_between(phase_grate_mass[::2], np.subtract(mean_grate_mass_at_fixed_phase[::2], stdev_grate_mass_at_fixed_phase[::2]), 
                np.add(mean_grate_mass_at_fixed_phase[::2], stdev_grate_mass_at_fixed_phase[::2]), color='red', alpha=0.2)
ax.plot(phase_grate_vol[::2], mean_grate_vol_at_fixed_phase[::2], 'o', color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=8, alpha=1, label='volume')
ax.fill_between(phase_grate_vol[::2], np.subtract(mean_grate_vol_at_fixed_phase[::2], stdev_grate_vol_at_fixed_phase[::2]), 
                np.add(mean_grate_vol_at_fixed_phase[::2], stdev_grate_vol_at_fixed_phase[::2]), color='blue', alpha=0.2) 
#ax.errorbar(phase_grate_mass[::2], mean_grate_mass_at_fixed_phase[::2], yerr=sterr_grate_mass_at_fixed_phase[::2],
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass') 
#ax.errorbar(phase_grate_vol[::2], mean_grate_vol_at_fixed_phase[::2], yerr=sterr_grate_vol_at_fixed_phase[::2],
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='volume')  
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean growth \n rate 'r'$\langle \lambda \rangle$ ($\mathrm{h}^{-1}$)', fontsize=10)
ax.set_yticks([-0.05, 0.05, 0.15])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
ax.legend(fontsize=2, loc= 'upper left' )
fig.savefig('Hela_pooled_growth_rate_vs_cell_cycle_phase_all_sliced.pdf', bbox_inches='tight')  
plt.close(fig) 


##without legend
width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=200)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_grate_mass[::2], mean_grate_mass_at_fixed_phase[::2], 'o', color=cool_red, markeredgecolor='black', markeredgewidth='1', markersize=10, alpha=1, label='mass')
ax.fill_between(phase_grate_mass[::2], np.subtract(mean_grate_mass_at_fixed_phase[::2], stdev_grate_mass_at_fixed_phase[::2]),
                np.add(mean_grate_mass_at_fixed_phase[::2], stdev_grate_mass_at_fixed_phase[::2]), color='red', alpha=0.2)
ax.plot(phase_grate_vol[::2], mean_grate_vol_at_fixed_phase[::2], 'o', color=cool_blue, markeredgecolor='black', markeredgewidth='1', markersize=8, alpha=1, label='volume')
ax.fill_between(phase_grate_vol[::2], np.subtract(mean_grate_vol_at_fixed_phase[::2], stdev_grate_vol_at_fixed_phase[::2]),
                np.add(mean_grate_vol_at_fixed_phase[::2], stdev_grate_vol_at_fixed_phase[::2]), color='blue', alpha=0.2)
#ax.errorbar(phase_grate_mass[::2], mean_grate_mass_at_fixed_phase[::2], yerr=sterr_grate_mass_at_fixed_phase[::2],
#    marker='o', color=cool_red, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='mass')
#ax.errorbar(phase_grate_vol[::2], mean_grate_vol_at_fixed_phase[::2], yerr=sterr_grate_vol_at_fixed_phase[::2],
#    marker='o', color=cool_blue, capsize=10, linestyle='', markeredgecolor='black', markeredgewidth='1', markersize=10,  alpha=1, label='volume')
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel('mean growth \n rate 'r'$\langle \lambda \rangle$ ($\mathrm{h}^{-1}$)', fontsize=10)
#ax.set_yticks([-0.05, 0.05, 0.15])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
#ax.legend(fontsize=2, loc= 'upper left' )
fig.savefig('Hela_pooled_growth_rate_vs_cell_cycle_phase_all_sliced_withoutlegend.pdf', bbox_inches='tight')
plt.close(fig)
