##Measure coefficient of variation of different variables

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
from scipy.stats.distributions import chi2
from scipy.stats import linregress
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from math import nan
from math import sqrt

#plt.style.use('mystyle')

############FUNCTION###############

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

def feltz_miller_test(cv_sample, size_sample):
    '''
    NB: the idea behind this test is that the dad variable as define below obeys chi-squared statistics with k - 1 degrees of freedom
    NB: k is the number of samples, remember that the only parameter in the chi-squared statistic is precisely k
    NB: then, the function return the pvalue of the dad variable
    '''

    sum_size_samples = sum([ size - 1 for size in size_sample ])

    est_true_cv = sum([ (size-1)*cv/sum_size_samples for size, cv in zip(size_sample, cv_sample) ])
    est_squared_cv = np.dot(np.subtract(size_sample, 1), np.multiply(cv_sample, cv_sample))/sum_size_samples

    dad = sum_size_samples/(est_true_cv*est_true_cv*(0.5 + est_true_cv*est_true_cv))*(est_squared_cv - est_true_cv*est_true_cv)

    return chi2.sf(dad, len(cv_sample)-1)

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
tau = args.tau ## tau has got to take the form k*0.25 because of the discrete nature of the experiment

# RELEVANT CELL CYCLE TIMES
division_times = cc.division_times(cells_mass)
mitotic_times = cc.mitotic_times(cells_vol, tau)
time_window, signa_delta_minus, signa_delta_plus = [0, 5], 1, 0
area_slope_change_times, area_slope_change_frame, delta_a_minus_at_slope_change, delta_a_plus_at_slope_change, der_ratio_at_slope_change = cc.observable_slope_change_times_ratio_method(cells_area, tau_slope_change,	
																																time_window, signa_delta_minus, signa_delta_plus)
spreadingphase_times = [area_slope_change_frame, area_slope_change_times]

######################TRACKS#######################

# END-OF-SPREADING SYNCH
cells_density_synchspread = cc.synch_cell_cycle_time_series_at_specific_timeframes(cells_density, spreadingphase_times[0])
cells_mass_synchspread = cc.synch_cell_cycle_time_series_at_specific_timeframes(cells_mass, spreadingphase_times[0])
cells_volume_synchspread = cc.synch_cell_cycle_time_series_at_specific_timeframes(cells_vol, spreadingphase_times[0])

# AFTER SPREADING BEFORE MITOSIS
cells_mass_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_mass, spreadingphase_times[1], mitotic_times[1]))
cells_vol_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_vol, spreadingphase_times[1], mitotic_times[1]))
cells_density_as_bm = cc.set_zero_t0(cc.cells_obs_between_timepoints(cells_density, spreadingphase_times[1], mitotic_times[1]))

###############CV ANALYSIS####################

# TARGET DENSITY
target_density = 0.14420501731519167

# DENSITY PDF AND CV AT BIRTH, END-OF-SPREADING AND MITOSIS 

## EOS
density_stat_around_eos = cc.value_in_time_range(cells_density, np.subtract(area_slope_change_times, 0.5), np.add(area_slope_change_times, 0.5))
cv_density_around_eos = stat.stdev(density_stat_around_eos)/stat.mean(density_stat_around_eos)

density_subtr_stat_around_eos = np.subtract(density_stat_around_eos, target_density)
mean_subtr_density_around_eos = stat.mean(density_subtr_stat_around_eos)
stdev_subtr_density_around_eos = stat.stdev(density_subtr_stat_around_eos)
sterr_subtr_density_around_eos = stdev_subtr_density_around_eos/sqrt(len(density_subtr_stat_around_eos))
cv_subtr_density_around_eos = stat.stdev(density_subtr_stat_around_eos)/stat.mean(density_subtr_stat_around_eos)

## MITOSIS
density_stat_around_mitosis = cc.value_in_time_range(cells_density, np.subtract(mitotic_times[1], 1), mitotic_times[1])
cv_density_around_mitosis = stat.stdev(density_stat_around_mitosis)/stat.mean(density_stat_around_mitosis)

density_subtr_stat_around_mitosis = np.subtract(density_stat_around_mitosis, target_density)
mean_subtr_density_around_mitosis = stat.mean(density_subtr_stat_around_mitosis)
stdev_subtr_density_around_mitosis = stat.stdev(density_subtr_stat_around_mitosis)
sterr_subtr_density_around_mitosis = stdev_subtr_density_around_mitosis/sqrt(len(density_subtr_stat_around_mitosis))
cv_subtr_density_around_mitosis = stat.stdev(density_subtr_stat_around_mitosis)/stat.mean(density_subtr_stat_around_mitosis)

## BIRTH
density_stat_around_birth = cc.value_in_time_range(cells_density, no_cells*[0], no_cells*[1])
cv_density_around_birth = stat.stdev(density_stat_around_birth)/stat.mean(density_stat_around_birth)

density_subtr_stat_around_birth = np.subtract(density_stat_around_birth, target_density)
mean_subtr_density_around_birth = stat.mean(density_subtr_stat_around_birth)
stdev_subtr_density_around_birth = stat.stdev(density_subtr_stat_around_birth)
sterr_subtr_density_around_birth = stdev_subtr_density_around_birth/sqrt(len(density_subtr_stat_around_birth))
cv_subtr_density_around_birth = stat.stdev(density_subtr_stat_around_birth)/stat.mean(density_subtr_stat_around_birth)

# FELTZ-MILLER TEST

feltz_miller_pvalue = feltz_miller_test([cv_density_around_birth, cv_density_around_eos, cv_density_around_mitosis], 
    [len(density_stat_around_birth), len(density_stat_around_eos), len(density_stat_around_mitosis)])
feltz_miller_pvalue_birth_eos = feltz_miller_test([cv_density_around_birth, cv_density_around_eos], [len(density_stat_around_birth), len(density_stat_around_eos)])
feltz_miller_pvalue_eos_mitosis = feltz_miller_test([cv_density_around_eos, cv_density_around_mitosis], [len(density_stat_around_eos), len(density_stat_around_mitosis)]) 
feltz_miller_pvalue_birth_mitosis = feltz_miller_test([cv_density_around_birth, cv_density_around_mitosis], [len(density_stat_around_birth), len(density_stat_around_birth)]) 

# REAL TIME
dt_binning = 1

# SYNCH TO BIRTH
mean_mass_at_fixed_time, sterr_mass_at_fixed_time, stdev_mass_at_fixed_time, cv_mass_at_fixed_time, time_mass = mean_stdev_cv_vs_time(cells_mass, dt_binning)
mean_volume_at_fixed_time, sterr_volume_at_fixed_time, stdev_volume_at_fixed_time, cv_volume_at_fixed_time, time_vol = mean_stdev_cv_vs_time(cells_vol, dt_binning)
mean_density_at_fixed_time, sterr_density_at_fixed_time, stdev_density_at_fixed_time, cv_density_at_fixed_time, time_density = mean_stdev_cv_vs_time(cells_density, dt_binning)

# SYNCH TO END-OF-SPREADING
#mean_mass_synchspread_at_fixed_time, sterr_mass_synchspread_at_fixed_time, stdev_mass_synchspread_at_fixed_time, cv_mass_synchspread_at_fixed_time, time_mass_synchspread = mean_stdev_cv_vs_time(cells_mass_synchspread, dt_binning)
#mean_volume_synchspread_at_fixed_time, sterr_volume_synchspread_at_fixed_time, stdev_volume_synchspread_at_fixed_time, cv_volume_synchspread_at_fixed_time, time_volume_synchspread = mean_stdev_cv_vs_time(cells_volume_synchspread, dt_binning)
#mean_density_synchspread_at_fixed_time, sterr_density_synchspread_at_fixed_time, stdev_density_synchspread_at_fixed_time, cv_density_synchspread_at_fixed_time, time_density_synchspread = mean_stdev_cv_vs_time(cells_density_synchspread, dt_binning)

# CELL CYCLE TIME
dphase_binning = 0.04

mean_mass_at_fixed_phase, sterr_mass_at_fixed_phase, stdev_mass_at_fixed_phase, cv_mass_at_fixed_phase, phase_mass = mean_stdev_cv_vs_phase(cells_mass, dphase_binning)
mean_volume_at_fixed_phase, sterr_volume_at_fixed_phase, stdev_volume_at_fixed_phase, cv_volume_at_fixed_phase, phase_vol = mean_stdev_cv_vs_phase(cells_vol, dphase_binning)
mean_density_at_fixed_phase, sterr_density_at_fixed_phase, stdev_density_at_fixed_phase, cv_density_at_fixed_phase, phase_density = mean_stdev_cv_vs_phase(cells_density, dphase_binning)

# INDEPENDENT NULL MODEL
cv_density_null_time = np.sqrt(np.add(
	np.multiply(cc.remove_nan(cv_mass_at_fixed_time), cc.remove_nan(cv_mass_at_fixed_time)), 
	np.multiply(cc.remove_nan(cv_volume_at_fixed_time), cc.remove_nan(cv_volume_at_fixed_time))
	)
	)
cv_density_null_phase = np.sqrt(np.add(
	np.multiply(cc.remove_nan(cv_mass_at_fixed_phase), cc.remove_nan(cv_mass_at_fixed_phase)), 
	np.multiply(cc.remove_nan(cv_volume_at_fixed_phase), cc.remove_nan(cv_volume_at_fixed_phase))
	)
	)

# DIFFERENCE WITH NULL MODEL
cv_diff_with_null_time = np.subtract(cv_density_null_time, cc.remove_nan(cv_density_at_fixed_time))
cv_diff_with_null_phase = np.subtract(cv_density_null_phase, cc.remove_nan(cv_density_at_fixed_phase))

###############PRINTING####################

print('birth, mean density minus target density = ', mean_subtr_density_around_birth)
print('birth, sterr density minus target density = ', sterr_subtr_density_around_birth)
print('birth, stdev density minus target density = ', stdev_subtr_density_around_birth)

print('eos, mean density minus target density = ', mean_subtr_density_around_eos)
print('eos, sterr density minus target density = ', sterr_subtr_density_around_eos)
print('eos, stdev density minus target density = ', stdev_subtr_density_around_eos)

print('mitosis, mean density minus target density = ', mean_subtr_density_around_mitosis)
print('mitosis, sterr density minus target density = ', sterr_subtr_density_around_mitosis)
print('mitosis, stdev density minus target density = ', stdev_subtr_density_around_mitosis)
df = pd.DataFrame({'birth, mean density minus target density = ': mean_subtr_density_around_birth, 'eos, mean density minus target density = ': mean_subtr_density_around_eos, 'mitosis, mean density minus target density = ': mean_subtr_density_around_mitosis }, index = [0])
df.to_csv('density_targetdensity_birth_EOS_mito.csv')
###############PLOTTING####################

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

# CV AT BIRTH, END-OF-SPREADING AND MITOSIS 
width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.bar(['Birth', 'EOS', 'Mitosis'], [stdev_subtr_density_around_birth, stdev_subtr_density_around_eos, stdev_subtr_density_around_mitosis], color='gray')
#ax.set_ylabel(r'$\sigma[\rho - \rho_{target}]$''\nstandard deviation', fontsize=10)
#ax.set_ylim(bottom=0.1)
ax.set_ylabel(r'$\sigma[\rho - \rho_{target}]$', fontsize=10)
ax.set_yticks([0, 0.01, 0.02, 0.03])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
fig.savefig('Hela_stdev_density_subtr_at_birth_eos_mitosis_tau%.2fh.pdf' % tau, bbox_inches='tight')   
plt.close(fig) 

# CV AT BIRTH, END-OF-SPREADING AND MITOSIS 
width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.bar(['Birth', 'EOS', 'Mitosis'], [cv_density_around_birth, cv_density_around_eos, cv_density_around_mitosis], color='gray')
#ax.text(2.7, 0.175, r'$p_{birth-eos}$ = %.2e' '\n' '$p_{eos-mitosis}$ = %.2e' '\n' r'$p_{birth-mitosis}$ = %.2e' % (feltz_miller_pvalue_birth_eos, feltz_miller_pvalue_eos_mitosis, feltz_miller_pvalue_birth_mitosis), color='gray', fontsize=7)
ax.set_ylabel(r'$CV[\rho]$''\ncoefficient of \n variation', fontsize=10)
ax.set_ylim(bottom=0.1)
#ax.set_yticks([1500, 2500, 3500, 4500])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
fig.savefig('Hela_CV_density_at_birth_eos_mitosis_tau%.2fh.pdf' % tau, bbox_inches='tight')   
plt.close(fig) 

# CV NULL VS DATA REAL TIME
width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(cv_density_null_time, cc.remove_nan(cv_density_at_fixed_time), 'o', color='purple', markeredgecolor='black', markeredgewidth='1', markersize=8)
ax.plot([0, 1], [0, 1], '--', color='black', linewidth=4, label='x=y')
ax.fill_between([0, 1], [0, 1], [1, 1], color='darkgrey', alpha=1)   
ax.fill_between([0, 1], [0, 0], [0, 1], color='whitesmoke', alpha=1)   
ax.set_xlabel(r'$\sqrt{CV[M]^2 + CV[V]^2 }$', fontsize=8)
ax.set_ylabel(r'$CV[\rho]$''\ncoefficient of \n variation', fontsize=8)
ax.set_xlim(0.15, 0.35)
ax.set_ylim(0.15, 0.35)
ax.set_xticks([0.15, 0.20, 0.25, 0.3, 0.35])
ax.set_yticks([0.15, 0.2, 0.25, 0.3, 0.35])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
#ax.legend(fontsize=7)
fig.savefig('Hela_CV_null_time_full_cell_cycle_tau%.2fh.pdf' % tau, bbox_inches='tight')  
plt.close(fig)

# CV NULL VS DATA CELL CYCLE TIME
width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(cv_density_null_phase, cc.remove_nan(cv_density_at_fixed_phase), 'o', color='purple', markeredgecolor='black', markeredgewidth='1', markersize=6)
ax.plot([0, 1], [0, 1], '--', color='black', linewidth=4, label='x=y')
ax.fill_between([0, 1], [0, 1], [1, 1], color='darkgrey', alpha=1)   
ax.fill_between([0, 1], [0, 0], [0, 1], color='whitesmoke', alpha=1)   
ax.set_xlabel(r'$\sqrt{CV[M]^2 + CV[V]^2 }$', fontsize=8)
ax.set_ylabel(r'$CV[\rho]$''\ncoefficient of \n variation', fontsize=8)
ax.set_xlim(0.15, 0.35)
ax.set_ylim(0.15, 0.35)
ax.set_xticks([0.15, 0.25, 0.35])
ax.set_yticks([0.15, 0.25, 0.35])
ax.tick_params(axis='x', length=4, width=1.5, labelsize=8)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=8)
#ax.legend(fontsize=7)
fig.savefig('Hela_CV_null_phase_full_cell_cycle_tau%.2fh.pdf' % tau, bbox_inches='tight')  
plt.close(fig)

# BOX PLOT CV NULL VS DATA BOX
width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.axhline(y=0, linestyle='--', linewidth=2, color='black', label='prediction\nno coupling')
boxplots = ax.boxplot([cv_diff_with_null_phase], patch_artist=True,
                whiskerprops=dict(color='black'), capprops=dict(color='black'), medianprops=dict(color='black'), 
                flierprops=dict(marker='o', markerfacecolor=cool_purple, markeredgecolor='gray', markeredgewidth=0.5, markersize=2, alpha=0.5, linestyle='none', label='data')) 
colors = [cool_purple]
for patch, color in zip(boxplots['boxes'], colors):
    patch.set_facecolor(color)
ax.set_ylabel(r'$CV_{null}[\rho] - CV[\rho]$''\ncoefficient of \n variation', fontsize=8)
ax.set_xticklabels(['HeLa \n adherent'], fontsize=8)
#ax.set_ylim(-0.1, 0.01)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
#ax.legend(fontsize=4)
fig.savefig('Hela_CV_null_vs_data_boxplot_phase_full_cell_cycle_tau%.2fh.pdf' % tau, bbox_inches='tight')  
plt.close(fig)

# CV REAL TIME
width = 2
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(time_density, cv_density_at_fixed_time, 'o', color='purple', markeredgecolor='black', markeredgewidth='1', markersize=8, label='density')
ax.plot(time_mass, cv_mass_at_fixed_time, 'o', color='red', markeredgecolor='black', markeredgewidth='1', markersize=8, label='mass')
ax.plot(time_vol, cv_volume_at_fixed_time, 'o', color='blue', markeredgecolor='black', markeredgewidth='1', markersize=8, label='volume')
ax.set_xlim(-1, 16)
ax.set_ylim(0.1, 0.3)
ax.set_yticks([0.1, 0.2, 0.3])
ax.set_xlabel('time from birth (h)', fontsize=10)
ax.set_ylabel(r'CV' '\ncoefficient of \n variation', fontsize=10)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=10)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=10)
#ax.legend(fontsize=5)
fig.savefig('Hela_CV_time_full_cell_cycle_tau%.2fh.pdf' % tau, bbox_inches='tight')  
plt.close(fig)

# CV CELL CYCLE TIME	
width = 2
height = width/1.618
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)
    ax.spines[['right', 'top']].set_visible(False)
ax.plot(phase_density, cv_density_at_fixed_phase, 'o', color=cool_purple, markeredgecolor='black',markeredgewidth='0.3', markersize=8, label='density')
ax.plot(phase_mass, cv_mass_at_fixed_phase, 'o', color=cool_orange, markeredgecolor='black', markeredgewidth='0.3', markersize=8, label='mass')
ax.plot(phase_vol, cv_volume_at_fixed_phase, 'o', color=cool_blue, markeredgecolor='black', markeredgewidth='0.3', markersize=8, label='volume')
ax.set_ylim(0.10, 0.25)
ax.set_xticks([0.0, 0.4, 0.8])
ax.set_yticks([0.10, 0.18, 0.25])
ax.set_xlabel(r'cell cycle phase', fontsize=10)
ax.set_ylabel(r'CV' '\ncoefficient of \n variation', fontsize=10)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=8)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=8)
ax.legend(fontsize=5)
fig.savefig('Hela_CV_phase_full_cell_cycle_tau%.2fh.pdf' % tau, bbox_inches='tight')  
plt.close(fig) 			

# BOX PLOT DENSITY DISTRIBUTION
width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
boxplots = ax.boxplot([density_stat_around_birth, density_stat_around_eos, density_stat_around_mitosis], patch_artist=True,
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                medianprops=dict(color='black'),
                flierprops=dict(marker='o', markerfacecolor='white', markeredgecolor='gray', markeredgewidth=0.5, markersize=2, alpha=0.5, linestyle='none'))  
colors = ['#810f7c', '#8856a7', '#8c96c6']
for patch, color in zip(boxplots['boxes'], colors):
    patch.set_facecolor(color)
ax.set_ylabel(r'$\rho$ (pg $\mathrm{\mu m}^{-3}$)', fontsize=8)
ax.set_xticklabels(['Birth', 'EOS', 'Mitosis'], fontsize=8)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=8)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=8)
fig.savefig('Hela_boxplot_density_at_birth_eos_mitosis_tau%.2fh.pdf' % tau, bbox_inches='tight')   
plt.close(fig)

# BOX PLOT SUBTRACTED DENSITY DISTRIBUTION
width = 1.5
height = width/1.618
fig, ax = plt.subplots(figsize=(width, height), dpi=100)
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)  
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.5)
    ax.spines[['right', 'top']].set_visible(False)
boxplots = ax.boxplot([density_subtr_stat_around_birth, density_subtr_stat_around_eos, density_subtr_stat_around_mitosis], patch_artist=True,
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                medianprops=dict(color='black'),
                flierprops=dict(marker='o', markerfacecolor='white', markeredgecolor='gray', markeredgewidth=0.5, markersize=2, alpha=0.5, linestyle='none'))  
colors = ['#810f7c', '#8856a7', '#8c96c6']
for patch, color in zip(boxplots['boxes'], colors):
    patch.set_facecolor(color)
ax.set_ylabel(r'$\rho - \rho_{target}$' ' \n (pg $\mathrm{\mu m}^{-3}$)', fontsize=8)
ax.set_xticklabels(['Birth', 'EOS', 'Mitosis'], fontsize=8)
ax.tick_params(axis='x', length=4, width=1.5, labelsize=8)
ax.tick_params(axis='y', length=4, width=1.5, labelsize=8)
fig.savefig('Hela_boxplot_density_subtr_at_birth_eos_mitosis_tau%.2fh.pdf' % tau, bbox_inches='tight')   
plt.close(fig)
