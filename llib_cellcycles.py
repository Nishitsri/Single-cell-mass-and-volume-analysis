# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
import itertools
import statistics as stat
import llib_figures as figs
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt
from math import nan
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde

###############################################	

def remove_nan(v):
	return [x for x in v if x is not None ]

def line(x, a, b):
	return np.add(a, np.multiply(b, x))

def linear_slope(x, y):
	if len(x) > 1 and len(y) > 1:
		popt, pcov = curve_fit(line, x, y)
		return popt[1]
	else:
		print('Data must contain at least two points.')
		return -1	

def num_stable_difference(x1, x2):
	significant_decimal_digits = max(len(str(x1).rsplit('.')[1]), len(str(x2).rsplit('.')[1]))
	return round(x1-x2, significant_decimal_digits)

def double_inequality(x, lower_bound, upper_bound):
	if x < upper_bound and x >= lower_bound:
		return True
	else:
		return False	

def scatter_to_func(x, y):
	ind = sorted(range(len(x)), key=lambda k: x[k])
	yy = []
	for k in range(len(y)):
		yy.append(y[ind[k]])
	return (sorted(x), yy)

	return (xx, yy, dev_xx, dev_yy) 

def delta_frame_given_tau(times, tau):
	delta_frame = None
	for k in range(len(times)):
		if times[k]-times[0] == tau :
			delta_frame = k
			return delta_frame
	if not delta_frame :		
		return 'No increment with such duration.'		

###############################################	

#######################
####BASIC TIMESERIES###
#######################

def cell_cycle_time_series(filename):
	with open(filename, 'r+') as file: 
		line_index = 0
		cells = []
		cell_time = []
		cell_obs = []	
		for line in file: 
			data_point = line.split()

			if line_index == 0:
				prev_cell_id = data_point[0]
				cell_time.append(float(data_point[1]))
				cell_obs.append(float(data_point[2]))
			else:
				cell_id = data_point[0]

				if  cell_id == prev_cell_id :
					cell_time.append(float(data_point[1]))
					cell_obs.append(float(data_point[2]))
				else :
					cells.append([cell_time, cell_obs])
					cell_time = []
					cell_obs = []	
					cell_time.append(float(data_point[1]))
					cell_obs.append(float(data_point[2]))
					prev_cell_id = cell_id
			line_index += 1
	cells.append([cell_time, cell_obs])			
	return cells

def normalize_cell_cycle_time_series(cells_obs_exp):
	return [ [cells_obs_exp[i][0], list(np.divide(cells_obs_exp[i][1], cells_obs_exp[i][1][0]))] if cells_obs_exp[i][1][0] != 0 else [cells_obs_exp[i][0], [np.nan]*len(cells_obs_exp[i][0])] for i in range(len(cells_obs_exp)) ]	

def normalize_cell_cycle_time_series_at_specific_timeframe(cells_obs_exp, event_frames):
	return [ [cells_obs_exp[i][0], list(np.divide(cells_obs_exp[i][1], cells_obs_exp[i][1][event_frames[i]]))] for i in range(len(cells_obs_exp))]	

def synch_cell_cycle_time_series_at_specific_timeframes(cells_obs_exp, event_frames):
	return [ [ list(np.subtract(cells_obs_exp[i][0], cells_obs_exp[i][0][event_frames[i]])), cells_obs_exp[i][1]] for i in range(len(cells_obs_exp))]	

def stack_cell_cycles_on_top(cells_obs_exp, event_frames):
	return [[ cells_obs_exp[i][0], [ cells_obs_exp[i][1][j] + cells_obs_exp[i][1][event_frames[i]] if j > event_frames[i] else cells_obs_exp[i][1][j] for j in range(len(cells_obs_exp[i][0])) ] ] for i in range(len(cells_obs_exp))]

#def set_zero_t0(cells_obs_exp):
#	return [ [ list(np.subtract(cells_obs_exp[i][0], cells_obs_exp[i][0][0])), cells_obs_exp[i][1]] for i in range(len(cells_obs_exp))]	

def set_zero_t0(cells_obs_exp):
	return [ [ [ num_stable_difference(cells_obs_exp[i][0][t], cells_obs_exp[i][0][0]) for t in range(len(cells_obs_exp[i][0]))], cells_obs_exp[i][1]] for i in range(len(cells_obs_exp))]	

def cell_index_original_file(filename):
	with open(filename, 'r+') as file: 
		line_index = 0
		cell_lineage = []
		for line in file: 
			data_point = line.split()

			if line_index == 0:
				prev_cell_id = data_point[0]
				cell_lineage.append(data_point[0])
			else:
				cell_id = data_point[0]
				if  cell_id != prev_cell_id :
					prev_cell_id = cell_id
					cell_lineage.append(data_point[0])
			line_index += 1	
	return cell_lineage			

def reduced_cell_cycle_time_series(cells_mass_exp, cells_vol_exp, cells_area_exp, cells_density_exp, cells_hgem_exp, cells_id_exp, cells_id_exp_hgem):	
	reduced_cells_mass_exp = []
	reduced_cells_vol_exp = []
	reduced_cells_area_exp = []
	reduced_cells_density_exp = []
	reduced_cells_id_exp = []
	for hgem_id in cells_id_exp_hgem:
		for cell_id_index in range(len(cells_id_exp)):
				if cells_id_exp[cell_id_index] == hgem_id:
					reduced_cells_mass_exp.append(cells_mass_exp[cell_id_index])
					reduced_cells_vol_exp.append(cells_vol_exp[cell_id_index])
					reduced_cells_area_exp.append(cells_area_exp[cell_id_index])
					reduced_cells_density_exp.append(cells_density_exp[cell_id_index])
					reduced_cells_id_exp.append(cells_id_exp[cell_id_index])
	reduced_cells_hgem_exp = cells_hgem_exp	
	reduced_cells_id_exp_hgem = cells_id_exp_hgem	

	return [reduced_cells_mass_exp, reduced_cells_vol_exp, reduced_cells_area_exp, reduced_cells_density_exp, reduced_cells_id_exp, reduced_cells_hgem_exp, reduced_cells_id_exp_hgem]					

def outlier_removing(outlier_vector, cells_mass_exp, cells_vol_exp, cells_area_exp, cells_density_exp, cells_hgem_exp, cells_id_exp, cells_id_exp_hgem):
	for n in outlier_vector:
		del cells_mass_exp[n]
		del cells_vol_exp[n]
		del cells_area_exp[n]
		del cells_density_exp[n]
		del cells_id_exp[n]
		del cells_id_exp_hgem[n]
		del cells_hgem_exp[n]

def cells_obs_between_timepoints_single_cell(cells_obs, timepoints1, timepoints2):
	
	obs_time, obs = [], []		
	for k in range(len(cells_obs[1])):
		if cells_obs[0][k] > timepoints1 and cells_obs[0][k] <= timepoints2:
			obs_time.append(cells_obs[0][k])
			obs.append(cells_obs[1][k])
	cells_obs_bat = [obs_time, obs]	
	return cells_obs_bat	

def cells_obs_before_timepoint(cells_obs_exp, timepoints):
	cells_obs_exp_bt = []		
	for i in range(len(cells_obs_exp)):	
		obs_time_i = []
		obs_i = []
		for k in range(len(cells_obs_exp[i][1])):
				if cells_obs_exp[i][0][k] <= timepoints[i]:
					obs_time_i.append(cells_obs_exp[i][0][k])
					obs_i.append(cells_obs_exp[i][1][k])
		cells_obs_exp_bt.append([obs_time_i, obs_i])	
	return cells_obs_exp_bt	

def cells_obs_after_timepoint(cells_obs_exp, timepoints):
	cells_obs_exp_at = []		
	for i in range(len(cells_obs_exp)):	
		obs_time_i = []
		obs_i = []
		for k in range(len(cells_obs_exp[i][1])):
				if cells_obs_exp[i][0][k] > timepoints[i]:
					obs_time_i.append(cells_obs_exp[i][0][k])
					obs_i.append(cells_obs_exp[i][1][k])
		cells_obs_exp_at.append([obs_time_i, obs_i])	
	return cells_obs_exp_at		

def cells_obs_between_timepoints(cells_obs_exp, timepoints1, timepoints2):
	cells_obs_exp_bat = []		
	for i in range(len(cells_obs_exp)):	
		obs_time_i = []
		obs_i = []
		for k in range(len(cells_obs_exp[i][1])):
				if cells_obs_exp[i][0][k] > timepoints1[i] and cells_obs_exp[i][0][k] <= timepoints2[i]:
					obs_time_i.append(cells_obs_exp[i][0][k])
					obs_i.append(cells_obs_exp[i][1][k])
		cells_obs_exp_bat.append([obs_time_i, obs_i])	
	return cells_obs_exp_bat		

def cells_obs_before_mitosis(cells_obs_exp, mitotic_times):
	cells_obs_exp_bm = []		
	for i in range(len(cells_obs_exp)):	
		obs_time_i = []
		obs_i = []
		for k in range(len(cells_obs_exp[i][1])):
				if cells_obs_exp[i][0][k] < mitotic_times[i]:
					obs_time_i.append(cells_obs_exp[i][0][k])
					obs_i.append(cells_obs_exp[i][1][k])
		cells_obs_exp_bm.append([obs_time_i, obs_i])	
	return cells_obs_exp_bm	

def cells_obs_before_g1s(cells_obs_exp, gs1_transition_times):
	cells_obs_exp_bg1s = []		
	for i in range(len(cells_obs_exp)):	
		obs_time_i = []
		obs_i = []
		for k in range(len(cells_obs_exp[i][1])):
				if cells_obs_exp[i][0][k] < gs1_transition_times[i]:
					obs_time_i.append(cells_obs_exp[i][0][k])
					obs_i.append(cells_obs_exp[i][1][k])
		cells_obs_exp_bg1s.append([obs_time_i, obs_i])	
	return cells_obs_exp_bg1s	

def cells_obs_before_mitosis_after_g1s(cells_obs_exp, mitotic_times, gs1_transition_times):
	cells_obs_exp_bm_ag1s = []		
	for i in range(len(cells_obs_exp)):	
		obs_time_i = []
		obs_i = []
		for k in range(len(cells_obs_exp[i][1])):
				if cells_obs_exp[i][0][k] > gs1_transition_times[i] and cells_obs_exp[i][0][k] < mitotic_times[i]:
					obs_time_i.append(cells_obs_exp[i][0][k])
					obs_i.append(cells_obs_exp[i][1][k])
		cells_obs_exp_bm_ag1s.append([obs_time_i, obs_i])	
	return cells_obs_exp_bm_ag1s
		
def cells_obs_in_mitosis(cells_obs_exp, mitotic_times):
	cells_obs_exp_in_mitosis = []		
	for i in range(len(cells_obs_exp)):	
		obs_time_i = []
		obs_i = []
		for k in range(len(cells_obs_exp[i][1])):
				if cells_obs_exp[i][0][k] > mitotic_times[i]:
					obs_time_i.append(cells_obs_exp[i][0][k])
					obs_i.append(cells_obs_exp[i][1][k])
		cells_obs_exp_in_mitosis.append([obs_time_i, obs_i])	
	return cells_obs_exp_in_mitosis	


############################
####OBSERVABLE TIMESERIES###
############################		

def increment(cells_obs_exp, tau):
	cells_incr_obs_exp = []
	for cell_index in range(len(cells_obs_exp)):
		time_incr = []
		obs_incr = []
		for i in range(len(cells_obs_exp[cell_index][0])):
			for k in range(len(cells_obs_exp[cell_index][0])-i):
				if num_stable_difference(cells_obs_exp[cell_index][0][i+k], cells_obs_exp[cell_index][0][i]) == tau :
					time_incr.append(cells_obs_exp[cell_index][0][i])
					obs_incr.append(cells_obs_exp[cell_index][1][i+k] - cells_obs_exp[cell_index][1][i])
		cells_incr_obs_exp.append([time_incr, obs_incr])
	return cells_incr_obs_exp

def total_increment(cells_obs_exp, tau): # output x[tau] - x[0] for each cell
	nocells = len(cells_obs_exp)

	increm_vec = []
	for n in range(nocells):
		increm = nan
		for k in range(len(cells_obs_exp[n][0])):
			if num_stable_difference(cells_obs_exp[n][0][k], cells_obs_exp[n][0][0]) == tau :
				increm = cells_obs_exp[n][1][k] - cells_obs_exp[n][1][0]
		increm_vec.append(increm)
	return increm_vec

def total_increment_linear_fit(cells_obs_exp):
	return [linear_slope(cell[0], cell[1])*cell[0][len(cell[0])-1] if len(cell[0]) > 1 else 0 for cell in cells_obs_exp]	

def value_at_timepoint_single_track(cell, tau): # output x[tau] for each cell
	value = None
	for k in range(len(cell[0])):
		if cell[0][k] == tau :
			value = cell[1][k]
			break
	return value
	
def find_frame_of_timepoint(cells_obs_exp, timepoints):
	nocells = len(cells_obs_exp)

	frames = []
	for n in range(nocells):
		for k in range(len(cells_obs_exp[n][0])):
			if timepoints[n] == timepoints[n] :
				if cells_obs_exp[n][0][k] == timepoints[n] :
					frames.append(k)
			else:
				frames.append(None)		
	return frames

def value_at_timepoint(cells_obs_exp, tau): # output x[tau] for each cell
	nocells = len(cells_obs_exp)

	value_vec = []
	for n in range(nocells):
		value = nan
		for k in range(len(cells_obs_exp[n][0])):
			if cells_obs_exp[n][0][k] == tau :
				value = cells_obs_exp[n][1][k]
		value_vec.append(value)
	return value_vec	

def value_in_time_range(cells_obs_exp, tmin, tmax): # output x[tau] for each cell
	nocells = len(cells_obs_exp)
	
	value_vec = []
	for n in range(nocells):
		value = nan
		for k in range(len(cells_obs_exp[n][0])):
			if tmin[n] <= cells_obs_exp[n][0][k] <= tmax[n] :
				value = cells_obs_exp[n][1][k]
		value_vec.append(value)
	return value_vec				

def two_obs_value_at_timepoint(cells_obs1_exp, cells_obs2_exp, tau): # output x[tau] for each cell
	nocells = len(cells_obs1_exp)

	value_vec_obs1 = []
	value_vec_obs2 = []
	for n in range(nocells):
		for k in range(len(cells_obs1_exp[n][0])):
			if cells_obs1_exp[n][0][k] == tau and cells_obs1_exp[n][1][k] ==  cells_obs1_exp[n][1][k]:
				for j in range(len(cells_obs2_exp[n][0])):
					if cells_obs2_exp[n][0][j] == tau and cells_obs2_exp[n][1][j] ==  cells_obs2_exp[n][1][j]:
						value_vec_obs1.append(cells_obs1_exp[n][1][k])
						value_vec_obs2.append(cells_obs2_exp[n][1][j])
	return [value_vec_obs1, value_vec_obs2]	

def all_values(cells_obs_exp):
	return remove_nan([cell[1][i] for cell in cells_obs_exp for i in range(len(cell[1]))])

def all_values_at_timeframe(cells_obs_exp, timeframe):
	no_cells = len(cells_obs_exp)
	return [cells_obs_exp[n][1][timeframe[n]] for n in range(no_cells)]

def all_values_until_timeframe(cells_obs_exp, timeframe):
	no_cells = len(cells_obs_exp)
	return [cells_obs_exp[n][1][i] for n in range(no_cells) for i in range(len(cells_obs_exp[n][1])) if i < timeframe ]

def all_values_until_timepoint(cells_obs_exp, timepoint):
	no_cells = len(cells_obs_exp)
	return all_values([ [cells_obs_exp[n][0][:timepoint[n]],cells_obs_exp[n][1][:timepoint[n]] ] for n in range(no_cells)])

def timeseries_until_timepoint(cells_obs_exp, timepoint):
	no_cells = len(cells_obs_exp)
	return [ [cells_obs_exp[n][0][:timepoint[n]],cells_obs_exp[n][1][:timepoint[n]] ] for n in range(no_cells) ]

def timeseries_after_timepoint(cells_obs_exp, timepoint):
	no_cells = len(cells_obs_exp)
	return [ [cells_obs_exp[n][0][:timepoint[n]],cells_obs_exp[n][1][timepoint[n]:] ] for n in range(no_cells) ]

def ratio_tracks(cells_obs_exp1, cells_obs_exp2):
	ratio = []
	for cell_index in range(len(cells_obs_exp1)):
		cell_ratio_track = [[], []]
		for t in range(min(len(cells_obs_exp1[cell_index][0]), len(cells_obs_exp2[cell_index][0]))):
			if cells_obs_exp1[cell_index][0][t] == cells_obs_exp1[cell_index][0][t]:
				cell_ratio_track[0].append(cells_obs_exp1[cell_index][0][t])
				cell_ratio_track[1].append(cells_obs_exp1[cell_index][1][t]/cells_obs_exp2[cell_index][1][t])
		ratio.append(cell_ratio_track)		
	return ratio

def der(cells_obs_exp, tau):
	cells_der_obs_exp = []
	for cell_index in range(len(cells_obs_exp)):
		time_der = []
		obs_der = []
		for i in range(len(cells_obs_exp[cell_index][0])):
			for k in range(len(cells_obs_exp[cell_index][0])-i):
				if num_stable_difference(cells_obs_exp[cell_index][0][i+k], cells_obs_exp[cell_index][0][i]) == tau :
					time_der.append(cells_obs_exp[cell_index][0][i])
					obs_der.append((cells_obs_exp[cell_index][1][i+k] - cells_obs_exp[cell_index][1][i])/tau)
		cells_der_obs_exp.append([time_der, obs_der])
	return cells_der_obs_exp	

def der_linear_fit(cells_obs_exp, tau):
	cells_der_obs_exp = []
	for cell_index in range(len(cells_obs_exp)):
		time_der = []
		obs_der = []
		for i in range(len(cells_obs_exp[cell_index][0])):
			for k in range(len(cells_obs_exp[cell_index][0])-i):
				if num_stable_difference(cells_obs_exp[cell_index][0][i+k], cells_obs_exp[cell_index][0][i]) == tau :
					time_der.append(cells_obs_exp[cell_index][0][i])
					x = [cells_obs_exp[cell_index][0][j] for j in range(i, i+k+1)] 
					y = [cells_obs_exp[cell_index][1][j] for j in range(i, i+k+1)]
					obs_der.append(linear_slope(x, y))
		cells_der_obs_exp.append([time_der, obs_der])
	return cells_der_obs_exp		

def growth_rate(cells_obs_exp, tau):
	cells_growth_rate_obs_exp = []
	for cell_index in range(len(cells_obs_exp)):
		time_growth_rate = []
		obs_growth_rate = []
		for i in range(len(cells_obs_exp[cell_index][0])):
			for k in range(len(cells_obs_exp[cell_index][0])-i):
				if num_stable_difference(cells_obs_exp[cell_index][0][i+k], cells_obs_exp[cell_index][0][i]) == tau :
					time_growth_rate.append(cells_obs_exp[cell_index][0][i])
					obs_growth_rate.append(((cells_obs_exp[cell_index][1][i+k] - cells_obs_exp[cell_index][1][i])/cells_obs_exp[cell_index][1][i])/tau)
		cells_growth_rate_obs_exp.append([time_growth_rate, obs_growth_rate])
	return cells_growth_rate_obs_exp	

def growth_rate_linear_fit(cells_obs_exp, tau):
	cells_growth_rate_obs_exp = []
	for cell_index in range(len(cells_obs_exp)):
		time_growth_rate = []
		obs_growth_rate = []
		for i in range(len(cells_obs_exp[cell_index][0])):
			for k in range(len(cells_obs_exp[cell_index][0])-i):
				if num_stable_difference(cells_obs_exp[cell_index][0][i+k], cells_obs_exp[cell_index][0][i]) == tau :
					time_growth_rate.append(cells_obs_exp[cell_index][0][i])
					x = [cells_obs_exp[cell_index][0][j] for j in range(i, i+k+1)] 
					y = [cells_obs_exp[cell_index][1][j] for j in range(i, i+k+1)]
					if cells_obs_exp[cell_index][1][i] != 0 :
						obs_growth_rate.append(linear_slope(x, y)/cells_obs_exp[cell_index][1][i])
					else :
						obs_growth_rate.append(np.nan)	
		cells_growth_rate_obs_exp.append([time_growth_rate, obs_growth_rate])
	return cells_growth_rate_obs_exp	

def growth_rate_linear_fit_without_smoothing(cells_obs_exp, tau): # this simply means that it's not a moving window, i.e., I associate a derivative to a time interval,
	cells_growth_rate_obs_exp = []							      # not to a single point (as I'm smoothing around the point anyway when I take a finite step)
	for cell_index in range(len(cells_obs_exp)):
		time_growth_rate = []
		obs_growth_rate = []
		i = 0
		while i < len(cells_obs_exp[cell_index][0]):
			found_right_increment = 0 
			for k in range(len(cells_obs_exp[cell_index][0])-i):
				if num_stable_difference(cells_obs_exp[cell_index][0][i+k], cells_obs_exp[cell_index][0][i]) == tau :
					found_right_increment = 1 # this is flag that tells you whether such two points do exist
					time_growth_rate.append(cells_obs_exp[cell_index][0][i])
					x = [cells_obs_exp[cell_index][0][j] for j in range(i, i+k+1)]
					y = [cells_obs_exp[cell_index][1][j] for j in range(i, i+k+1)]
					obs_growth_rate.append(linear_slope(x, y)/cells_obs_exp[cell_index][1][i])
					new_i = i + k

			if cells_obs_exp[cell_index][0][len(cells_obs_exp[cell_index][0])-1] -  cells_obs_exp[cell_index][0][new_i]	and found_right_increment >= tau:
				i  = new_i
			else : 
				i += 1								

		cells_growth_rate_obs_exp.append([time_growth_rate, obs_growth_rate])

	return cells_growth_rate_obs_exp

###############################
####STATISTICS AT FIXED TIME###
###############################
		
def obs_stat_at_fixed_time(cells_obs_exp, dt): # this will compute all the statistics for each n*dt time point
	tmax_obs = max([cells_obs_exp[i][0][j] for i in range(len(cells_obs_exp)) for j in range(len(cells_obs_exp[i][0]))])
	tmin_obs = min(([cells_obs_exp[i][0][j] for i in range(len(cells_obs_exp)) for j in range(len(cells_obs_exp[i][0]))]))
	time_obs = [dt*i - tmin_obs for i in range(int((tmax_obs-tmin_obs)/dt))] ## MODIFY THIS 
	obs_stat_at_fixed_time = []
	for t in time_obs:
		obs_t = []
		for i in range(len(cells_obs_exp)): 
			for j in range(len(cells_obs_exp[i][0])):
				if cells_obs_exp[i][0][j] == t and cells_obs_exp[i][1][j] == cells_obs_exp[i][1][j]:
					obs_t.append(cells_obs_exp[i][1][j])		
		obs_stat_at_fixed_time.append(obs_t)	
	return time_obs, obs_stat_at_fixed_time	


def obs_stat_at_fixed_phase(cells_obs_exp, phase_step):
	phase_obs = [phase_step*i for i in range(1, int(1/phase_step))]
	#phase_obs_centres = (phase_obs[:-1] + phase_obs[1:])/2. 
	obs_stat_at_fixed_phase = []
	for phase in phase_obs:
		obs_phase = []
		for i in range(len(cells_obs_exp)): 
			maxframe = len(cells_obs_exp[i][0])-1
			for j in range(len(cells_obs_exp[i][0])):
				phase_ij = cells_obs_exp[i][0][j]/cells_obs_exp[i][0][maxframe]
				if phase_ij <= phase and phase_ij > phase - phase_step:
					obs_phase.append(cells_obs_exp[i][1][j])		
		obs_stat_at_fixed_phase.append(obs_phase)	
	return phase_obs, obs_stat_at_fixed_phase	

def obs_stat_at_fixed_time_fixed_bins(cells_obs_exp, b):
	t = []
	x = []
	for cell in cells_obs_exp:
		t += [cell[0][i] for i in range(len(cell[0])) if cell[1][i] == cell[1][i]]
		x += [cell[1][i] for i in range(len(cell[1])) if cell[1][i] == cell[1][i]]	
	(t, x) = scatter_to_func(t, x)	

	n = len(x)
	tt_intervals = []
	xx_intervals = []
	k = 0
	while k*b+b < n :
		tt_bin = []
		xx_bin = []           
		for i in range(k*b, k*b+b):
			tt_bin.append(t[i])
			xx_bin.append(x[i]) 
		tt_intervals.append(stat.mean(tt_bin))      
		xx_intervals.append(xx_bin)   
		k += 1

	return tt_intervals, xx_intervals      

def obs_stat_at_fixed_phase_with_division_times(cells_obs_exp, division_times, phase_step):
	phase_obs = [phase_step*i for i in range(1, int(1/phase_step))]
	obs_stat_at_fixed_phase = []
	for phase in phase_obs:
		obs_phase = []
		for i in range(len(cells_obs_exp)): 
			maxframe = len(cells_obs_exp[i][0])-1
			for j in range(len(cells_obs_exp[i][0])):
				phase_ij = cells_obs_exp[i][0][j]/division_times[i]
				if phase_ij <= phase and phase_ij > phase - phase_step:
					obs_phase.append(cells_obs_exp[i][1][j])		
		obs_stat_at_fixed_phase.append(obs_phase)	
	phase_nonone = [phase_obs[i] for i in range(len(phase_obs)) if len(obs_stat_at_fixed_phase[i])>2]
	obs_stat_at_fixed_phase_nonone = [obs_stat_at_fixed_phase[i] for i in range(len(phase_obs)) if len(obs_stat_at_fixed_phase[i])>2]	
	return phase_nonone, obs_stat_at_fixed_phase_nonone 

def obs_stat(cells_obs_exp):
	return list(itertools.chain.from_iterable([cells_obs_exp[nocell][1] for nocell in range(len(cells_obs_exp))]))
		
def obs_timeavg_stat(cells_obs_exp):	
	return [sum(cells_obs_exp[nocell][1])/len(cells_obs_exp[nocell][1]) for nocell in range(len(cells_obs_exp))]

def mean_in_real_time(cells_obs_exp, dt):
	time, stat_at_fixed_time = obs_stat_at_fixed_time(cells_obs_exp, dt)
	return [time, [stat.mean(stat_at_fixed_time[k]) for k in range(len(time))]]	

def obs_track_time_avg(cells_obs_exp):
	nocells = len(cells_obs_exp)
	obs_time_avg_vec = []
	for n in range(nocells):
		noframes = len(cells_obs_exp[n][1])
		if noframes != 0 :
			obs_time_avg = sum([ cells_obs_exp[n][1][j] for j in range(noframes)])/noframes
			obs_time_avg_vec.append(obs_time_avg)	
	return obs_time_avg_vec	

def obs_track_variability(cells_obs_exp):
	nocells = len(cells_obs_exp)
	track_variability = []
	for n in range(nocells):
		noframes = len(cells_obs_exp[n][1])
		if noframes != 0 :
			obs_time_avg = sum([ cells_obs_exp[n][1][j] for j in range(noframes)])/noframes
			obs_time_var = sqrt(sum([ (cells_obs_exp[n][1][j]-obs_time_avg)*(cells_obs_exp[n][1][j]-obs_time_avg) for j in range(len(cells_obs_exp[n][1])) ])/noframes)/abs(obs_time_avg)
			track_variability.append(obs_time_var)
	return track_variability		

#############################
####IDENTIFY SPECIAL TIMES###
#############################

def division_times(cells_obs_exp):
	division_frames = [len(cells_obs_exp[cell_index][0])-1 for cell_index in range(len(cells_obs_exp))]
	division_times = [cells_obs_exp[cell_index][0][len(cells_obs_exp[cell_index][0])-1] for cell_index in range(len(cells_obs_exp)) ]				
	return division_frames, division_times

def cells_obs_at_division(cells_obs_exp):
	div_times, no_cells = division_times(cells_obs_exp), len(cells_obs_exp)
	return [cells_obs_exp[n][1][div_times[0][n]] for n in range(no_cells)]

def osmotic_shock_times(cells_vol_exp, tau): # here we take the maximum value as shock time;
	no_cells = len(cells_vol_exp)			 # this is not the real start of the shock; use the der ratio func for that
	delta_frame_tau = delta_frame_given_tau(cells_vol_exp[0][0], tau)	

	cells_vol_growth_rate = growth_rate(cells_vol_exp, tau)
	max_growth_rate = [ max(cells_vol_growth_rate[i][1]) for i in range(no_cells) ]

	index_max_growth_rate = [ cells_vol_growth_rate[i][1].index(max_growth_rate[i]) + delta_frame_tau for i in range(no_cells)]
	osmotic_shock_times = [ cells_vol_growth_rate[i][0][index_max_growth_rate[i]] for i in range(no_cells) ]

	return index_max_growth_rate, osmotic_shock_times

def mitotic_times(cells_vol_exp, tau):
	no_cells = len(cells_vol_exp)
	cells_vol_growth_rate = growth_rate(cells_vol_exp, tau)
	div_times = division_times(cells_vol_exp)[1]

	time_threShold = 2/3
	adjusted_cells_vol_growth_rate = [
	[		[cells_vol_growth_rate[i][0][frame] for frame in range(len(cells_vol_growth_rate[i][1])) 
				if cells_vol_growth_rate[i][0][frame]/div_times[i] >= time_threShold],
			[cells_vol_growth_rate[i][1][frame] for frame in range(len(cells_vol_growth_rate[i][1])) 
				if cells_vol_growth_rate[i][0][frame]/div_times[i] >= time_threShold]] 
		for i in range(no_cells)	
	]	
	max_growth_rate = [ max(adjusted_cells_vol_growth_rate[i][1]) for i in range(no_cells) ]
	index_max_growth_rate = [ adjusted_cells_vol_growth_rate[i][1].index(max_growth_rate[i]) for i in range(no_cells)]
	mitotic_times = [ adjusted_cells_vol_growth_rate[i][0][index_max_growth_rate[i]] for i in range(no_cells) ]

	# NB: you need this index vector because you've used an adjusted vector before
	index_vol_at_mitotic_times = [i for n in range(no_cells) 
					for i in range(len(cells_vol_exp[n][1])) if cells_vol_exp[n][0][i] == mitotic_times[n]]
				
	return index_vol_at_mitotic_times, mitotic_times

def mother_mitotic_times(cells_vol_exp, tau):
	# NB: this function finds mitosis of the mother cell; input data must be syncronized at division, so that it looks for mitotis only for times before zero
	no_cells = len(cells_vol_exp)
	cells_vol_growth_rate = growth_rate(cells_vol_exp, tau)
	div_times = division_times(cells_vol_exp)[1]

	adjusted_cells_vol_growth_rate = [
	[		[cells_vol_growth_rate[i][0][frame] for frame in range(len(cells_vol_growth_rate[i][1])) 
				if cells_vol_growth_rate[i][0][frame] <= 0 ],
			[cells_vol_growth_rate[i][1][frame] for frame in range(len(cells_vol_growth_rate[i][1])) 
				if cells_vol_growth_rate[i][0][frame] <= 0 ]] 
		for i in range(no_cells)	
	]	

	max_growth_rate = [ max(adjusted_cells_vol_growth_rate[i][1]) for i in range(no_cells) ]
	index_max_growth_rate = [ adjusted_cells_vol_growth_rate[i][1].index(max_growth_rate[i]) for i in range(no_cells)]
	mitotic_times = [ adjusted_cells_vol_growth_rate[i][0][index_max_growth_rate[i]] for i in range(no_cells) ]

	# NB: you need this index vector because you've used an adjusted vector before
	index_vol_at_mitotic_times = [i for n in range(no_cells) 
					for i in range(len(cells_vol_exp[n][1])) if cells_vol_exp[n][0][i] == mitotic_times[n]]
				
	return index_vol_at_mitotic_times, mitotic_times	

def cells_obs_at_mitosis(cells_obs_exp, cells_vol_exp, tau):
	mit_times, no_cells = mitotic_times(cells_vol_exp, tau), len(cells_obs_exp)
	'''occasionally the mitotic frame ientified by the volume does not exist in the other observables; we put None then'''
	return [cells_obs_exp[n][1][mit_times[0][n]] if mit_times[0][n] < len(cells_obs_exp[n][1]) else None for n in range(no_cells)]

def g1s_transition_times(cells_hgem_exp, tau):
	nocells = len(cells_hgem_exp)
	cells_hgem_growth_rate = growth_rate(cells_hgem_exp, tau)

	max_growth_rate = [ max(cells_hgem_growth_rate[i][1]) for i in range(nocells) ]
	index_max_growth_rate = [ cells_hgem_growth_rate[i][1].index(max_growth_rate[i]) for i in range(nocells)]
	g1s_transition_times = [ cells_hgem_growth_rate[i][0][index_max_growth_rate[i]] for i in range(nocells) ]

	return index_max_growth_rate, g1s_transition_times	

def cells_obs_at_g1s(cells_obs_exp, cells_hgem_exp, tau):
	g1s_times, no_cells = g1s_transition_times(cells_hgem_exp, tau), len(cells_obs_exp)
	return [cells_obs_exp[n][1][g1s_times[0][n]] for n in range(no_cells)]	

def minimum_volume_times(cells_vol_exp):
	index_min_vol = [ cell_vol[1].index(min(cell_vol[1])) for cell_vol in cells_vol_exp]
	time_min_vol = [ cell_vol[0][cell_vol[1].index(min(cell_vol[1]))] for cell_vol in cells_vol_exp]
	return index_min_vol, time_min_vol

def observable_slope_change_times_ratio_method_single_cell(cells_obs, tau, time_window, signa_delta_minus, signa_delta_plus):

	eps = 0.00000001 # this is to regularize stuff

	delta_obs_minus_vec, delta_obs_plus_vec, der_ratio_vec, time_vec, frame_vec = [], [], [], [], []
	for i in range(len(cells_obs[0])): # go through the entire timeseries
		for k in range(len(cells_obs[0])-i): # this should automatically take care of the border
			delta_t_minus = num_stable_difference(cells_obs[0][i], cells_obs[0][i-k]) 
			delta_t_plus = num_stable_difference(cells_obs[0][i+k], cells_obs[0][i])
			if i-k > 0 and delta_t_minus == tau and delta_t_plus == tau:
				x_plus = [cells_obs[0][j] for j in range(i, i+k+1)] 
				y_plus = [cells_obs[1][j] for j in range(i, i+k+1)]
				delta_obs_plus = linear_slope(x_plus, y_plus)*tau # estimating the increment with a linear fit is more robus to noisy single points
				x_minus = [cells_obs[0][j] for j in range(i-k, i+1)] 
				y_minus = [cells_obs[1][j] for j in range(i-k, i+1)]
				delta_obs_minus = linear_slope(x_minus, y_minus)*tau	
				if delta_obs_minus*signa_delta_minus >= 0 and delta_obs_plus*signa_delta_plus >= 0 and delta_obs_minus > delta_obs_plus: # check both increment are positive and that the 1st one is bigger than the 2nd one
					der_ratio = delta_obs_plus/delta_obs_minus
					time_der_ratio = cells_obs[0][i]
					delta_obs_minus_vec.append(delta_obs_minus)
					delta_obs_plus_vec.append(delta_obs_plus)
					der_ratio_vec.append(der_ratio)
					time_vec.append(time_der_ratio)
					frame_vec.append(i)

	der_ratio_in_time_window = [der_ratio_vec[i] for i in range(len(der_ratio_vec)) if time_vec[i] <= time_window[1] and time_vec[i] >= time_window[0]]	
	if len(der_ratio_in_time_window) > 0 :
		min_der_ratio = min(der_ratio_in_time_window)
		index_min_der_ratio = der_ratio_vec.index(min_der_ratio)
		time_min_der_ratio = time_vec[index_min_der_ratio]
		frame_min_der_ratio = frame_vec[index_min_der_ratio]
		delta_obs_minus_at_slope_change = delta_obs_minus_vec[index_min_der_ratio]
		delta_obs_plus_at_slope_change = delta_obs_plus_vec[index_min_der_ratio]
		der_ratio_at_slope_change = min_der_ratio
		slope_change_time = time_min_der_ratio
		slope_change_frame = frame_min_der_ratio
	else :
		delta_obs_minus_at_slope_change = None
		delta_obs_plus_at_slope_change = None
		der_ratio_at_slope_change = None
		slope_change_time = None
		slope_change_frame = None

	return slope_change_time, slope_change_frame, delta_obs_minus_at_slope_change, delta_obs_plus_at_slope_change, der_ratio_at_slope_change

def observable_slope_change_times_ratio_method(cells_obs_exp, tau, time_window, signa_delta_minus, signa_delta_plus):	# look for the ratio of the derivative to find the slope change
	nocells = len(cells_obs_exp)
	eps = 0.00000001 # this is to regularize stuff

	slope_change_time = []
	slope_change_frame = []
	delta_obs_minus_at_slope_change = []
	delta_obs_plus_at_slope_change = []
	der_ratio_at_slope_change = []

	for cell_index in range(nocells): # go through each each individually
		delta_obs_minus_vec = []
		delta_obs_plus_vec = []
		der_ratio_vec = []
		time_vec = []
		frame_vec = []

		for i in range(len(cells_obs_exp[cell_index][0])): # go through the entire timeseries
			for k in range(len(cells_obs_exp[cell_index][0])-i): # this should automatically take care of the border
				delta_t_minus = num_stable_difference(cells_obs_exp[cell_index][0][i], cells_obs_exp[cell_index][0][i-k]) 
				delta_t_plus = num_stable_difference(cells_obs_exp[cell_index][0][i+k], cells_obs_exp[cell_index][0][i])
				if i-k > 0 and delta_t_minus == tau and delta_t_plus == tau:
					x_plus = [cells_obs_exp[cell_index][0][j] for j in range(i, i+k+1)] 
					y_plus = [cells_obs_exp[cell_index][1][j] for j in range(i, i+k+1)]
					delta_obs_plus = linear_slope(x_plus, y_plus)*tau # estimating the increment with a linear fit is more robus to noisy single points
					x_minus = [cells_obs_exp[cell_index][0][j] for j in range(i-k, i+1)] 
					y_minus = [cells_obs_exp[cell_index][1][j] for j in range(i-k, i+1)]
					delta_obs_minus = linear_slope(x_minus, y_minus)*tau	
					if delta_obs_minus*signa_delta_minus >= 0 and delta_obs_plus*signa_delta_plus >= 0 :
						der_ratio = delta_obs_plus/delta_obs_minus
						time_der_ratio = cells_obs_exp[cell_index][0][i]
						delta_obs_minus_vec.append(delta_obs_minus)
						delta_obs_plus_vec.append(delta_obs_plus)
						der_ratio_vec.append(der_ratio)
						time_vec.append(time_der_ratio)
						frame_vec.append(i)

		der_ratio_in_time_window = [der_ratio_vec[i] for i in range(len(der_ratio_vec)) if time_vec[i] <= time_window[1] and time_vec[i] >= time_window[0]]	
		if len(der_ratio_in_time_window) > 0 :
			min_der_ratio = min(der_ratio_in_time_window)
			index_min_der_ratio = der_ratio_vec.index(min_der_ratio)
			time_min_der_ratio = time_vec[index_min_der_ratio]
			frame_min_der_ratio = frame_vec[index_min_der_ratio]
			delta_obs_minus_at_slope_change.append(delta_obs_minus_vec[index_min_der_ratio])
			delta_obs_plus_at_slope_change.append(delta_obs_plus_vec[index_min_der_ratio])
			der_ratio_at_slope_change.append(min_der_ratio)
			slope_change_time.append(time_min_der_ratio)
			slope_change_frame.append(frame_min_der_ratio)
		else :
			delta_obs_minus_at_slope_change.append(None)
			delta_obs_plus_at_slope_change.append(None)
			der_ratio_at_slope_change.append(None)
			slope_change_time.append(None)	
			slope_change_frame.append(None)

	return slope_change_time, slope_change_frame, delta_obs_minus_at_slope_change, delta_obs_plus_at_slope_change, der_ratio_at_slope_change

def observable_slope_change_times_diff_method(cells_obs_exp, tau, time_window, signa_delta_minus, signa_delta_plus):	# look for the difference of the derivative to find the slope change
	nocells = len(cells_obs_exp)
	eps = 0.00000001 # this is to regularize stuff

	slope_change_time = []
	slope_change_frame = []
	delta_obs_minus_at_slope_change = []
	delta_obs_plus_at_slope_change = []
	der_diff_at_slope_change = []

	for cell_index in range(nocells): # go through each each individually
		delta_obs_minus_vec = []
		delta_obs_plus_vec = []
		der_diff_vec = []
		time_vec = []
		frame_vec = []

		for i in range(len(cells_obs_exp[cell_index][0])): # go through the entire timeseries
			for k in range(len(cells_obs_exp[cell_index][0])-i): # this should automatically take care of the border
				delta_t_minus = num_stable_difference(cells_obs_exp[cell_index][0][i], cells_obs_exp[cell_index][0][i-k]) 
				delta_t_plus = num_stable_difference(cells_obs_exp[cell_index][0][i+k], cells_obs_exp[cell_index][0][i])
				if i-k > 0 and delta_t_minus == tau and delta_t_plus == tau:
					x_plus = [cells_obs_exp[cell_index][0][j] for j in range(i, i+k+1)] 
					y_plus = [cells_obs_exp[cell_index][1][j] for j in range(i, i+k+1)]
					delta_obs_plus = linear_slope(x_plus, y_plus)*tau # estimating the increment with a linear fit is more robus to noisy single points
					x_minus = [cells_obs_exp[cell_index][0][j] for j in range(i-k, i+1)] 
					y_minus = [cells_obs_exp[cell_index][1][j] for j in range(i-k, i+1)]
					delta_obs_minus = linear_slope(x_minus, y_minus)*tau	
					if delta_obs_minus*signa_delta_minus >= 0 and delta_obs_plus*signa_delta_plus >= 0 :
						der_diff = delta_obs_minus - delta_obs_plus
						time_der_diff = cells_obs_exp[cell_index][0][i]
						delta_obs_plus_vec.append(delta_obs_plus)
						der_diff_vec.append(der_diff)
						delta_obs_minus_vec.append(delta_obs_minus)
						time_vec.append(time_der_diff)
						frame_vec.append(i)

		der_diff_in_time_window = [der_diff_vec[i] for i in range(len(der_diff_vec)) if time_vec[i] <= time_window[1] and time_vec[i] >= time_window[0]]	
		if len(der_diff_in_time_window) > 0 :
			max_der_diff = max(der_diff_in_time_window)
			index_max_der_diff = der_diff_vec.index(max_der_diff)
			time_max_der_diff = time_vec[index_max_der_diff]
			frame_max_der_diff = frame_vec[index_max_der_diff]
			delta_obs_minus_at_slope_change.append(delta_obs_minus_vec[index_max_der_diff])
			delta_obs_plus_at_slope_change.append(delta_obs_plus_vec[index_max_der_diff])
			der_diff_at_slope_change.append(max_der_diff)
			slope_change_time.append(time_max_der_diff)
			slope_change_frame.append(frame_max_der_diff)
		else :
			delta_obs_minus_at_slope_change.append(None)
			delta_obs_plus_at_slope_change.append(None)
			der_diff_at_slope_change.append(None)
			slope_change_time.append(None)	
			slope_change_frame.append(None)

	return slope_change_time, slope_change_frame, delta_obs_minus_at_slope_change, delta_obs_plus_at_slope_change, der_diff_at_slope_change		

def scatter_plot_vectors_single_cell(cells_obs_x, cells_obs_y):
	x_scatter_plot = []
	y_scatter_plot = []
	for i in range(min(len(cells_obs_x[0]), len(cells_obs_y[0]))):
		if cells_obs_y[0][i] == cells_obs_x[0][i] :
			x_scatter_plot.append(cells_obs_x[1][i])
			y_scatter_plot.append(cells_obs_y[1][i])						
	correlation = np.corrcoef(x_scatter_plot, y_scatter_plot)[0,1]
	return [x_scatter_plot, y_scatter_plot], correlation	

def scatter_plot_vectors(cells_obs_exp_x, cells_obs_exp_y):
	x_scatter_plot = []
	y_scatter_plot = []
	nocells = len(cells_obs_exp_x)
	for nocell in range(nocells):
		for i in range(min(len(cells_obs_exp_y[nocell][0]), len(cells_obs_exp_x[nocell][0]))):
			if cells_obs_exp_y[nocell][0][i] == cells_obs_exp_x[nocell][0][i] :
			#if num_stable_difference(cells_obs_exp_y[nocell][0][i], cells_obs_exp_x[nocell][0][i]) == 0:
				x_scatter_plot.append(cells_obs_exp_x[nocell][1][i])
				y_scatter_plot.append(cells_obs_exp_y[nocell][1][i])						
	correlation = np.corrcoef(x_scatter_plot, y_scatter_plot)[0,1]
	return [x_scatter_plot, y_scatter_plot], correlation	

def scatter_plot_vectors_delayed(cells_obs_exp_x, cells_obs_exp_y, tau):
	x_scatter_plot = []
	y_scatter_plot = []
	nocells = len(cells_obs_exp_x)
	for nocell in range(nocells):
		for i in range(min(len(cells_obs_exp_y[nocell][0]), len(cells_obs_exp_x[nocell][0]))):
			for j in range(i, max(len(cells_obs_exp_y[nocell][0]), len(cells_obs_exp_x[nocell][0]))):
				if cells_obs_exp_y[nocell][0][j] == cells_obs_exp_x[nocell][0][i] + tau :
					x_scatter_plot.append(cells_obs_exp_x[nocell][1][i])
					y_scatter_plot.append(cells_obs_exp_y[nocell][1][j])
					break						
	correlation = np.corrcoef(x_scatter_plot, y_scatter_plot)[0,1]
	return [x_scatter_plot, y_scatter_plot], correlation


def scatter_plot_3D_vectors(cells_obs_exp_x, cells_obs_exp_y, cells_obs_exp_z):
	x_scatter_plot = []
	y_scatter_plot = []
	z_scatter_plot = []
	nocells = len(cells_obs_exp_x)
	for nocell in range(nocells):
		for i in range(min(len(cells_obs_exp_x[nocell][0]), len(cells_obs_exp_y[nocell][0]), len(cells_obs_exp_z[nocell][0]))):
			if cells_obs_exp_z[nocell][0][i] == cells_obs_exp_y[nocell][0][i] == cells_obs_exp_x[nocell][0][i] :
				x_scatter_plot.append(cells_obs_exp_x[nocell][1][i])
				y_scatter_plot.append(cells_obs_exp_y[nocell][1][i])	
				z_scatter_plot.append(cells_obs_exp_z[nocell][1][i])					
	return [x_scatter_plot, y_scatter_plot, z_scatter_plot]

def scatter_plot_vectors_with_name(cells_obs_exp_x, cells_obs_exp_y, cells_obs_id):
	x_scatter_plot, y_scatter_plot, name_scatter_plot = [], [], []
	nocells = len(cells_obs_exp_x)
	for nocell in range(nocells):
		for i in range(min(len(cells_obs_exp_y[nocell][0]), len(cells_obs_exp_x[nocell][0]))):
			if cells_obs_exp_y[nocell][0][i] == cells_obs_exp_x[nocell][0][i] :
				x_scatter_plot.append(cells_obs_exp_x[nocell][1][i])
				y_scatter_plot.append(cells_obs_exp_y[nocell][1][i])	
				name_scatter_plot.append(cells_obs_id[nocell])			
	return [x_scatter_plot, y_scatter_plot, name_scatter_plot]	

def scatter_plot_vectors_3D(cells_obs_exp_x, cells_obs_exp_y, cells_obs_exp_z):
	x_scatter_plot, y_scatter_plot, z_scatter_plot = [], [], []
	nocells = len(cells_obs_exp_x)
	for nocell in range(nocells):
		for i in range(min(len(cells_obs_exp_x[nocell][0]), len(cells_obs_exp_y[nocell][0]), len(cells_obs_exp_z[nocell][0]))):
			if len(set([cells_obs_exp_z[nocell][0][i], cells_obs_exp_y[nocell][0][i], cells_obs_exp_x[nocell][0][i]])) == 1 :
				x_scatter_plot.append(cells_obs_exp_x[nocell][1][i])
				y_scatter_plot.append(cells_obs_exp_y[nocell][1][i])
				z_scatter_plot.append(cells_obs_exp_z[nocell][1][i])						
	return [x_scatter_plot, y_scatter_plot, z_scatter_plot]		

def scatter_plot_vectors_in_time_interval(cells_obs_exp_x, cells_obs_exp_y, time_lim):
	x_scatter_plot = []
	y_scatter_plot = []
	nocells = len(cells_obs_exp_x)
	for nocell in range(nocells):
		for i in range(min(len(cells_obs_exp_y[nocell][0]), len(cells_obs_exp_x[nocell][0]))):
			if cells_obs_exp_y[nocell][0][i] == cells_obs_exp_x[nocell][0][i] :
				if double_inequality(cells_obs_exp_x[nocell][0][i], time_lim[0], time_lim[1]) and double_inequality(cells_obs_exp_x[nocell][0][i], time_lim[0], time_lim[1]):
					x_scatter_plot.append(cells_obs_exp_x[nocell][1][i])
					y_scatter_plot.append(cells_obs_exp_y[nocell][1][i])		
	if len(x_scatter_plot) > 1 and len(y_scatter_plot) > 1:									
		correlation = np.corrcoef(x_scatter_plot, y_scatter_plot)[0,1]
		return [x_scatter_plot, y_scatter_plot], correlation	
	else :
		return False, False ## this is just a quick and dirty way to flag if the scatter plot fails		

def scatter_plot_vectors_with_time(cells_obs_exp_x, cells_obs_exp_y):
	x_scatter_plot = []
	y_scatter_plot = []
	time = []
	nocells = len(cells_obs_exp_x)
	for nocell in range(nocells):
		for i in range(min(len(cells_obs_exp_y[nocell][0]), len(cells_obs_exp_x[nocell][0]))):
			if cells_obs_exp_y[nocell][0][i] == cells_obs_exp_x[nocell][0][i] :
				x_scatter_plot.append(cells_obs_exp_x[nocell][1][i])
				y_scatter_plot.append(cells_obs_exp_y[nocell][1][i])	
				time.append(cells_obs_exp_y[nocell][0][i])
	return [x_scatter_plot, y_scatter_plot, time]		

##################
####PRINT STUFF###
##################


def print_cell_cycle_curves(cells_mass_exp, cells_vol_exp, cells_area_exp, cells_density_exp, cells_id_exp):
	nocells_exp = len(cells_mass_exp)
	pdf_pages = PdfPages('full_cell_cycles_curves.pdf')
	for i in range(nocells_exp):
		lineage_id = cells_id_exp[i]
		print(lineage_id, i)
		fig = figs.single_cell_panel(lineage_id,
			cells_mass[i], 
			cells_vol[i], 
			cells_area[i], 
			cells_density[i]
		 )
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 	

def print_normalized_cell_cycle_curves(cells_mass_exp, cells_vol_exp, cells_area_exp, cells_density_exp, cells_id_exp):
	nocells_exp = len(cells_mass_exp)
	pdf_pages = PdfPages('full_cell_cycles_curves.pdf')
	for i in range(nocells_exp):
		lineage_id = cells_id_exp[i]
		print(lineage_id, i)
		fig = figs.single_cell_panel(lineage_id,
			[cells_mass_exp[i][0], np.divide(cells_mass_exp[i][1], cells_mass_exp[i][1][0])],
			[cells_vol_exp[i][0], np.divide(cells_vol_exp[i][1], cells_vol_exp[i][1][0])],
			[cells_area_exp[i][0], np.divide(cells_area_exp[i][1], cells_area_exp[i][1][0])],
			[cells_density_exp[i][0], np.divide(cells_density_exp[i][1], cells_density_exp[i][1][0])] 
		 )
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 		

#######################
####BINNING FUNCTIONS###
#######################

def binningdata_bins_mean(x, y, b): #this one bins data in a fixed # of bins interval    
	(x, y) = scatter_to_func(x, y)
	
	n = len(x)
	xx_mean = []
	yy_mean = []
	dev_xx = []
	dev_yy = []
	k = 0
	while k*b+b < n :
		xx_bin = []           
		yy_bin = []
		for i in range(k*b, k*b+b):
			xx_bin.append(x[i])       
			yy_bin.append(y[i])
		xx_mean.append(stat.mean(xx_bin))
		yy_mean.append(stat.mean(yy_bin))   
		dev_xx.append(stat.stdev(xx_bin))
		dev_yy.append(stat.stdev(yy_bin)) 
		k += 1           

	return (xx_mean, yy_mean, dev_xx, dev_yy, np.divide(dev_xx, sqrt(b)), np.divide(dev_yy, sqrt(b))) 

def binningdata_bins_median(x, y, b): #this one bins data in a fixed # of bins interval    
	(x, y) = scatter_to_func(x, y)
	
	n = len(x)
	xx_median = []
	yy_median = []
	k = 0
	while k*b+b < n :
		xx_bin = []           
		yy_bin = []
		for i in range(k*b, k*b+b):
			xx_bin.append(x[i])       
			yy_bin.append(y[i])
		xx_median.append(stat.median(xx_bin))
		yy_median.append(stat.median(yy_bin))    
		k += 1    
			 
	return (xx_median, yy_median) 

############################
####BINNING FUNCTION LOOP###
############################

def scatter_plot_der_y_vs_variable_x(ensemble_x, ensemble_y, tau_vec, window_size):
	no_tau = len(tau_vec)

	# COMPUTE GROWTH RATE
	der_y_vs_tau = [ der_linear_fit(ensemble_y, tau) for tau in tau_vec ]

	# COMPUTE SCATTER PLOT
	der_y_vs_variable_x_vs_tau = [ scatter_plot_vectors(ensemble_x, der_y) for der_y in der_y_vs_tau ]

	# COMPUTE BINNED PLOT
	binned_der_y_vs_variable_x_vs_tau = [ binningdata_bins_mean(der_y_vs_variable_x_vs_tau[i][0], der_y_vs_variable_x_vs_tau[i][1], window_size) for i in range(no_tracks)]
	
	# COMPUTE LINEAR REGRESSION
	der_y_vs_variable_x_linreg_vs_tau = [ linregress(binned_der_y_vs_variable_x_vs_tau[i][0], binned_der_y_vs_variable_x_vs_tau[i][1]) for i in range(no_tau) ]
	intercept_vs_tau = [der_y_vs_variable_x_linreg_vs_tau[i].intercept for i in range(no_tau)]
	slope_vs_tau = [der_y_vs_variable_x_linreg_vs_tau[i].slope for i in range(no_tau)]
	rvalue_vs_tau = [der_y_vs_variable_x_linreg_vs_tau[i].rvalue for i in range(no_tau)]

	# PLOTTING
	figs.simple_plot('grate_windowsize%g_intercept.pdf' % (window_size), tau_vec, intercept_vs_tau, 'black')
	figs.simple_plot('grate_windowsize%g_slope.pdf' % (window_size), tau_vec, slope_vs_tau, 'black')
	figs.simple_plot('grate_windowsize%g_rvalue.pdf' % (window_size), tau_vec, rvalue_vs_tau, 'black')

	return der_y_vs_variable_x_vs_tau, binned_der_y_vs_variable_x_vs_tau, der_y_vs_variable_x_linreg_vs_tau


def scatter_plot_growth_rate_y_vs_variable_x_vs_der_step(ensemble_x, ensemble_y, tau_vec, window_size):
	no_tau = len(tau_vec)

	# COMPUTE GROWTH RATE
	growth_rate_y_vs_tau = [ growth_rate_linear_fit(ensemble_y, tau) for tau in tau_vec ]

	# COMPUTE SCATTER PLOT
	growth_rate_y_vs_variable_x_vs_tau = [ scatter_plot_vectors(ensemble_x, growth_rate_y) for growth_rate_y in growth_rate_y_vs_tau ]

	# COMPUTE BINNED PLOT
	binned_growth_rate_y_vs_variable_x_vs_tau = [ binningdata_bins_mean(growth_rate_y_vs_variable_x_vs_tau[i][0], growth_rate_y_vs_variable_x_vs_tau[i][1], window_size) for i in range(no_tracks)]
	
	# COMPUTE LINEAR REGRESSION
	growth_rate_y_vs_variable_x_linreg_vs_tau = [ linregress(binned_growth_rate_y_vs_variable_x_vs_tau[i][0], binned_growth_rate_y_vs_variable_x_vs_tau[i][1]) for i in range(no_tau) ]
	intercept_vs_tau = [growth_rate_y_vs_variable_x_linreg_vs_tau[i].intercept for i in range(no_tau)]
	slope_vs_tau = [growth_rate_y_vs_variable_x_linreg_vs_tau[i].slope for i in range(no_tau)]
	rvalue_vs_tau = [growth_rate_y_vs_variable_x_linreg_vs_tau[i].rvalue for i in range(no_tau)]

	# PLOTTING
	figs.simple_plot('grate_windowsize%g_intercept.pdf' % (window_size), tau_vec, intercept_vs_tau, 'black')
	figs.simple_plot('grate_windowsize%g_slope.pdf' % (window_size), tau_vec, slope_vs_tau, 'black')
	figs.simple_plot('grate_windowsize%g_rvalue.pdf' % (window_size), tau_vec, rvalue_vs_tau, 'black')

	return growth_rate_y_vs_variable_x_vs_tau, binned_growth_rate_y_vs_variable_x_vs_tau, growth_rate_y_vs_variable_x_linreg_vs_tau
