# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from math import pi
from math import sqrt

def closest_index(v, x):
	distance =[abs(v[i]-x) for i in range(len(v))]
	return distance.index(min(distance))

#######################################################

def simple_plot(filename, x, y, x_labelname, y_labelname, colorname):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.plot(x, y, 'o', color=colorname, markeredgecolor='black', markeredgewidth='1', markersize=15)
	ax.set_xlabel(r'%s' % x_labelname, fontsize=15)
	ax.set_ylabel(r'%s' % y_labelname, fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	

def simple_plot_with_label(filename, x, y, x_labelname, y_labelname, colorname, label):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.plot(x, y, 'o', color=colorname, markeredgecolor='black', markeredgewidth='1', markersize=15, label=label)
	ax.set_xlabel(r'%s' % x_labelname, fontsize=15)
	ax.set_ylabel(r'%s' % y_labelname, fontsize=15)
	ax.legend(frameon=True, fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 			

def simple_histogram(filename, data, dataname, colorname):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.hist(data, color=colorname)
	ax.set_xlabel(r'%s' % dataname, fontsize=15)
	ax.set_ylabel(r'occurrence', fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 

def simple_histogram_with_label(filename, data, dataname, colorname, label):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.hist(data, color=colorname, label=label)
	ax.set_xlabel(r'%s' % dataname, fontsize=15)
	ax.set_ylabel(r'occurrence', fontsize=15)
	ax.legend(frameon=True, fontsize=10)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 

def double_simple_histogram_with_label(filename, data1, data2, dataname, colorname1, colorname2, label1, label2):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.hist(data1, bins=20, density=True, color=colorname1, alpha=0.9, label=label1)
	ax.hist(data2, bins=20, density=True, color=colorname2, alpha=0.9, label=label2)
	ax.set_xlabel(r'%s' % dataname, fontsize=15)
	ax.set_ylabel(r'occurrence', fontsize=15)
	ax.legend(frameon=True, fontsize=10)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	

def triple_simple_histogram_with_label(filename, data1, data2, data3, dataname, colorname1, colorname2, colorname3, label1, label2, label3):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.hist(data1, bins=20, density=True, color=colorname1, alpha=0.9, label=label1)
	ax.hist(data2, bins=20, density=True, color=colorname2, alpha=0.9, label=label2)
	ax.hist(data3, bins=20, density=True, color=colorname3, alpha=0.9, label=label3)
	ax.set_xlabel(r'%s' % dataname, fontsize=15)
	ax.set_ylabel(r'occurrence', fontsize=15)
	ax.legend(frameon=True, fontsize=10)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 						

def double_boxplot_withlabel(filename, data1, data2, dataname, colorname1, colorname2, label1, label2):
	colors = [colorname1, colorname2]
	data = [data1, data2]

	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	bp = ax.boxplot(data, sym='',labels=[label1, label2])
	for median in bp['medians']: # change median bar
	    median.set(color='black', linewidth = 1)	
	for i in range(len(data)): # insert point inside the rectangle
		y = data[i]
		x = np.random.normal(i+1, 0.02, len(y))
		plt.plot(x, y, 'o', color=colors[i], markersize=5, markeredgewidth='0.5', markeredgecolor='black', alpha=0.5)     
	ax.set_ylabel(r'%s' % dataname, fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 

def bar_histogram(filename, x, y, x_labelname, y_labelname, colorname):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.bar(x, y, 5, color=colorname)
	ax.set_xlabel(r'%s' % x_labelname, fontsize=15)
	ax.set_ylabel(r'%s' % y_labelname, fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 		

##############
####PANELS####
##############

def simple_panel(lineage_id, x, y, x_labelname, y_labelname, colorname):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	ax.plot(x, y, 'o', color=colorname, markeredgecolor='black',markeredgewidth='1', markersize=15)
	ax.set_xlabel(r'%s' % x_labelname, fontsize=15)
	ax.set_ylabel(r'%s' % y_labelname, fontsize=15)
	return fig

def single_lineage_panel(lineage_id,mass_curve, vol_curve, area_curve, density_curve, division_time):
	nocells_in_lineage = len(mass_curve)

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	mass_cmap = cm.Reds_r
	for j in range(nocells_in_lineage):
			ax[0][0].plot(mass_curve[j][0], mass_curve[j][1], 'o', color=mass_cmap(j/nocells_in_lineage), 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
				
	for k in range(len(division_time)): 
		if k == 0:              
			ax[0][0].axvline(x=division_time[k], color='black', ls='--', lw=2, label='division')
		else :  
			ax[0][0].axvline(x=division_time[k], color='black', ls='--', lw=2)
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	ax[0][0].legend(frameon=True, fontsize=10)

	vol_cmap = cm.Blues_r
	for j in range(nocells_in_lineage):
		ax[0][1].plot(vol_curve[j][0], vol_curve[j][1], 'o', color=vol_cmap(j/nocells_in_lineage), 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
			
	for k in range(len(division_time)):                 
		ax[0][1].axvline(x=division_time[k], color='black', ls='--', lw=2)
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	#ax[0][1].legend(frameon=True, fontsize=15)

	area_cmap = cm.Greens_r
	for j in range(nocells_in_lineage):
		ax[1][0].plot(area_curve[j][0], area_curve[j][1], 'o', color=area_cmap(j/nocells_in_lineage), 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
			
	for k in range(len(division_time)):                 
		ax[1][0].axvline(x=division_time[k], color='black', ls='--', lw=2) 
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)

	density_cmap = cm.Purples_r
	for j in range(nocells_in_lineage):
		ax[1][1].plot(density_curve[j][0], density_curve[j][1], 'o', color=density_cmap(j/nocells_in_lineage), 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
			
	for k in range(len(division_time)):                 
		ax[1][1].axvline(x=division_time[k], color='black', ls='--', lw=2) 
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

	return fig

def single_lineage_panel_with_mitosis(lineage_id,mass_curve, vol_curve, area_curve, density_curve, division_time, mitotic_time):
	nocells_in_lineage = len(mass_curve)

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	mass_cmap = cm.Reds_r
	for j in range(nocells_in_lineage):
			ax[0][0].plot(mass_curve[j][0], mass_curve[j][1], 'o', color=mass_cmap(j/nocells_in_lineage), 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
				
	for k in range(len(division_time)): 
		if k == 0:              
			ax[0][0].axvline(x=division_time[k], color='black', ls='--', lw=2, label='division')
		else :  
			ax[0][0].axvline(x=division_time[k], color='black', ls='--', lw=2)
	for k in range(len(mitotic_time)):                 
		ax[0][0].axvline(x=mitotic_time[k], color='red', ls='--', lw=2)		
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	ax[0][0].legend(frameon=True, fontsize=10)

	vol_cmap = cm.Blues_r
	for j in range(nocells_in_lineage):
		ax[0][1].plot(vol_curve[j][0], vol_curve[j][1], 'o', color=vol_cmap(j/nocells_in_lineage), 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
			
	for k in range(len(division_time)):                 
		ax[0][1].axvline(x=division_time[k], color='black', ls='--', lw=2)
	for k in range(len(mitotic_time)):   
		if k == 1:              
			ax[0][1].axvline(x=mitotic_time[k], color='red', ls='--', lw=2, label='mitosis')		
		else:
			ax[0][1].axvline(x=mitotic_time[k], color='red', ls='--', lw=2)		
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	ax[0][1].legend(frameon=True, fontsize=15)

	area_cmap = cm.Greens_r
	for j in range(nocells_in_lineage):
		ax[1][0].plot(area_curve[j][0], area_curve[j][1], 'o', color=area_cmap(j/nocells_in_lineage), 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
			
	for k in range(len(division_time)):                 
		ax[1][0].axvline(x=division_time[k], color='black', ls='--', lw=2) 
	for k in range(len(mitotic_time)):                 
		ax[1][0].axvline(x=mitotic_time[k], color='red', ls='--', lw=2)		
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)

	density_cmap = cm.Purples_r
	for j in range(nocells_in_lineage):
		ax[1][1].plot(density_curve[j][0], density_curve[j][1], 'o', color=density_cmap(j/nocells_in_lineage), 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
			
	for k in range(len(division_time)):                 
		ax[1][1].axvline(x=division_time[k], color='black', ls='--', lw=2) 
	for k in range(len(mitotic_time)):                 
		ax[1][1].axvline(x=mitotic_time[k], color='red', ls='--', lw=2)		
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

	return fig

def generic_two_by_two_panel(filename, mass_curve, vol_curve, area_curve, density_curve, x_labelname, y_labelname_vector):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.45, hspace=0.4)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][1].set_xlabel(x_labelname, fontsize=15)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=15)
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(area_curve[0], area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][0].set_xlabel(x_labelname, fontsize=15)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][1].set_xlabel(x_labelname, fontsize=15)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	

def generic_two_by_two_panel_with_err(filename, mass_curve, vol_curve, area_curve, density_curve, mass_err_curve, vol_err_curve, area_err_curve, density_err_curve, x_labelname, y_labelname_vector):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.8, hspace=0.4)
	
	ax[0][0].errorbar(mass_curve[0], mass_curve[1], yerr=mass_err_curve[1], marker='o', color='red', capsize=10,
			markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][0].set_xlabel(x_labelname, fontsize=15)
	ax[0][0].set_xlabel(x_labelname, fontsize=15)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].errorbar(vol_curve[0], vol_curve[1], yerr=vol_err_curve[1], marker='o', color='blue', capsize=10,
		markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][1].set_xlabel(x_labelname, fontsize=15)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=15)
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].errorbar(area_curve[0], area_curve[1], yerr=area_err_curve[1], marker='o', color='green', capsize=10,
			markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][0].set_xlabel(x_labelname, fontsize=15)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].errorbar(density_curve[0], density_curve[1], yerr=density_err_curve[1], marker='o', color='purple', capsize=10,
		markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][1].set_xlabel(x_labelname, fontsize=15)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)		

def generic_two_by_two_panel_with_limits(filename, mass_curve, vol_curve, area_curve, density_curve, x_labelname_vector, y_labelname_vector, time_xlim, mass_ylim, vol_ylim, area_ylim, density_ylim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.65, hspace=0.4)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][0].set_xlabel(x_labelname_vector[0], fontsize=15)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=15)
	ax[0][0].set_xlim(time_xlim[0], time_xlim[1])
	ax[0][0].set_ylim(mass_ylim[0], mass_ylim[1])
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][1].set_xlabel(x_labelname_vector[1], fontsize=15)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=15)
	ax[0][1].set_xlim(time_xlim[0], time_xlim[1])
	ax[0][1].set_ylim(vol_ylim[0], vol_ylim[1])
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(area_curve[0], area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][0].set_xlabel(x_labelname_vector[2], fontsize=15)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=15)
	ax[1][0].set_xlim(time_xlim[0], time_xlim[1])
	ax[1][0].set_ylim(area_ylim[0], area_ylim[1])
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][1].set_xlabel(x_labelname_vector[3], fontsize=15)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=15)
	ax[1][1].set_xlim(time_xlim[0], time_xlim[1])
	ax[1][1].set_ylim(density_ylim[0], density_ylim[1])
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	

def generic_two_by_two_panel_with_limits_with_err(filename, mass_curve, vol_curve, area_curve, density_curve, mass_err_curve, vol_err_curve, area_err_curve, density_err_curve, 
	x_labelname, y_labelname_vector, time_xlim, mass_ylim, vol_ylim, area_ylim, density_ylim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.65, hspace=0.4)
	
	ax[0][0].errorbar(mass_curve[0], mass_curve[1], yerr=mass_err_curve[1], marker='o', color='red', capsize=10,
			markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][0].set_xlabel(x_labelname[0], fontsize=15)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=15)
	ax[0][0].set_xlim(time_xlim[0], time_xlim[1])
	ax[0][0].set_ylim(mass_ylim[0], mass_ylim[1])
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].errorbar(vol_curve[0], vol_curve[1], yerr=vol_err_curve[1], marker='o', color='blue', capsize=10,
		markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][1].set_xlabel(x_labelname[1], fontsize=15)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=15)
	ax[0][1].set_xlim(time_xlim[0], time_xlim[1])
	ax[0][1].set_ylim(vol_ylim[0], vol_ylim[1])
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].errorbar(area_curve[0], area_curve[1], yerr=area_err_curve[1], marker='o', color='green', capsize=10,
			markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][0].set_xlabel(x_labelname[2], fontsize=15)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=15)
	ax[1][0].set_xlim(time_xlim[0], time_xlim[1])
	ax[1][0].set_ylim(area_ylim[0], area_ylim[1])
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].errorbar(density_curve[0], density_curve[1], yerr=density_err_curve[1], marker='o', color='purple', capsize=10,
		markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][1].set_xlabel(x_labelname[3], fontsize=15)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=15)
	ax[1][1].set_xlim(time_xlim[0], time_xlim[1])
	ax[1][1].set_ylim(density_ylim[0], density_ylim[1])
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)			

def super_generic_two_by_two_panel_no_limits(filename, curve1, curve2, curve3, curve4, color_vector, x_labelname_vector, y_labelname_vector):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.45, hspace=0.4)
	
	ax[0][0].plot(curve1[0], curve1[1], 'o', color=color_vector[0], 
			markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][0].set_xlabel(x_labelname_vector[0], fontsize=20)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=20)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(curve2[0], curve2[1], 'o', color=color_vector[1], 
		markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][1].set_xlabel(x_labelname_vector[1], fontsize=20)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=20)
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(curve3[0], curve3[1], 'o', color=color_vector[2], 
			markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][0].set_xlabel(x_labelname_vector[2], fontsize=20)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=20)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(curve4[0], curve4[1], 'o', color=color_vector[3], 
		markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][1].set_xlabel(x_labelname_vector[3], fontsize=20)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=20)
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	 	

def super_generic_two_by_two_panel_no_limits_with_err(filename, curve1, curve2, curve3, curve4, err_curve1, err_curve2, err_curve3, err_curve4, color_vector, x_labelname_vector, y_labelname_vector):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.45, hspace=0.4)
	
	ax[0][0].errorbar(curve1[0], curve1[1], yerr=err_curve1,  marker='o', color=color_vector[0], capsize=10, 
			markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][0].set_xlabel(x_labelname_vector[0], fontsize=20)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=20)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].errorbar(curve2[0], curve2[1], yerr=err_curve2, marker='o', color=color_vector[1], capsize=10, 
		markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][1].set_xlabel(x_labelname_vector[1], fontsize=20)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=20)
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].errorbar(curve3[0], curve3[1], yerr=err_curve3, marker='o', color=color_vector[2], capsize=10, 
			markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][0].set_xlabel(x_labelname_vector[2], fontsize=20)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=20)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].errorbar(curve4[0], curve4[1], yerr=err_curve4, marker='o', color=color_vector[3], capsize=10, 
		markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][1].set_xlabel(x_labelname_vector[3], fontsize=20)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=20)
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	 	


def super_generic_two_by_two_panel_with_limits(filename, curve1, curve2, curve3, curve4, color_vector, x_labelname_vector, y_labelname_vector, curve1_ylim, curve2_ylim, curve3_ylim, curve4_ylim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.45, hspace=0.4)
	
	ax[0][0].plot(curve1[0], curve1[1], 'o', color=color_vector[0], 
			markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][0].set_xlabel(x_labelname_vector[0], fontsize=20)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=20)
	#ax[0][0].set_ylim(curve1_ylim[0], curve1_ylim[1])
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(curve2[0], curve2[1], 'o', color=color_vector[1], 
		markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][1].set_xlabel(x_labelname_vector[1], fontsize=20)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=20)
	#ax[0][1].set_ylim(curve2_ylim[0], curve2_ylim[1])
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(curve3[0], curve3[1], 'o', color=color_vector[2], 
			markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][0].set_xlabel(x_labelname_vector[2], fontsize=20)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=20)
	ax[1][0].set_ylim(curve3_ylim[0], curve3_ylim[1])
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(curve4[0], curve4[1], 'o', color=color_vector[3], 
		markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][1].set_xlabel(x_labelname_vector[3], fontsize=20)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=20)
	ax[1][1].set_ylim(curve4_ylim[0], curve4_ylim[1])
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	

def super_generic_two_by_two_panel_with_limits_with_err(filename, curve1, curve2, curve3, curve4, err_curve1, err_curve2, err_curve3, err_curve4, color_vector, x_labelname_vector, y_labelname_vector, curve1_ylim, curve2_ylim, curve3_ylim, curve4_ylim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.45, hspace=0.4)
	
	ax[0][0].errorbar(curve1[0], curve1[1], yerr=err_curve1, marker='o', color=color_vector[0], capsize=10, 
			markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][0].set_xlabel(x_labelname_vector[0], fontsize=20)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=20)
	#ax[0][0].set_ylim(curve1_ylim[0], curve1_ylim[1])
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].errorbar(curve2[0], curve2[1], yerr=err_curve2, marker='o', color=color_vector[1], capsize=10,
		markeredgecolor='black', markeredgewidth='2', markersize=20) 
	ax[0][1].set_xlabel(x_labelname_vector[1], fontsize=20)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=20)
	#ax[0][1].set_ylim(curve2_ylim[0], curve2_ylim[1])
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].errorbar(curve3[0], curve3[1], yerr=err_curve3, marker='o', color=color_vector[2], capsize=10,
			markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][0].set_xlabel(x_labelname_vector[2], fontsize=20)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=20)
	ax[1][0].set_ylim(curve3_ylim[0], curve3_ylim[1])
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].errorbar(curve4[0], curve4[1], yerr=err_curve4, marker='o', color=color_vector[3], capsize=10,
		markeredgecolor='black', markeredgewidth='2', markersize=20)
	ax[1][1].set_xlabel(x_labelname_vector[3], fontsize=20)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=20)
	ax[1][1].set_ylim(curve4_ylim[0], curve4_ylim[1])
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)		

def super_generic_two_by_two_panel_with_limits_with_err_with_colorcoding(filename, curve1, curve2, curve3, curve4, err_curve1, err_curve2, err_curve3, err_curve4, barlabel, color_vector, x_labelname_vector, y_labelname_vector, curve1_ylim, curve2_ylim, curve3_ylim, curve4_ylim, xlim):
	
	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.45, hspace=0.4)
	
	Z = [[0,0],[0,0]]
	levels = sorted(list(set(color_vector[0])))
	parcolorcoding = plt.contourf(Z, levels, cmap=cm.inferno)
	cnorm = [(c - min(color_vector[0]))/(max(color_vector[0]) - min(color_vector[0])) for c in color_vector[0]]
	ax[0][0].scatter(curve1[0], curve1[1], marker='o', s=300, color=cm.inferno(cnorm), edgecolors='black') 
	ax[0][0].errorbar(curve1[0], curve1[1], yerr=err_curve1, color='black', linestyle='None', linewidth=1, capsize=10) 
	ax[0][0].plot(curve1[0], curve1[1], '-', color='black', linewidth=1) 
	ax[0][0].set_xlabel(x_labelname_vector[0], fontsize=20)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=20)
	ax[0][0].set_xlim(xlim[0], xlim[1])
	#ax[0][0].set_ylim(curve1_ylim[0], curve1_ylim[1])

	
	Z = [[0,0],[0,0]]
	levels = sorted(list(set(color_vector[1])))
	parcolorcoding = plt.contourf(Z, levels, cmap=cm.inferno)
	cnorm = [(c - min(color_vector[1]))/(max(color_vector[1]) - min(color_vector[1])) for c in color_vector[1]]
	ax[0][1].scatter(curve2[0], curve2[1], marker='o', s=300, color=cm.inferno(cnorm), edgecolors='black') 
	ax[0][1].errorbar(curve2[0], curve2[1], yerr=err_curve2, color='black', linestyle='None', linewidth=1, capsize=10) 
	ax[0][1].plot(curve2[0], curve2[1], '-', color='black', linewidth=1) 
	ax[0][1].set_xlabel(x_labelname_vector[1], fontsize=20)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=20)
	ax[0][1].set_xlim(xlim[0], xlim[1])
	#ax[0][1].set_ylim(curve2_ylim[0], curve2_ylim[1])
	
	
	Z = [[0,0],[0,0]]
	levels = sorted(list(set(color_vector[2])))
	parcolorcoding = plt.contourf(Z, levels, cmap=cm.inferno)
	cnorm = [(c - min(color_vector[2]))/(max(color_vector[2]) - min(color_vector[2])) for c in color_vector[2]]
	ax[1][0].scatter(curve3[0], curve3[1], marker='o', s=300, color=cm.inferno(cnorm), edgecolors='black') 
	ax[1][0].errorbar(curve3[0], curve3[1], yerr=err_curve3, color='black', linestyle='None', linewidth=1, capsize=10) 
	ax[1][0].plot(curve3[0], curve3[1], '-', color='black', linewidth=1) 
	ax[1][0].set_xlabel(x_labelname_vector[2], fontsize=20)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=20)
	ax[1][0].set_xlim(xlim[0], xlim[1])
	ax[1][0].set_ylim(curve3_ylim[0], curve3_ylim[1])

	
	Z = [[0,0],[0,0]]
	levels = sorted(list(set(color_vector[3])))
	parcolorcoding = plt.contourf(Z, levels, cmap=cm.inferno)
	cnorm = [(c - min(color_vector[3]))/(max(color_vector[3]) - min(color_vector[3])) for c in color_vector[3]]
	ax[1][1].scatter(curve4[0], curve4[1], marker='o', s=300, color=cm.inferno(cnorm), edgecolors='black') 
	ax[1][1].errorbar(curve4[0], curve4[1], yerr=err_curve4, color='black', linestyle='None', linewidth=1, capsize=10) 
	ax[1][1].plot(curve4[0], curve4[1], '-', color='black', linewidth=1) 
	ax[1][1].set_xlabel(x_labelname_vector[3], fontsize=20)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=20)
	ax[1][1].set_xlim(xlim[0], xlim[1])
	ax[1][1].set_ylim(curve4_ylim[0], curve4_ylim[1])
	
	cax = plt.axes([1.05, 0.125, 0.03, 0.8])
	clb = plt.colorbar(parcolorcoding, cax=cax, format='%.1f')
	clb.set_label(barlabel, fontsize=20)
	clb.ax.tick_params(axis='both',which='minor', bottom=False, top=False, left=False, right=False)
	clb.ax.tick_params(axis='both',which='major', length=23)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)		

def single_cell_panel(lineage_id,mass_curve, vol_curve, area_curve, density_curve, x_labelname):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.4, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].set_xlabel(x_labelname, fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].set_xlabel(x_labelname, fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(area_curve[0], area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].set_xlabel(x_labelname, fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].set_xlabel(x_labelname, fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

	return fig    

def single_cell_panel_norm(lineage_id,mass_curve, vol_curve, area_curve, density_curve, x_labelname):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].set_xlabel(x_labelname, fontsize=15)
	ax[0][0].set_ylabel(r'normalized mass', fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].set_xlabel(x_labelname, fontsize=15)
	ax[0][1].set_ylabel(r'normalized volume', fontsize=15)
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(area_curve[0], area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].set_xlabel(x_labelname, fontsize=15)
	ax[1][0].set_ylabel(r'norm spreading area', fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].set_xlabel(x_labelname, fontsize=15)
	ax[1][1].set_ylabel(r'norm density', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

	return fig   	

def single_cell_panel_with_smoothed_curves(lineage_id,mass_curve, vol_curve, area_curve, density_curve, 
	mass_smoothedcurve, vol_smoothedcurve, area_smoothedcurve, density_smoothedcurve):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], '-o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='raw') 
	ax[0][0].plot(mass_smoothedcurve[0], mass_smoothedcurve[1], '-o', color='coral', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='smoothed') 
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], '-o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='raw') 
	ax[0][1].plot(vol_smoothedcurve[0], vol_smoothedcurve[1], '-o', color='lightblue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='smoothed') 	
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	ax[0][1].legend(frameon=True, fontsize=10)

	ax[1][0].plot(area_curve[0], area_curve[1], '-o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='raw')
	ax[1][0].plot(area_smoothedcurve[0], area_smoothedcurve[1], '-o', color='lime', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='smoothed')
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	ax[1][0].legend(frameon=True, fontsize=10)

	ax[1][1].plot(density_curve[0], density_curve[1], '-o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='raw')
	ax[1][1].plot(density_smoothedcurve[0], density_smoothedcurve[1], '-o', color='violet', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='smoothed')
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	ax[1][1].legend(frameon=True, fontsize=10)

	return fig    	

def single_cell_panel_with_mitosis(lineage_id,mass_curve, vol_curve, area_curve, density_curve, mitotic_time):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].axvline(x=mitotic_time, color='red', ls='--', lw=2)
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].axvline(x=mitotic_time, color='red', ls='--', lw=2, label='mitosis')
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	ax[0][1].legend(frameon=True, fontsize=10)
	
	ax[1][0].plot(area_curve[0], area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].axvline(x=mitotic_time, color='red', ls='--', lw=2)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].axvline(x=mitotic_time, color='red', ls='--', lw=2)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

def single_cell_panel_with_mitosis_and_g1s(lineage_id,mass_curve, vol_curve, area_curve, density_curve, hgem_curve, mitotic_time, g1s_transition_time):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].axvline(x=mitotic_time, color='grey', ls='--', lw=2)
	ax[0][0].axvline(x=g1s_transition_time, color='black', ls='--', lw=2)
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	left, bottom, width, height = [0.22, 0.78, 0.1, 0.1]
	ax_inset = fig.add_axes([left, bottom, width, height])
	for axis in ['top','bottom','left','right']:
  		ax_inset.spines[axis].set_linewidth(1)
	ax_inset.xaxis.set_tick_params(width=1, length= 4, which='minor', labelsize=10)
	ax_inset.xaxis.set_tick_params(width=1, length= 8, which='major',labelsize=10)
	ax_inset.yaxis.set_tick_params(width=1, length= 4, which='minor', labelsize=10)
	ax_inset.yaxis.set_tick_params(width=1, length= 8, which='major',labelsize=10)
	ax_inset.plot(hgem_curve[0], hgem_curve[1], 'o', color='black', markeredgecolor='black',markeredgewidth='0.1', markersize=5)
	ax_inset.axvline(x=g1s_transition_time, color='black', ls='--', lw=2)
	#ax_inset.set_xlabel(r'time (h)', fontsize=15)
	ax_inset.set_ylabel(r'hgemini (A.U.)', fontsize=7)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].axvline(x=mitotic_time, color='grey', ls='--', lw=2, label='mitosis')
	ax[0][1].axvline(x=g1s_transition_time, color='black', ls='--', lw=2, label='G1-S transition')
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	ax[0][1].legend(frameon=True, fontsize=10)
	
	ax[1][0].plot(area_curve[0], area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].axvline(x=mitotic_time, color='grey', ls='--', lw=2)
	ax[1][0].axvline(x=g1s_transition_time, color='black', ls='--', lw=2)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].axvline(x=mitotic_time, color='grey', ls='--', lw=2)
	ax[1][1].axvline(x=g1s_transition_time, color='black', ls='--', lw=2)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

	return fig   

def single_cell_panel_with_timepoint(lineage_id,mass_curve, vol_curve, area_curve, density_curve, timepoint, timepoint_name, xlabel):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], '-o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	if timepoint is not None:
		ax[0][0].axvline(x=timepoint, color='red', ls='--', lw=2)
	ax[0][0].set_xlabel(xlabel, fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)

	ax[0][1].plot(vol_curve[0], vol_curve[1], '-o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	if timepoint is not None:
		ax[0][1].axvline(x=timepoint, color='red', ls='--', lw=2)
	ax[0][1].set_xlabel(xlabel, fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	
	ax[1][0].plot(area_curve[0], area_curve[1], '-o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	if timepoint is not None:
		ax[1][0].axvline(x=timepoint, color='red', ls='--', lw=2, label=timepoint_name)
		ax[1][0].legend(frameon=True, fontsize=10)

	ax[1][0].set_xlabel(xlabel, fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], '-o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	if timepoint is not None:
		ax[1][1].axvline(x=timepoint, color='red', ls='--', lw=2)
	ax[1][1].set_xlabel(xlabel, fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)

def single_cell_panel_with_two_timepoint(lineage_id,mass_curve, vol_curve, area_curve, density_curve, timepoint1, timepoint_name1, timepoint2, timepoint_name2, xlabel):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.4, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], '-o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].axvline(x=timepoint1, color='grey', ls='--', lw=2)
	ax[0][0].axvline(x=timepoint2, color='black', ls='--', lw=2)
	ax[0][0].set_xlabel(xlabel, fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)

	ax[0][1].plot(vol_curve[0], vol_curve[1], '-o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].axvline(x=timepoint1, color='grey', ls='--', lw=2, label=timepoint_name1)
	ax[0][1].axvline(x=timepoint2, color='black', ls='--', lw=2, label=timepoint_name2)
	ax[0][1].set_xlabel(xlabel, fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(area_curve[0], area_curve[1], '-o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].axvline(x=timepoint1, color='grey', ls='--', lw=2)
	ax[1][0].axvline(x=timepoint2, color='black', ls='--', lw=2)
	ax[1][0].set_xlabel(xlabel, fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], '-o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].axvline(x=timepoint1, color='grey', ls='--', lw=2)
	ax[1][1].axvline(x=timepoint2, color='black', ls='--', lw=2)
	ax[1][1].set_xlabel(xlabel, fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)	 

def single_cell_panel_with_area_slope_change(lineage_id,mass_curve, vol_curve, area_curve, density_curve, area_slope_change_time, xlabel):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], '-o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].axvline(x=area_slope_change_time, color='red', ls='--', lw=2)
	ax[0][0].set_xlabel(xlabel, fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], '-o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].axvline(x=area_slope_change_time, color='red', ls='--', lw=2)
	ax[0][1].set_xlabel(xlabel, fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	
	ax[1][0].plot(area_curve[0], area_curve[1], '-o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].axvline(x=area_slope_change_time, color='red', ls='--', lw=2, label='area slope change')
	ax[1][0].set_xlabel(xlabel, fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	ax[1][0].legend(frameon=True, fontsize=10)
	
	ax[1][1].plot(density_curve[0], density_curve[1], '-o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].axvline(x=area_slope_change_time, color='red', ls='--', lw=2)
	ax[1][1].set_xlabel(xlabel, fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)	

def single_cell_panel_with_area_slope_change_and_fit(lineage_id,mass_curve, vol_curve, area_curve, density_curve, area_slope_change_time, fit_par):

	area_slope_change_frame = area_curve[0].index(area_slope_change_time)

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].axvline(x=area_slope_change_time, color='red', ls='--', lw=1)
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].axvline(x=area_slope_change_time, color='red', ls='--', lw=1)
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	
	ax[1][0].plot(area_curve[0], area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].plot(area_curve[0][:area_slope_change_frame], np.add(fit_par[1], np.multiply(fit_par[0], area_curve[0][:area_slope_change_frame])), '--', linewidth=3, color='black', label='linear fit')
	ax[1][0].axvline(x=area_slope_change_time, color='red', ls='--', lw=1, label='area slope change')
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	ax[1][0].legend(frameon=True, fontsize=10)
	
	ax[1][1].plot(density_curve[0], density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].axvline(x=area_slope_change_time, color='red', ls='--', lw=1)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)		

def single_cell_panel_with_initial_spreading_phase(lineage_id,mass_curve, vol_curve, area_curve, density_curve, 
	mass_slope_change_time, volume_slope_change_time, area_slope_change_time, density_slope_change_time, ):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cells %s' % lineage_id, fontsize=20)
	
	
	ax[0][0].plot(mass_curve[0], mass_curve[1], 'o', color='red', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	if mass_slope_change_time is not None :
		ax[0][0].axvline(x=mass_slope_change_time, color='red', ls='-', lw=2)
	if area_slope_change_time is not None :
		ax[0][0].axvline(x=area_slope_change_time, color='green', ls='--', lw=2)
	if density_slope_change_time is not None :
		ax[0][0].axvline(x=density_slope_change_time, color='purple', ls=':', lw=2)
	if volume_slope_change_time is not None :
		ax[0][0].axvline(x=volume_slope_change_time, color='blue', ls='-.', lw=2)
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve[0], vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	if mass_slope_change_time is not None :
		ax[0][1].axvline(x=mass_slope_change_time, color='red', ls='-', lw=2, label='mass slope change')	
	if area_slope_change_time is not None :
		ax[0][1].axvline(x=area_slope_change_time, color='green', ls='--', lw=2, label='area slope change')
	if density_slope_change_time is not None :
		ax[0][1].axvline(x=density_slope_change_time, color='purple', ls=':', lw=2, label='density slope change')
	if volume_slope_change_time is not None :
		ax[0][1].axvline(x=volume_slope_change_time, color='blue', ls='-.', lw=2, label='volume slope change')
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu m^{3}$)', fontsize=15)
	ax[0][1].legend(frameon=True, fontsize=10)
	
	ax[1][0].plot(area_curve[0], area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	if mass_slope_change_time is not None :
		ax[1][0].axvline(x=mass_slope_change_time, color='red', ls='-', lw=2)
	if area_slope_change_time is not None :
		ax[1][0].axvline(x=area_slope_change_time, color='green', ls='--', lw=2)
	if density_slope_change_time is not None :
		ax[1][0].axvline(x=density_slope_change_time, color='purple', ls=':', lw=2)
	if volume_slope_change_time is not None :
		ax[1][0].axvline(x=volume_slope_change_time, color='blue', ls='-.', lw=2)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu m^{2}$)', fontsize=15)
	
	ax[1][1].plot(density_curve[0], density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	if mass_slope_change_time is not None :
		ax[1][1].axvline(x=mass_slope_change_time, color='red', ls='-', lw=2)
	if area_slope_change_time is not None :
		ax[1][1].axvline(x=area_slope_change_time, color='green', ls='--', lw=2)
	if density_slope_change_time is not None :
		ax[1][1].axvline(x=density_slope_change_time, color='purple', ls=':', lw=2)
	if volume_slope_change_time is not None :
		ax[1][1].axvline(x=volume_slope_change_time, color='blue', ls='-.', lw=2)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)	

def single_cell_growth_rate_panel(lineage_id,growth_rate_mass_curve, growth_rate_vol_curve, growth_rate_area_curve, growth_rate_density_curve):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	fig.suptitle('Cell %s' % lineage_id, fontsize=20)
	
	ax[0][0].plot(growth_rate_mass_curve[0], growth_rate_mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass growth rate ($\mathrm{h}^{-1}$)', fontsize=15)
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(growth_rate_vol_curve[0], growth_rate_vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'vol growth rate ($\mathrm{h}^{-1}$)', fontsize=15)
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(growth_rate_area_curve[0], growth_rate_area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'area growth rate ($\mathrm{h}^{-1}$)', fontsize=15)
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(growth_rate_density_curve[0], growth_rate_density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density growth rate ($\mathrm{h}^{-1}$)', fontsize=15)
	#ax[1][1].legend(frameon=True, fontsize=15)

	return fig   	  

def compare_twobytwo_panel(filename,mass_curve_cond, vol_curve_cond, area_curve_cond, density_curve_cond, 
										mass_curve_ctr, vol_curve_ctr, area_curve_ctr, density_curve_ctr, x_labelname, y_labelname_vector):
	width = 4.6*2+5
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.3, hspace=0.4)
	
	ax[0][0].plot(mass_curve_cond[0], mass_curve_cond[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='peg') 
	ax[0][0].plot(mass_curve_ctr[0], mass_curve_ctr[1], '^', color='coral', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='control') 
	ax[0][0].set_xlabel(x_labelname, fontsize=15)
	ax[0][0].set_ylabel(y_labelname_vector[0], fontsize=15)
	ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(vol_curve_cond[0], vol_curve_cond[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='peg') 
	ax[0][1].plot(vol_curve_ctr[0], vol_curve_ctr[1], '^', color='deepskyblue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='control') 
	ax[0][1].set_xlabel(x_labelname, fontsize=15)
	ax[0][1].set_ylabel(y_labelname_vector[1], fontsize=15)
	ax[0][1].legend(frameon=True, fontsize=10)
	
	ax[1][0].plot(area_curve_cond[0], area_curve_cond[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='peg')
	ax[1][0].plot(area_curve_ctr[0], area_curve_ctr[1], '^', color='lime', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='control')
	ax[1][0].set_xlabel(x_labelname, fontsize=15)
	ax[1][0].set_ylabel(y_labelname_vector[2], fontsize=15)
	ax[1][0].legend(frameon=True, fontsize=10)
	
	ax[1][1].plot(density_curve_cond[0], density_curve_cond[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='peg')
	ax[1][1].plot(density_curve_ctr[0], density_curve_ctr[1], '^', color='mediumorchid', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5, label='control')
	ax[1][1].set_xlabel(x_labelname, fontsize=15)
	ax[1][1].set_ylabel(y_labelname_vector[3], fontsize=15)
	ax[1][1].legend(frameon=True, fontsize=10)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	   		

############
####MEANS###
############

def mean_phase(filename, time, meancurve, labelname, colorname, xlim):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.plot(time, meancurve, 'o', color=colorname, markeredgecolor='black',markeredgewidth='1', markersize=15)
	ax.set_xlabel(r'phase [0, 1]', fontsize=15)
	ax.set_ylabel(r'%s' % labelname, fontsize=15)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(top=1.05*meancurve[closest_index(time, xlim[1])])
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 		

def mean_time(filename, time, meancurve, x_labelname, y_labelname, colorname):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.plot(time, meancurve, 'o', color=colorname, markeredgecolor='black',markeredgewidth='1', markersize=15)
	ax.set_xlabel(r'%s' % x_labelname, fontsize=15)
	ax.set_ylabel(r'%s' % y_labelname, fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	

def mean_time_with_errors(filename, time, meancurve, standard_err, x_labelname, y_labelname, colorname):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.errorbar(time, meancurve, yerr=standard_err, color=colorname, linestyle='None', marker = 'o', capsize=5, markersize=20, mec='black', mew='1')
	ax.set_xlabel(r'%s' % x_labelname, fontsize=15)
	ax.set_ylabel(r'%s' % y_labelname, fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 		

def mean_panel_time_generic(filename, mean_curve1, mean_curve2, mean_curve3, mean_curve4, ylabel1, ylabel2, ylabel3, ylabel4, xlim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.6, hspace=0.4)
	
	ax[0][0].plot(mean_curve1[0], mean_curve1[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=15) 
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(ylabel1, fontsize=15)
	ax[0][0].set_xlim(xlim[0], xlim[1])
	ax[0][0].set_ylim(top=1.05*max(mean_curve1[1]))

	ax[0][1].plot(mean_curve2[0], mean_curve2[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=15) 
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(ylabel2, fontsize=15)
	ax[0][1].set_xlim(xlim[0], xlim[1])
	ax[0][1].set_ylim(top=1.05*max(mean_curve2[1]))
	
	ax[1][0].plot(mean_curve3[0], mean_curve3[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=15)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(ylabel3, fontsize=15)
	ax[1][0].set_xlim(xlim[0], xlim[1])
	ax[1][0].set_ylim(top=1.05*max(mean_curve3[1]))
	
	ax[1][1].plot(mean_curve4[0], mean_curve4[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=15)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(ylabel4, fontsize=15)
	ax[1][1].set_xlim(xlim[0], xlim[1])
	ax[1][1].set_ylim(top=1.05*max(mean_curve4[1]))

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	 	


def mean_panel_time(filename, mean_mass_curve, mean_vol_curve, mean_area_curve, mean_density_curve, xlim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.4, hspace=0.4)
	
	ax[0][0].plot(mean_mass_curve[0], mean_mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	ax[0][0].set_xlim(xlim[0], xlim[1])
	ax[0][0].set_ylim(top=1.05*max(mean_mass_curve[1]))
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(mean_vol_curve[0], mean_vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu \mathrm{m}^{3}$)', fontsize=15)
	ax[0][1].set_xlim(xlim[0], xlim[1])
	ax[0][1].set_ylim(top=1.05*max(mean_vol_curve[1]))
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(mean_area_curve[0], mean_area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu \mathrm{m}^{2}$)', fontsize=15)
	ax[1][0].set_xlim(xlim[0], xlim[1])
	ax[1][0].set_ylim(top=1.05*max(mean_area_curve[1]))
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(mean_density_curve[0], mean_density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu \mathrm{m}^{-3}$)', fontsize=15)
	ax[1][1].set_xlim(xlim[0], xlim[1])
	ax[1][1].set_ylim(top=1.05*max(mean_density_curve[1]))
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	 	

def mean_panel_phase(filename, mean_mass_curve, mean_vol_curve, mean_area_curve, mean_density_curve, xlim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.4, hspace=0.4)
	
	ax[0][0].plot(mean_mass_curve[0], mean_mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].set_xlabel(r'phase [0, 1]', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	ax[0][0].set_xlim(xlim[0], xlim[1])
	ax[0][0].set_ylim(top=1.05*max(mean_mass_curve[1]))
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(mean_vol_curve[0], mean_vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].set_xlabel(r'phase [0, 1]', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu \mathrm{m}^{3}$)', fontsize=15)
	ax[0][1].set_xlim(xlim[0], xlim[1])
	ax[0][1].set_ylim(top=1.05*max(mean_vol_curve[1]))
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(mean_area_curve[0], mean_area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].set_xlabel(r'phase [0, 1]', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu \mathrm{m}^{2}$)', fontsize=15)
	ax[1][0].set_xlim(xlim[0], xlim[1])
	ax[1][0].set_ylim(top=1.05*max(mean_area_curve[1]))
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(mean_density_curve[0], mean_density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].set_xlabel(r'phase [0, 1]', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu \mathrm{m}^{-3}$)', fontsize=15)
	ax[1][1].set_xlim(xlim[0], xlim[1])
	ax[1][1].set_ylim(top=1.05*max(mean_density_curve[1]))
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	 

def mean_norm_panel_time(filename, mean_mass_curve, mean_vol_curve, mean_area_curve, mean_density_curve, xlim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.4, hspace=0.4)
	
	ax[0][0].plot(mean_mass_curve[0], mean_mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'norm mass $M/M_0$', fontsize=15)
	ax[0][0].set_xlim(xlim[0], xlim[1])
	ax[0][0].set_ylim(top=1.05*max(mean_mass_curve[1]))

	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(mean_vol_curve[0], mean_vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'norm volume $V/V_0$', fontsize=15)
	ax[0][1].set_xlim(xlim[0], xlim[1])
	ax[0][1].set_ylim(top=1.05*max(mean_vol_curve[1]))
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(mean_area_curve[0], mean_area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'norm spreading area $A/A_0$', fontsize=15)
	ax[1][0].set_xlim(xlim[0], xlim[1])
	ax[1][0].set_ylim(top=1.05*max(mean_area_curve[1]))
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(mean_density_curve[0], mean_density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'norm density $\rho/\rho_0$', fontsize=15)
	ax[1][1].set_xlim(xlim[0], xlim[1])
	ax[1][1].set_ylim(top=1.05*max(mean_density_curve[1]))
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 			

def mean_norm_panel_phase(filename, mean_mass_curve, mean_vol_curve, mean_area_curve, mean_density_curve, xlim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.4, hspace=0.4)
	
	ax[0][0].plot(mean_mass_curve[0], mean_mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][0].set_xlabel(r'phase [0, 1]', fontsize=15)
	ax[0][0].set_ylabel(r'norm mass $M/M_0$', fontsize=15)
	ax[0][0].set_xlim(xlim[0], xlim[1])
	ax[0][0].set_ylim(top=1.05*max(mean_mass_curve[1]))

	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(mean_vol_curve[0], mean_vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5) 
	ax[0][1].set_xlabel(r'phase [0, 1]', fontsize=15)
	ax[0][1].set_ylabel(r'norm volume $V/V_0$', fontsize=15)
	ax[0][1].set_xlim(xlim[0], xlim[1])
	ax[0][1].set_ylim(top=1.05*max(mean_vol_curve[1]))
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(mean_area_curve[0], mean_area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][0].set_xlabel(r'phase [0, 1]', fontsize=15)
	ax[1][0].set_ylabel(r'norm spreading area $A/A_0$', fontsize=15)
	ax[1][0].set_xlim(xlim[0], xlim[1])
	ax[1][0].set_ylim(top=1.05*max(mean_area_curve[1]))
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(mean_density_curve[0], mean_density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='0.5', markersize=5)
	ax[1][1].set_xlabel(r'phase [0, 1]', fontsize=15)
	ax[1][1].set_ylabel(r'norm density $\rho/\rho_0$', fontsize=15)
	ax[1][1].set_xlim(xlim[0], xlim[1])
	ax[1][1].set_ylim(top=1.05*max(mean_density_curve[1]))
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	

def mean_panel_time_forzoom(filename, mean_mass_curve, mean_vol_curve, mean_area_curve, mean_density_curve, xlim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.4, hspace=0.4)
	
	ax[0][0].plot(mean_mass_curve[0], mean_mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='1', markersize=10) 
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'mass (pg)', fontsize=15)
	ax[0][0].set_xlim(xlim[0], xlim[1])
	ax[0][0].set_ylim(top=1.05*max([mean_mass_curve[1][i] for i in range(len(mean_mass_curve[1])) if mean_mass_curve[0][i] <= xlim[1] and mean_mass_curve[0][i] >= xlim[0]]))
	ax[0][0].set_ylim(bottom=0.95*min([mean_mass_curve[1][i] for i in range(len(mean_mass_curve[1])) if mean_mass_curve[0][i] <= xlim[1] and mean_mass_curve[0][i] >= xlim[0]]))
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(mean_vol_curve[0], mean_vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='1', markersize=10) 
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'volume ($\mu \mathrm{m}^{3}$)', fontsize=15)
	ax[0][1].set_xlim(xlim[0], xlim[1])
	ax[0][1].set_ylim(top=1.05*max([mean_vol_curve[1][i] for i in range(len(mean_vol_curve[1])) if mean_vol_curve[0][i] <= xlim[1] and mean_vol_curve[0][i] >= xlim[0]]))
	ax[0][1].set_ylim(bottom=0.95*min([mean_vol_curve[1][i] for i in range(len(mean_vol_curve[1])) if mean_vol_curve[0][i] <= xlim[1] and mean_vol_curve[0][i] >= xlim[0]]))
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(mean_area_curve[0], mean_area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='1', markersize=10)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'spreading area ($\mu \mathrm{m}^{2}$)', fontsize=15)
	ax[1][0].set_xlim(xlim[0], xlim[1])
	ax[1][0].set_ylim(top=1.05*max([mean_area_curve[1][i] for i in range(len(mean_area_curve[1])) if mean_area_curve[0][i] <= xlim[1] and mean_area_curve[0][i] >= xlim[0]]))
	ax[1][0].set_ylim(bottom=0.85*min([mean_area_curve[1][i] for i in range(len(mean_area_curve[1])) if mean_area_curve[0][i] <= xlim[1] and mean_area_curve[0][i] >= xlim[0]]))
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(mean_density_curve[0], mean_density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='1', markersize=10)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'density (pg $\mu \mathrm{m}^{-3}$)', fontsize=15)
	ax[1][1].set_xlim(xlim[0], xlim[1])
	ax[1][1].set_ylim(top=1.05*max([mean_density_curve[1][i] for i in range(len(mean_density_curve[1])) if mean_density_curve[0][i] <= xlim[1] and mean_density_curve[0][i] >= xlim[0]]))
	ax[1][1].set_ylim(bottom=0.95*min([mean_density_curve[1][i] for i in range(len(mean_density_curve[1])) if mean_density_curve[0][i] <= xlim[1] and mean_density_curve[0][i] >= xlim[0]]))
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 		


def mean_norm_panel_time_forzoom(filename, mean_mass_curve, mean_vol_curve, mean_area_curve, mean_density_curve, xlim):

	width = 4.6*2+1
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height), nrows=2, ncols=2, dpi=100)
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.92, wspace=0.4, hspace=0.4)
	
	ax[0][0].plot(mean_mass_curve[0], mean_mass_curve[1], 'o', color='red', 
			markeredgecolor='black', markeredgewidth='1', markersize=10) 
	ax[0][0].set_xlabel(r'time (h)', fontsize=15)
	ax[0][0].set_ylabel(r'norm mass $M/M_0$', fontsize=15)
	ax[0][0].set_xlim(xlim[0], xlim[1])
	ax[0][0].set_ylim(top=1.05*max([mean_mass_curve[1][i] for i in range(len(mean_mass_curve[1])) if mean_mass_curve[0][i] <= xlim[1] and mean_mass_curve[0][i] >= xlim[0]]))
	ax[0][0].set_ylim(bottom=0.95*min([mean_mass_curve[1][i] for i in range(len(mean_mass_curve[1])) if mean_mass_curve[0][i] <= xlim[1] and mean_mass_curve[0][i] >= xlim[0]]))
	#ax[0][0].legend(frameon=True, fontsize=10)

	ax[0][1].plot(mean_vol_curve[0], mean_vol_curve[1], 'o', color='blue', 
		markeredgecolor='black', markeredgewidth='1', markersize=10) 
	ax[0][1].set_xlabel(r'time (h)', fontsize=15)
	ax[0][1].set_ylabel(r'norm volume $V/V_0$', fontsize=15)
	ax[0][1].set_xlim(xlim[0], xlim[1])
	ax[0][1].set_ylim(top=1.05*max([mean_vol_curve[1][i] for i in range(len(mean_vol_curve[1])) if mean_vol_curve[0][i] <= xlim[1] and mean_vol_curve[0][i] >= xlim[0]]))
	ax[0][1].set_ylim(bottom=0.95*min([mean_vol_curve[1][i] for i in range(len(mean_vol_curve[1])) if mean_vol_curve[0][i] <= xlim[1] and mean_vol_curve[0][i] >= xlim[0]]))
	#ax[0][1].legend(frameon=True, fontsize=15)
	
	ax[1][0].plot(mean_area_curve[0], mean_area_curve[1], 'o', color='green', 
			markeredgecolor='black', markeredgewidth='1', markersize=10)
	ax[1][0].set_xlabel(r'time (h)', fontsize=15)
	ax[1][0].set_ylabel(r'norm spreading area $A/A_0$', fontsize=15)
	ax[1][0].set_xlim(xlim[0], xlim[1])
	ax[1][0].set_ylim(top=1.05*max([mean_area_curve[1][i] for i in range(len(mean_area_curve[1])) if mean_area_curve[0][i] <= xlim[1] and mean_area_curve[0][i] >= xlim[0]]))
	ax[1][0].set_ylim(bottom=0.85*min([mean_area_curve[1][i] for i in range(len(mean_area_curve[1])) if mean_area_curve[0][i] <= xlim[1] and mean_area_curve[0][i] >= xlim[0]]))
	#ax[1][0].legend(frameon=True, fontsize=15)
	
	ax[1][1].plot(mean_density_curve[0], mean_density_curve[1], 'o', color='purple', 
		markeredgecolor='black', markeredgewidth='1', markersize=10)
	ax[1][1].set_xlabel(r'time (h)', fontsize=15)
	ax[1][1].set_ylabel(r'norm density $\rho/\rho_0$', fontsize=15)
	ax[1][1].set_xlim(xlim[0], xlim[1])
	ax[1][1].set_ylim(top=1.05*max([mean_density_curve[1][i] for i in range(len(mean_density_curve[1])) if mean_density_curve[0][i] <= xlim[1] and mean_density_curve[0][i] >= xlim[0]]))
	ax[1][1].set_ylim(bottom=0.95*min([mean_density_curve[1][i] for i in range(len(mean_density_curve[1])) if mean_density_curve[0][i] <= xlim[1] and mean_density_curve[0][i] >= xlim[0]]))
	#ax[1][1].legend(frameon=True, fontsize=15)

	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 				

def mean_two_growth_rates_with_density(filename, time_mass, mean_growth_rate_mass, time_vol, mean_growth_rate_vol, time_density, mean_density, x_labelname):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.plot(time_mass, mean_growth_rate_mass, 'o', color='red', 
	markeredgecolor='black',markeredgewidth='1', markersize=15, label='mass growth rate')
	ax.plot(time_vol, mean_growth_rate_vol, 'o', color='blue', 
		markeredgecolor='black',markeredgewidth='1', markersize=15, label='volume growth rate')
	
	left, bottom, width, height = [0.4, 0.7, 0.35, 0.25]
	ax2 = fig.add_axes([left, bottom, width, height]) ## for some reason adding an inset this way does not play well with defining a width and height of the figure
	for axis in ['top','bottom','left','right']:
	  ax2.spines[axis].set_linewidth(1)
	ax2.xaxis.set_tick_params(width=1, length= 4, which='minor', labelsize=10)
	ax2.xaxis.set_tick_params(width=1, length= 8, which='major',labelsize=10)
	ax2.yaxis.set_tick_params(width=1, length= 4, which='minor', labelsize=10)
	ax2.yaxis.set_tick_params(width=1, length= 8, which='major',labelsize=10)
	ax2.plot(time_density, mean_density, 'o', color='purple', markeredgecolor='black',markeredgewidth='0.1', markersize=5)
	ax2.set_xlabel(x_labelname, fontsize=10)
	ax2.set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=10)
	ax.set_xlabel(x_labelname, fontsize=15)
	ax.set_ylabel(r'average growth rate ($\mathrm{min}^{-1}$)', fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	

def mean_three_growth_rates_with_density(filename, time_mass, mean_growth_rate_mass, time_vol, mean_growth_rate_vol, time_area, mean_growth_rate_area, time_density, mean_density, x_labelname):
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.plot(time_mass, mean_growth_rate_mass, 'o', color='red', 
	markeredgecolor='black',markeredgewidth='1', markersize=15, label='mass growth rate')
	ax.plot(time_vol, mean_growth_rate_vol, 'o', color='blue', 
		markeredgecolor='black',markeredgewidth='1', markersize=15, label='volume growth rate')
	ax.plot(time_area, mean_growth_rate_area, 'o', color='green', 
		markeredgecolor='black',markeredgewidth='1', markersize=15, label='area growth rate')
	
	left, bottom, width, height = [0.4, 0.7, 0.35, 0.25]
	ax2 = fig.add_axes([left, bottom, width, height]) ## for some reason adding an inset this way does not play well with defining a width and height of the figure
	for axis in ['top','bottom','left','right']:
	  ax2.spines[axis].set_linewidth(1)
	ax2.xaxis.set_tick_params(width=1, length= 4, which='minor', labelsize=10)
	ax2.xaxis.set_tick_params(width=1, length= 8, which='major',labelsize=10)
	ax2.yaxis.set_tick_params(width=1, length= 4, which='minor', labelsize=10)
	ax2.yaxis.set_tick_params(width=1, length= 8, which='major',labelsize=10)
	ax2.plot(time_density, mean_density, 'o', color='purple', markeredgecolor='black',markeredgewidth='0.1', markersize=5)
	ax2.set_xlabel(x_labelname, fontsize=10)
	ax2.set_ylabel(r'density (pg $\mu m^{-3}$)', fontsize=10)
	ax.set_xlabel(x_labelname, fontsize=15)
	ax.set_ylabel(r'average growth rate ($\mathrm{min}^{-1}$)', fontsize=15)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)			

#############
##ALLTRACKS##
#############

def all_single_cells(filename, cells_id_all_exp, cells_mass_all_exp, cells_vol_all_exp, cells_area_all_exp, cells_density_all_exp, x_labelname):
	nocells = len(cells_id_all_exp)
	pdf_pages = PdfPages(filename)
	for i in range(nocells):
		lineage_id = cells_id_all_exp[i]
		fig = single_cell_panel(lineage_id,cells_mass_all_exp[i], cells_vol_all_exp[i], cells_area_all_exp[i], cells_density_all_exp[i], x_labelname)
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close()

def all_cell_cycles_one_observables(filename, cells_obs_all_exp, x_labelname, y_labelname):
	pdf_pages = PdfPages(filename)
	for cell in cells_obs_all_exp:
		width = 4.6
		height = width/1.618
		fig, ax = plt.subplots(figsize=(width, height))
		ax.plot(cell[0], cell[1], '-', linewidth=1)
		ax.set_xlabel(x_labelname, fontsize=20)
		ax.set_ylabel(y_labelname, fontsize=20)
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig) 
	pdf_pages.close()
	
def all_cell_cycles_one_observables_single_plot(filename, cells_obs_all_exp, x_labelname, y_labelname):
	width = 8
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height))
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	for cell in cells_obs_all_exp:
		ax.plot(cell[0], cell[1], '-', linewidth=1)
	ax.set_xlabel(x_labelname, fontsize=20)
	ax.set_ylabel(y_labelname, fontsize=20)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 

def all_cell_cycles_one_observables_single_plot_with_lim(filename, cells_obs_all_exp, xlim, ylim, x_labelname, y_labelname):
	width = 8
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height))
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	for cell in cells_obs_all_exp:
		ax.plot(cell[0], cell[1], '-', linewidth=1)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)		
	ax.set_xlabel(x_labelname, fontsize=20)
	ax.set_ylabel(y_labelname, fontsize=20)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	

def all_cell_cycles_one_observables_single_plot_with_lim_linlog(filename, cells_obs_all_exp, xlim, ylim, x_labelname, y_labelname):
	width = 8
	height = width/1.618
	fig, ax = plt.subplots(figsize=(width, height))
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.set_yscale('log')
	for cell in cells_obs_all_exp:
		ax.plot(cell[0], cell[1], '-', linewidth=1)
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)		
	ax.set_xlabel(x_labelname, fontsize=20)
	ax.set_ylabel(y_labelname, fontsize=20)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig) 	

def all_cell_cycles_one_observable(filename, cells_id_all_exp, cells_obs_all_exp, x_labelname, y_labelname, colorname):
	nocells = len(cells_id_all_exp)
	pdf_pages = PdfPages(filename)
	for i in range(nocells):
		lineage_id = cells_id_all_exp[i]
		fig = simple_panel(lineage_id, cells_obs_all_exp[i][0], cells_obs_all_exp[i][1], x_labelname, y_labelname, colorname)
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 	

def all_single_cell_cycles(filename, cells_id_all_exp, cells_mass_all_exp, cells_vol_all_exp, cells_area_all_exp, cells_density_all_exp, x_labelname):
	nocells = len(cells_id_all_exp)
	pdf_pages = PdfPages(filename)
	for i in range(nocells):
		lineage_id = cells_id_all_exp[i]
		fig = single_cell_panel(lineage_id,cells_mass_all_exp[i], cells_vol_all_exp[i], cells_area_all_exp[i], cells_density_all_exp[i], x_labelname)
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 

def all_single_cell_cycles_norm(filename, cells_id_all_exp, cells_mass_all_exp, cells_vol_all_exp, cells_area_all_exp, cells_density_all_exp, x_labelname):
	nocells = len(cells_id_all_exp)
	pdf_pages = PdfPages(filename)
	for i in range(nocells):
		lineage_id = cells_id_all_exp[i]
		fig = single_cell_panel_norm(lineage_id,cells_mass_all_exp[i], cells_vol_all_exp[i], cells_area_all_exp[i], cells_density_all_exp[i], x_labelname)
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 	

def all_single_cell_cycles_with_timepoints(filename, cells_id_all_exp, cells_mass_all_exp, cells_vol_all_exp, cells_area_all_exp, cells_density_all_exp, timepoint, timepoint_name, xlabel):
	nocells = len(cells_id_all_exp)
	pdf_pages = PdfPages(filename)
	for i in range(nocells):
		lineage_id = cells_id_all_exp[i]
		fig = single_cell_panel_with_timepoint(lineage_id,cells_mass_all_exp[i], cells_vol_all_exp[i], cells_area_all_exp[i], cells_density_all_exp[i], timepoint[i], timepoint_name, xlabel)
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 

def all_single_cell_cycles_with_two_timepoints(filename, cells_id_all_exp, cells_mass_all_exp, cells_vol_all_exp, cells_area_all_exp, cells_density_all_exp, timepoint1, timepoint_name1, timepoint2, timepoint_name2, xlabel):
	nocells = len(cells_id_all_exp)
	pdf_pages = PdfPages(filename)
	for i in range(nocells):
		lineage_id = cells_id_all_exp[i]
		fig = single_cell_panel_with_two_timepoint(lineage_id,cells_mass_all_exp[i], cells_vol_all_exp[i], cells_area_all_exp[i], cells_density_all_exp[i], timepoint1[i], timepoint_name1, timepoint2[i], timepoint_name2, xlabel)
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 		

def all_single_cell_cycles_with_slope_change(filename, cells_id_all_exp, cells_mass_all_exp, cells_vol_all_exp, cells_area_all_exp, cells_density_all_exp, area_slope_change_time, xlabel):
	nocells = len(cells_id_all_exp)
	pdf_pages = PdfPages(filename)
	for i in range(nocells):
		lineage_id = cells_id_all_exp[i]
		fig = single_cell_panel_with_area_slope_change(lineage_id,cells_mass_all_exp[i], cells_vol_all_exp[i], cells_area_all_exp[i], cells_density_all_exp[i], area_slope_change_time[i], xlabel)
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 	

def all_single_cell_cycles_with_slope_change_and_fit(filename, cells_id_all_exp, cells_mass_all_exp, cells_vol_all_exp, cells_area_all_exp, cells_density_all_exp, area_slope_change_time, fit_par):
	nocells = len(cells_id_all_exp)
	pdf_pages = PdfPages(filename)
	for i in range(nocells):
		lineage_id = cells_id_all_exp[i]
		fig = single_cell_panel_with_area_slope_change_and_fit(lineage_id,
															cells_mass_all_exp[i], cells_vol_all_exp[i], cells_area_all_exp[i], cells_density_all_exp[i], area_slope_change_time[i], fit_par[i])
		pdf_pages.savefig(fig, bbox_inches='tight')
		plt.close(fig)	
	pdf_pages.close() 	

##############
###SCATTERS###
##############

def scatter_plot(filename, x, y, xlabel, ylabel, xlim, ylim, colorname, corr_coef):
	width=4.6
	height=width/1.618 
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.scatter(x, y, color=colorname, s=50, alpha=0.5, label='r = %.3f' % corr_coef)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	ax.legend(frameon=True, fontsize=8, loc='upper right')
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)

def scatter_plot_with_binning(filename, x, y, x_bin, y_bin, xlabel, ylabel, xlim, ylim, colorname, corr_coef):
	width=4.6
	height=width/1.618 
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.scatter(x, y, color=colorname, s=20, alpha=0.1, label='r = %.3f' % corr_coef)
	ax.plot(x_bin, y_bin, '-o', color=colorname, markeredgecolor='black')
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	ax.legend(frameon=True, fontsize=8, loc='upper right')
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)			

def binned_scatter_plot(filename, x_bin, y_bin, xlabel, ylabel, xlim, ylim, colorname):
	width=4.6
	height=width/1.618 
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.plot(x_bin, y_bin, '-o', color=colorname, markeredgecolor='black')
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)		

def binned_scatter_plot_with_error(filename, x_bin, y_bin, x_err, y_err, xlabel, ylabel, xlim, ylim, colorname):
	width=4.6
	height=width/1.618 
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.errorbar(x_bin, y_bin, xerr=x_err, yerr=y_err, color=colorname, linestyle='-', marker = 'o', capsize=5, markersize=15, mec='black', mew='1')
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)		

def double_binned_scatter_plot(filename, x1_bin, y1_bin, x2_bin, y2_bin, xlabel, ylabel, color1, color2, label1, label2, xlim, ylim):
	width=4.6
	height=width/1.618 
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.plot(x1_bin, y1_bin, '-o', color=color1, markeredgecolor='black', label=label1)
	ax.plot(x2_bin, y2_bin, '-o', color=color2, markeredgecolor='black', label=label2)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	ax.legend(frameon=True, fontsize=15)
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)		

def double_binned_scatter_plot_with_error(filename, x1_bin, y1_bin, x1_err, y1_err, x2_bin, y2_bin, x2_err, y2_err, xlabel, ylabel, color1, color2, label1, label2, xlim, ylim):
	width=4.6
	height=width/1.618 
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	ax.errorbar(x1_bin, y1_bin, xerr=x1_err, yerr=y1_err, color=color1, linestyle='None', linewidth=1, marker = 'o', capsize=5, markersize=10, mec='black', mew='1', label=label1)
	ax.errorbar(x2_bin, y2_bin, xerr=x2_err, yerr=y2_err, color=color2, linestyle='None', linewidth=1, marker = 'o', capsize=5, markersize=10, mec='black', mew='1', label=label2)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	ax.legend(frameon=True, fontsize=15)
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	

def scatter_densityplot(filename, x, y, xlabel, ylabel, legend_string):
	width=4.6
	height=width/1.618 
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	xy = np.vstack([x, y])
	z = gaussian_kde(xy)(xy) 
	znorm = [(p - z.min())/(z.max() - z.min()) for p in z]
	ax.scatter(x, y, color=cm.inferno(znorm), s=50, label=legend_string)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.legend(frameon=True, fontsize=8, loc='upper right')
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	

def scatter_densityplot_with_binning(filename, x, y, x_bin, y_bin, xlabel, ylabel, colorname, legend_string):
	width=4.6
	height=width/1.618 
	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	xy = np.vstack([x, y])
	z = gaussian_kde(xy)(xy) 
	znorm = [(p - z.min())/(z.max() - z.min()) for p in z]
	ax.plot(x_bin, y_bin, '-o', color=colorname, markeredgecolor='black', label=legend_string)
	ax.scatter(x, y, color=cm.inferno(znorm), s=50)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.legend(frameon=True, fontsize=8, loc='upper right')
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)		

def scatter_colorcoded_with_vector(filename, x, y, color_vector, xlabel, ylabel, barlabel):
	width=4.6
	height=width/1.618 

	Z = [[0,0],[0,0]]
	levels = sorted(list(set(color_vector)))
	parcolorcoding = plt.contourf(Z, levels, cmap=cm.inferno)

	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	cnorm = [(c - min(color_vector))/(max(color_vector) - min(color_vector)) for c in color_vector]
	ax.scatter(x, y, color=cm.inferno(cnorm), s=1, alpha=1)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	#ax.set_xlim(-1.5, 1.5)
	#ax.set_ylim(-10, 10)
	clb = plt.colorbar(parcolorcoding)
	clb.set_label(barlabel, fontsize=15)
	clb.ax.tick_params(axis='both',which='minor',bottom=False, top=False, left=False, right=False)
	clb.ax.tick_params(axis='both',which='major', length=8)
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	

def scatter_colorcoded_with_vector_withlim(filename, x, y, color_vector, xlim, ylim, xlabel, ylabel, barlabel):
	width=4.6
	height=width/1.618 

	Z = [[0,0],[0,0]]
	levels = sorted(list(set(color_vector)))
	parcolorcoding = plt.contourf(Z, levels, cmap=cm.inferno)

	fig, ax = plt.subplots()
	fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.99)
	#ax.set_xscale('symlog', linthreshx=0.1)
	#ax.set_yscale('symlog', linthreshx=0.1)
	cnorm = [(c - min(color_vector))/(max(color_vector) - min(color_vector)) for c in color_vector]
	ax.scatter(x, y, color=cm.inferno(cnorm), s=1, alpha=1)
	ax.set_xlabel(xlabel, fontsize=15)
	ax.set_ylabel(ylabel, fontsize=15)
	ax.set_xlim(xlim[0], xlim[1])
	ax.set_ylim(ylim[0], ylim[1])
	clb = plt.colorbar(parcolorcoding)
	clb.set_label(barlabel, fontsize=15)
	clb.ax.tick_params(axis='both',which='minor',bottom=False, top=False, left=False, right=False)
	clb.ax.tick_params(axis='both',which='major', length=8)
	fig.set_size_inches(width, height)
	fig.savefig(filename, bbox_inches='tight')  
	plt.close(fig)	
