# -*- coding: utf-8 -*-


import statistics
from math import sqrt

########################
###PRINTING FUNCTIONS###
########################

def print_to_file(filename, column_vector):
    no_of_columns = len(column_vector)
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

#######################
###BINNING FUNCTIONS###
#######################

def scatter_to_func(x, y):
    ind = sorted(range(len(x)), key=lambda k: x[k])
    yy = []
    for k in range(len(y)):
        yy.append(y[ind[k]])
    return (sorted(x), yy)

def binningdata_mean(x, y, b): #this one bins data in intervals containing b datapoints     
    ordered_data = scatter_to_func(x, y)    
    x = ordered_data[0]
    y = ordered_data[1]
                                  
    n = len(x)
    xx = []
    dev_xx = []
    yy = []
    dev_yy = []
    
    k = 0
    while k*b+b < n :
        xv = []
        yv = []  
        for i in range(k*b, k*b+b) :
            xv.append(x[i])
            yv.append(y[i]) 
        xx.append(statistics.mean(xv))
        yy.append(statistics.mean(yv))
        if b > 1 :
            dev_xx.append(statistics.stdev(xv))
            dev_yy.append(statistics.stdev(yv))
        else :
            dev_xx.append(0)
            dev_yy.append(0)
        k += 1    
        
    return (xx, yy, dev_xx, dev_yy)

def binningdata_median(x, y, b): #this one bins data in intervals containing b datapoints     
    ordered_data = scatter_to_func(x, y)    
    x = ordered_data[0]
    y = ordered_data[1]
                                  
    n = len(x)
    xx = []
    dev_xx = []
    yy = []
    dev_yy = []
    
    k = 0
    while k*b+b < n :
        xv = []
        yv = []  
        for i in range(k*b, k*b+b) :
            xv.append(x[i])
            yv.append(y[i]) 
        xx.append(statistics.median(xv))
        yy.append(statistics.median(yv))
        if b > 1 :
            dev_xx.append(statistics.stdev(xv))
            dev_yy.append(statistics.stdev(yv))
        else :
            dev_xx.append(0)
            dev_yy.append(0)
        k += 1    
        
    return (xx, yy, dev_xx, dev_yy)

def binningdata_mode(x, y, b): #this one bins data in intervals containing b datapoints     
    ordered_data = scatter_to_func(x, y)    
    x = ordered_data[0]
    y = ordered_data[1]
                                  
    n = len(x)
    xx = []
    dev_xx = []
    yy = []
    dev_yy = []
    
    k = 0
    while k*b+b < n :
        xv = []
        yv = []  
        for i in range(k*b, k*b+b) :
            xv.append(x[i])
            yv.append(y[i]) 
        print(xv)    
        xx.append(statistics.mode(xv))
        yy.append(statistics.mode(yv))
        if b > 1 :
            dev_xx.append(statistics.stdev(xv))
            dev_yy.append(statistics.stdev(yv))
        else :
            dev_xx.append(0)
            dev_yy.append(0)
        k += 1    
        
    return (xx, yy, dev_xx, dev_yy)

#######################
###BINNING FUNCTIONS###
#######################

def mothers(cells_mass, cells_vol, cells_area, cells_density, cells_id):
    cells_mass_mothers, cells_vol_mothers, cells_area_mothers, cells_density_mothers, cells_id_mothers = [], [], [], [], []
    for i in range(len(cells_id)) :
        if cells_id[i][len(cells_id[i])-1] == '1' : ## we built the id names so that if mothers, the last letter is the number 1
            cells_mass_mothers.append(cells_mass[i])
            cells_vol_mothers.append(cells_vol[i])
            cells_area_mothers.append(cells_area[i])
            cells_density_mothers.append(cells_density[i])
            cells_id_mothers.append(cells_id[i])
    return cells_mass_mothers, cells_vol_mothers, cells_area_mothers, cells_density_mothers, cells_id_mothers       

def read_plateaux_file(plateaux_file):
    plateaux_times = []
    with open(plateaux_file, 'r+') as file: 
                line_index = 0
                for line in file: 
                    data_point = line.split()
                    if data_point[0] != 'None':
                        plateaux_times.append([float(data_point[0]), float(data_point[1])])
                    else:
                        plateaux_times.append(None) 
    return plateaux_times           

def scatter_without_plateaux_in_time_interval(plateaux_times, cells_obs1, cells_obs2, time_interval):
    tmin, tmax = time_interval[0], time_interval[1]
    no_cells, scatter_plot = len(cells_obs1), [[], []]  
    for n in range(no_cells):
        for t in range(len(cells_obs1[n][0])):

                if plateaux_times[n] != None:
                    if (cells_obs1[n][0][t] < plateaux_times[n][0] or cells_obs1[n][0][t] > plateaux_times[n][1]) and cells_obs1[n][0][t] >= tmin and cells_obs1[n][0][t] <= tmax:
                        for k in range(len(cells_obs2[n][0])):
                            if cells_obs1[n][0][t] == cells_obs2[n][0][k] :
                                scatter_plot[0].append(cells_obs1[n][1][t])
                                scatter_plot[1].append(cells_obs2[n][1][k])
                else:       
                    if cells_obs1[n][0][t] >= tmin and cells_obs1[n][0][t] <= tmax :
                        for k in range(len(cells_obs2[n][0])):
                                if cells_obs1[n][0][t] == cells_obs2[n][0][k] :
                                    scatter_plot[0].append(cells_obs1[n][1][t]) 
                                    scatter_plot[1].append(cells_obs2[n][1][k])
    return scatter_plot 