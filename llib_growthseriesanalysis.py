# -*- coding: utf-8 -*-


import re
import llib_customfunction
import statistics
import numpy as np
#list(dict.fromkeys(keywords)

##################

class CellDataFrame():
	def __init__(self, s, dt):
		cells = [] #cells[id][frame][time, mass, volume, area]
		cell = []
		id_original = []
		with open(s, 'r+') as file: 
			line_index = 0
			for line in file: 
				data_point = line.split()
				#data_point[0] = re.sub('\D','',data_point[0]) # this trick removes anything that is not a number from the 1st col
				#cell_id = int(data_point[0])
				cell_id = data_point[0]
				if line_index == 0:
					prev_cell_id = data_point[0]
				if cell_id != prev_cell_id:
					id_original.append(prev_cell_id)
					cells.append(cell)
					cell = []
					prev_cell_id = cell_id
					cell_instant_vec = [float(data_point[i]) if data_point[i] != 'None' else None for i in range(1, len(data_point))]
					cell.append(cell_instant_vec)
				else :    
					cell_instant_vec = [float(data_point[i]) if data_point[i] != 'None' else None for i in range(1, len(data_point))]
					cell.append(cell_instant_vec)
				line_index += 1   
			id_original.append(prev_cell_id)     
		cells.append(cell)
		
		max_no_frame = 0
		for nocells in range(len(cells)):
			if len(cells[nocells]) > max_no_frame:
				max_no_frame = len(cells[nocells])       
		time = []
		for n in range(max_no_frame):
			time.append(n*dt)    
	
		self.Data = cells
		self.NoCell = len(cells)
		self.Time = time
		self.TimeStep = dt
		self.MaxNoFrame = max_no_frame 
		self.CellsOriginalId = id_original   

	def single_cell_mitotictime_within_a_lineage(self, nocell):
		mitotic_time_vector = []

		if nocell < len(self.Data):

			cell_timevol = [self.TimeStep*self.Data[nocell][i][0]          
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][2] is not None and self.Data[nocell][i][0] is not None]
			cell_vol = [self.Data[nocell][i][2] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][2] is not None and self.Data[nocell][i][0] is not None]
			cell_linid = [self.Data[nocell][i][5] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][5] is not None and self.Data[nocell][i][0] is not None]         
		  
			linid_vec = list(dict.fromkeys(cell_linid)) # we use dictionaries, not sets to preserve order
			for linid in linid_vec:
				time_linid = []
				vol_linid = []
				for i in range(len(cell_timevol)):
					if cell_linid[i] == linid:
						time_linid.append(cell_timevol[i])
						vol_linid.append(cell_vol[i])

				if len(time_linid) > 2:		
					vol_linid_increm = [ vol_linid[i+2] - vol_linid[i] for i in range(len(time_linid)-2)] 		
					time_linid_increm = [time_linid[i+2] - time_linid[i] for i in range(len(time_linid)-2)]
					der_vol_linid =  list(np.divide(vol_linid_increm, time_linid_increm)) 
					max_der = max(der_vol_linid)
					if max_der > 600:
						mitotic_time_vector.append(time_linid[der_vol_linid.index(max_der)])
			return mitotic_time_vector 
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)           

	def single_cell_divisiontime_within_a_lineage(self, nocell):
		division_time_vector = []

		if nocell < len(self.Data):

			cell_timemass = [self.TimeStep*self.Data[nocell][i][0]          
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][1] is not None and self.Data[nocell][i][0] is not None]
			cell_mass = [self.Data[nocell][i][1] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][1] is not None and self.Data[nocell][i][0] is not None]
			cell_linid = [self.Data[nocell][i][5] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][5] is not None and self.Data[nocell][i][0] is not None]         

			max_time = max(cell_timemass) # this trick is to avoid counting as division time the last experimental point           
			linid_vec = list(dict.fromkeys(cell_linid)) # we use dictionaries, not sets to preserve order
			for linid in linid_vec:
				time_linid = []
				mass_linid = []
				for i in range(len(cell_timemass)):
					if cell_linid[i] == linid:
						time_linid.append(cell_timemass[i])
						mass_linid.append(cell_mass[i])
				if linid == linid_vec[0]:       #this singles out the very first cell of the lineage
					if time_linid[len(time_linid)-1] < max_time:        
						division_time_vector.append(time_linid[len(time_linid)-1])  
				else : #condition excludes the last point
					if time_linid[len(time_linid)-1] < max_time and mass_linid[len(time_linid)-1]/mass_linid[0] > 1.4:        
						division_time_vector.append(time_linid[len(time_linid)-1])       
			return division_time_vector 
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)         

	def single_cell_G1Stime_within_a_lineage(self, nocell):

		frame_increment = 4
		transition_threshold = 0.3

		G1S_transition_vector = []
		if nocell < len(self.Data):

			cell_timehgem = [self.TimeStep*self.Data[nocell][i][0]          
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][4] is not None]
			cell_hgem = [self.Data[nocell][i][4] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][4] is not None]
			cell_linid = [self.Data[nocell][i][5] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][5] is not None]

			linid_vec = list(dict.fromkeys(cell_linid))       
			for linid in linid_vec:
				time_linid = []
				hgem_linid = []
				for i in range(len(cell_timehgem)):
					if cell_linid[i] == linid:
						time_linid.append(cell_timehgem[i])
						hgem_linid.append(cell_hgem[i]) 
				hgem_increments = [(hgem_linid[i+frame_increment] - hgem_linid[i])/hgem_linid[i] for i in range(len(hgem_linid)-frame_increment)] 
				
				if len(hgem_increments) > 0 : # this one checks that the increment vector if well defined
					max_hgem_increment = max(hgem_increments)
					if max_hgem_increment > transition_threshold and hgem_increments.index(max_hgem_increment) > 3: 
						# this one check that the maximum makes sense and is not too low and not too close to the beginning
						G1S_transition_vector.append(time_linid[hgem_increments.index(max_hgem_increment)])  
					else:
						G1S_transition_vector.append('None')  
				else:
					G1S_transition_vector.append('None')        

			return G1S_transition_vector  
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)

	def single_cell_hgem_curves_within_a_lineage(self, nocell):
		hgem_curves_within_a_lineage = []

		if nocell < len(self.Data):

			cell_timehgem = [self.TimeStep*self.Data[nocell][i][0]          
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][4] is not None]
			cell_hgem = [self.Data[nocell][i][4] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][4] is not None]
			cell_linid = [self.Data[nocell][i][5] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][5] is not None]

			linid_vec = list(dict.fromkeys(cell_linid))         
			for linid in linid_vec:
				time_linid = []
				hgem_linid = []
				for i in range(len(cell_timehgem)):
					if cell_linid[i] == linid:
						time_linid.append(cell_timehgem[i])
						hgem_linid.append(cell_hgem[i])    
				hgem_curves_within_a_lineage.append([time_linid, hgem_linid])

			return hgem_curves_within_a_lineage  
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)          
		
						  
	def single_cell_masscurve(self, nocell):
		if nocell < len(self.Data):

			cell_timemass = [self.TimeStep*self.Data[nocell][i][0]          
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][1] is not None and self.Data[nocell][i][0] is not None]
			cell_mass = [self.Data[nocell][i][1] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][1] is not None and self.Data[nocell][i][0] is not None]       
			
			return([cell_timemass, cell_mass])
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)     

	def single_cell_masscurves_within_a_lineage(self, nocell): 
		mass_curves_within_a_lineage = []    

		if nocell < len(self.Data):

			cell_timemass = [self.TimeStep*self.Data[nocell][i][0]          
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][1] is not None and self.Data[nocell][i][0] is not None]
			cell_mass = [self.Data[nocell][i][1] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][1] is not None and self.Data[nocell][i][0] is not None]
			cell_linid = [self.Data[nocell][i][5] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][5] is not None and self.Data[nocell][i][0] is not None]
					
			linid_vec = list(dict.fromkeys(cell_linid))         
			for linid in linid_vec:
				time_linid = []
				mass_linid = []
				for i in range(len(cell_timemass)):
					if cell_linid[i] == linid:
						time_linid.append(cell_timemass[i])
						mass_linid.append(cell_mass[i]) 
				mass_curves_within_a_lineage.append([time_linid, mass_linid])     

			return mass_curves_within_a_lineage    
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)   


	def single_cell_volcurve(self, nocell):
		if nocell < len(self.Data):
			cell_timevol = [self.TimeStep*self.Data[nocell][i][0]           
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][2] is not None and self.Data[nocell][i][0] is not None]        
			cell_vol = [self.Data[nocell][i][2] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][2] is not None and self.Data[nocell][i][0] is not None]
			
			return([cell_timevol, cell_vol])
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)

	def single_cell_volcurves_within_a_lineage(self, nocell):
		vol_curves_within_a_lineage = []

		if nocell < len(self.Data):

			cell_timevol = [self.TimeStep*self.Data[nocell][i][0]          
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][2] is not None and self.Data[nocell][i][0] is not None]
			cell_vol = [self.Data[nocell][i][2] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][2] is not None and self.Data[nocell][i][0] is not None]
			cell_linid = [self.Data[nocell][i][5] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][5] is not None and self.Data[nocell][i][0] is not None]

			linid_vec = list(dict.fromkeys(cell_linid))         
			for linid in linid_vec:
				time_linid = []
				vol_linid = []
				for i in range(len(cell_timevol)):
					if cell_linid[i] == linid:
						time_linid.append(cell_timevol[i])
						vol_linid.append(cell_vol[i]) 
				vol_curves_within_a_lineage.append([time_linid, vol_linid])

			return vol_curves_within_a_lineage    
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)  


	def single_cell_areacurve(self, nocell):
		if nocell < len(self.Data):
			cell_timearea = [self.TimeStep*self.Data[nocell][i][0]            
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][3] is not None and self.Data[nocell][i][0] is not None]        
			cell_area = [self.Data[nocell][i][3] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][3] is not None and self.Data[nocell][i][0] is not None]
			
			return([cell_timearea, cell_area])
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells) 

	def single_cell_areacurves_within_a_lineage(self, nocell):
		area_curves_within_a_lineage = []

		if nocell < len(self.Data):

			cell_timearea = [self.TimeStep*self.Data[nocell][i][0]          
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][3] is not None and self.Data[nocell][i][0] is not None]
			cell_area = [self.Data[nocell][i][3] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][3] is not None and self.Data[nocell][i][0] is not None]
			cell_linid = [self.Data[nocell][i][5] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][5] is not None and self.Data[nocell][i][0] is not None]

			linid_vec = list(dict.fromkeys(cell_linid))         
			for linid in linid_vec:
				time_linid = []
				area_linid = []
				for i in range(len(cell_timearea)):
					if cell_linid[i] == linid:
						time_linid.append(cell_timearea[i])
						area_linid.append(cell_area[i]) 
				area_curves_within_a_lineage.append([time_linid, area_linid])

			return area_curves_within_a_lineage    
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells) 
						
	
	def single_cell_densitycurve(self, nocell):
		if nocell < len(self.Data):
			cell_timedensity = [self.TimeStep*self.Data[nocell][i][0]           
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][1] is not None and self.Data[nocell][i][2] is not None and self.Data[nocell][i][0] is not None]
			cell_density = [self.Data[nocell][i][1]/self.Data[nocell][i][2] for i in range(len(self.Data[nocell]))       
					if self.Data[nocell][i][1] is not None and self.Data[nocell][i][2] is not None and self.Data[nocell][i][0] is not None]
			
			return([cell_timedensity, cell_density])
		else :
			return('Cell ID must be between 0 and %g' % len(self.Data))    

	def single_cell_densitycurves_within_a_lineage(self, nocell): 
		density_curves_within_a_lineage = []

		if nocell < len(self.Data):

			cell_timedensity = [self.TimeStep*self.Data[nocell][i][0]           
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][1] is not None and self.Data[nocell][i][2] is not None]
			cell_density = [self.Data[nocell][i][1]/self.Data[nocell][i][2] for i in range(len(self.Data[nocell]))       
					if self.Data[nocell][i][1] is not None and self.Data[nocell][i][2] is not None]
			cell_linid = [self.Data[nocell][i][5] for i in range(len(self.Data[nocell])) 
					if self.Data[nocell][i][5] is not None and self.Data[nocell][i][0] is not None]

			linid_vec = list(dict.fromkeys(cell_linid))         
			for linid in linid_vec:
				time_linid = []
				density_linid = []
				for i in range(len(cell_timedensity)):
					if cell_linid[i] == linid:
						time_linid.append(cell_timedensity[i])
						density_linid.append(cell_density[i]) 
				density_curves_within_a_lineage.append([time_linid, density_linid])

			return density_curves_within_a_lineage    
		else :
			return('Cell ID must be between 0 and %g' % self.NoCells)      

#################################################################################################                  

	def single_cell_svratio(self, nocell):
		if nocell < len(self.Data):
			cell_timesvratio = [self.Data[nocell][i][0]            
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][3] is not None and self.Data[nocell][i][2] is not None]
			cell_svratio = [self.Data[nocell][i][3]/self.Data[nocell][i][2] for i in range(len(self.Data[nocell]))       
					if self.Data[nocell][i][3] is not None and self.Data[nocell][i][2] is not None]
			
			return([cell_timesvratio, cell_svratio])
		else :
			return('Cell ID must be between 0 and %g' % len(self.Data))   
			

	def single_cell_smratio(self, nocell):
		if nocell < len(self.Data):
			cell_timesmratio = [self.Data[nocell][i][0]            
			for i in range(len(self.Data[nocell]))
			if self.Data[nocell][i][3] is not None and self.Data[nocell][i][1] is not None]
			cell_smratio = [self.Data[nocell][i][3]/self.Data[nocell][i][1] for i in range(len(self.Data[nocell]))       
					if self.Data[nocell][i][3] is not None and self.Data[nocell][i][1] is not None]
			
			return([cell_timesmratio, cell_smratio])
		else :
			return('Cell ID must be between 0 and %g' % len(self.Data))               
			
	def allmasses(self):
		mass = [self.Data[nocell][i][1] 
						for nocell in range(self.NoCell) for i in range(len(self.Data[nocell])) 
						if self.Data[nocell][i][1] is not None]
		return mass
	
	def allvolumes(self):
		vol = [self.Data[nocell][i][2]
						for nocell in range(self.NoCell) for i in range(len(self.Data[nocell])) 
						if self.Data[nocell][i][2] is not None]
		return vol
	
	def alldensities(self):
		density = [self.Data[nocell][i][1]/self.Data[nocell][i][2]
						for nocell in range(self.NoCell) for i in range(len(self.Data[nocell])) 
						if self.Data[nocell][i][1] is not None and self.Data[nocell][i][2] is not None]
		return density
	
	def allmassderivatives(self, tau):
		massder = [self.Data[nocell][i+tau][1] - self.Data[nocell][i][1]
						for nocell in range(self.NoCell) for i in range(len(self.Data[nocell])-tau) 
						if self.Data[nocell][i+tau][1] is not None and self.Data[nocell][i][1] is not None]
		return massder
	
	def allvolderivatives(self, tau):
		volder = [self.Data[nocell][i+tau][2] - self.Data[nocell][i][2]
						for nocell in range(self.NoCell) for i in range(len(self.Data[nocell])-tau) 
						if self.Data[nocell][i+tau][2] is not None and self.Data[nocell][i][2] is not None]
		return volder
							
	def ensembleavgstd_cellmass_time_series(self): 
		dt = self.TimeStep                
		meanmass = []
		dstdmass = []              
		time = []      
		for frame in range(self.MaxNoFrame):
			framemass = []
			for nocell in range(self.NoCell):
				if frame < len(self.Data[nocell]) :
					if self.Data[nocell][frame][1] is not None:
						framemass.append(self.Data[nocell][frame][1])
			if len(framemass) > 1 :            
				meanm = statistics.mean(framemass)
				dvstdm = statistics.stdev(framemass)             
				meanmass.append(meanm)  
				dstdmass.append(dvstdm)
				time.append(frame*dt)
			
		return(time, meanmass, dstdmass)  
		
	def ensembleavgstd_cellvol_time_series(self): 
		dt = self.TimeStep                
		meanvol = []
		dstdvol = []              
		time = []      
		for frame in range(self.MaxNoFrame):
			framevol = []
			for nocell in range(self.NoCell):
				if frame < len(self.Data[nocell]) :
					if self.Data[nocell][frame][2] is not None:
						framevol.append(self.Data[nocell][frame][2])
			if len(framevol) > 1 :            
				meanv = statistics.mean(framevol)            
				dvstdv = statistics.stdev(framevol)             
				meanvol.append(meanv)  
				dstdvol.append(dvstdv)
				time.append(frame*dt)
			
		return(time, meanvol, dstdvol)     
		
	def ensembleavgstd_celldensity_time_series(self):   
		dt = self.TimeStep                
		meandensity = []
		dstddensity = []  
		time = []                  
		for frame in range(self.MaxNoFrame):
			framedensity = []
			for nocell in range(self.NoCell):
				if frame < len(self.Data[nocell]) :
					if self.Data[nocell][frame][1] is not None and self.Data[nocell][frame][2] is not None:
						framedensity.append(self.Data[nocell][frame][1]/self.Data[nocell][frame][2])
			if len(framedensity) > 1 :  
				mean = statistics.mean(framedensity)            
				dvstd = statistics.stdev(framedensity) 
				meandensity.append(mean)  
				dstddensity.append(dvstd)
				time.append(frame*dt)
			
		return(time, meandensity, dstddensity)     

	def ensembleavgstd_massgrowthrate_timeseries(self, tau): 
		dt = self.TimeStep                
		meanmass = []
		dstdmass = []             
		time = []      
		for frame in range(self.MaxNoFrame):
			framemass = []
			for nocell in range(self.NoCell):
				if frame < len(self.Data[nocell])-tau :
					if self.Data[nocell][frame+tau][1] is not None and self.Data[nocell][frame][1] is not None and self.Data[nocell][frame+tau][2] is not None and self.Data[nocell][frame][2] is not None:
						mass_growthrate = (self.Data[nocell][frame+tau][1]-self.Data[nocell][frame][1])/(dt*tau*self.Data[nocell][frame][1])
						framemass.append(mass_growthrate)
			if len(framemass) > 1 :  
				meanm = statistics.mean(framemass)
				dvstdm = statistics.stdev(framemass)             
				meanmass.append(meanm)  
				dstdmass.append(dvstdm)
				time.append(frame*dt)
			
		return(time, meanmass, dstdmass) 
		
	def ensembleavgstd_volgrowthrate_timeseries(self, tau): 
		dt = self.TimeStep                
		meanvol = []
		dstdvol = []             
		time = []      
		for frame in range(self.MaxNoFrame):
			framevol = []
			for nocell in range(self.Noself.Data):
				if frame < len(self.Data[nocell])-tau :
					if self.Data[nocell][frame+tau][1] is not None and self.Data[nocell][frame][1] is not None and self.Data[nocell][frame+tau][2] is not None and self.Data[nocell][frame][2] is not None:
						vol_growthrate = (self.Data[nocell][frame+tau][2]-self.Data[nocell][frame][2])/(dt*tau*self.Data[nocell][frame][2])
						framevol.append(vol_growthrate)
			if len(framevol) > 1 :  
				meanv = statistics.mean(framevol)
				dvstdv = statistics.stdev(framevol)             
				meanvol.append(meanv)  
				dstdvol.append(dvstdv)
				time.append(frame*dt)
			
		return(time, meanvol, dstdvol)      
		
	def ensembleavgstd_massgrowthspeed_vs_mass(self, tau, b):
		dt = self.TimeStep 
		massderivative = []
		mass = []
		
		for nocell in range(self.NoCell):
			for frame in range(len(self.Data[nocell])-tau):
				if self.Data[nocell][frame+tau][1] is not None and self.Data[nocell][frame][1] is not None :
					massderivative.append((self.Data[nocell][frame+tau][1] - self.Data[nocell][frame][1])/(dt*tau))
					mass.append(self.Data[nocell][frame][1])               
		masscurve = customfunction.binningdata_median(mass, massderivative, b)   
		
		return(masscurve) #(xx, yy, dev_xx, dev_yy)
		
	def ensembleavgstd_volgrowthspeed_vs_vol(self, tau, b):
		dt = self.TimeStep 
		volderivative = []
		vol = []
		
		for nocell in range(self.NoCell):
			for frame in range(len(self.Data[nocell])-tau):
				if self.Data[nocell][frame+tau][2] is not None and self.Data[nocell][frame][2] is not None :
					volderivative.append((self.Data[nocell][frame+tau][2] - self.Data[nocell][frame][2])/(dt*tau))
					vol.append(self.Data[nocell][frame][2])               
		volcurve = customfunction.binningdata_median(vol, volderivative, b)   
		
		return(volcurve)  #(xx, yy, dev_xx, dev_yy) 
		
	def ensembleavgstd_massgrowthspeed_vs_density(self, tau, b):
		dt = self.TimeStep 
		massderivative = []
		density = []
		
		for nocell in range(self.NoCell):
			for frame in range(len(self.Data[nocell])-tau):
				if self.Data[nocell][frame+tau][1] is not None and self.Data[nocell][frame][1] is not None and self.Data[nocell][frame][2] is not None :
					massderivative.append((self.Data[nocell][frame+tau][1] - self.Data[nocell][frame][1])/(dt*tau))
					density.append(self.Data[nocell][frame][1]/self.Data[nocell][frame][2])               
		masscurve = customfunction.binningdata_median(density, massderivative, b)   
		
		return(masscurve)
		
	def ensembleavgstd_volgrowthspeed_vs_density(self, tau, b):
		dt = self.TimeStep 
		volderivative = []
		density = []
		
		for nocell in range(self.NoCell):
			for frame in range(len(self.Data[nocell])-tau):
				if self.Data[nocell][frame+tau][2] is not None and self.Data[nocell][frame][2] is not None and self.Data[nocell][frame][1] is not None:
					volderivative.append((self.Data[nocell][frame+tau][2] - self.Data[nocell][frame][2])/(dt*tau))
					density.append(self.Data[nocell][frame][1]/self.Data[nocell][frame][2])               
		volcurve = customfunction.binningdata_median(density, volderivative, b)   
		
		return(volcurve)    
		
	def ensembleavgstd_massgrowthrate_vs_density(self, tau, b):
		dt = self.TimeStep 
		masslogderivative = []
		density = []
		
		for nocell in range(self.NoCell):
			for frame in range(len(self.Data[nocell])-tau):
				if self.Data[nocell][frame+tau][1] is not None and self.Data[nocell][frame][1] is not None and self.Data[nocell][frame][2] is not None :
					masslogderivative.append((self.Data[nocell][frame+tau][1] - self.Data[nocell][frame][1])/(dt*tau*self.Data[nocell][frame][1]))
					density.append(self.Data[nocell][frame][1]/self.Data[nocell][frame][2])               
		masscurve = customfunction.binningdata_median(density, masslogderivative, b)   
		
		return(masscurve)
		
	def ensembleavgstd_volgrowthrate_vs_density(self, tau, b):
		dt = self.TimeStep 
		vollogderivative = []
		density = []
		
		for nocell in range(self.NoCell):
			for frame in range(len(self.Data[nocell])-tau):
				if self.Data[nocell][frame+tau][2] is not None and self.Data[nocell][frame][2] is not None and self.Data[nocell][frame][1] is not None:
					vollogderivative.append((self.Data[nocell][frame+tau][2] - self.Data[nocell][frame][2])/(dt*tau*self.Data[nocell][frame][2]))
					density.append(self.Data[nocell][frame][1]/self.Data[nocell][frame][2])               
		volcurve = customfunction.binningdata_median(density, vollogderivative, b)   
		
		return(volcurve)       