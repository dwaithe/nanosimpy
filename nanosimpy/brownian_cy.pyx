import cython
cimport cython
import numpy as np
cimport numpy as np
import sys

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
DTYPE = np.int32
ctypedef np.int32_t DTYPE_i

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)


def brownian(track_arr_py,num_of_steps_py, num_of_mol_py, id_list_py,width_py,map_img_py,phob_out_py,phob_in_py,R_py,R2_py,offset_py,scale_in_py,scale_out_py):
	#Variables that have been passed.
	cdef np.ndarray[DTYPE_t, ndim=2] rand_out
	cdef np.ndarray[DTYPE_t, ndim=2] rand_in
	cdef dict track_arr = track_arr_py
	cdef int num_of_steps = num_of_steps_py
	cdef int num_of_mol = num_of_mol_py
	cdef np.ndarray[DTYPE_i, ndim=1] id_list = id_list_py
	cdef int width = width_py
	cdef np.ndarray[DTYPE_i, ndim=1] map_img = map_img_py
	cdef int R = R_py
	cdef int R2 = R2_py
	cdef float phob_out = phob_out_py
	cdef float phob_in = phob_in_py
	cdef float offset = offset_py
	cdef float scale_in = scale_in_py
	cdef float scale_out = scale_out_py
	cdef float iper
	#New variables.
	cdef int y
	cdef int x
	cdef int y1
	cdef int x1
	cdef int yl
	cdef int xl
	cdef int rc
	cdef int rcl
	cdef float phi
	cdef int part_break = 0

	#Pregenerate the random moves.
	
	sys.stdout.write('\n')
	for b in range(0, num_of_mol):
		sys.stdout.write('\r')
		per = int(np.round(100*float(b)/float(num_of_mol),0))
		sys.stdout.write("Processing tracks: [%-20s] %d%% complete" % ('='*(per/5), per))
		sys.stdout.flush()

		
		rand_in  = np.random.normal(size=[2,num_of_steps])*scale_in
		rand_out = np.random.normal(size=[2,num_of_steps])*scale_out
		

		for i in range(1, num_of_steps):

			#This is for trapping. It will not function if ptrap_off and on are 0.0
			

			if id_list[b] > 0:
				track_arr[b][:,i] = track_arr[b][:,i-1] + rand_in[:,i]
			elif id_list[b] == 0:
				track_arr[b][:,i] = track_arr[b][:,i-1] + rand_out[:,i]

			#If has moved chamber then:
			y = int(round(track_arr[b][0,i]))
			x = int(round(track_arr[b][1,i]))
			yl = int(round(track_arr[b][0,i-1]))
			xl = int(round(track_arr[b][1,i-1]))
			
			rc = y*width + x
			rcl = yl*width + xl
			
			
			#If we have hopped.
			if id_list[b] != map_img[rc]:
				phob_actual = 0
				#Jumping to the outside
				if map_img[rc] == 0:
					phob_actual = phob_out
				#Jumping to inside
				elif map_img[rc] > 0:
					phob_actual = phob_in
				#Jumping outside then inside (across). Hopefully rare.
				if (id_list[b]> 0 and map_img[rc] > 0):
					phob_actual = phob_out*phob_in
				
				#The particle has passed through the boundary in or out.
				if np.random.random()<= phob_actual:
					#Free to pass
					id_list[b] = map_img[rc]
				#The particle did not pass through boundary but it still moves withing
				else:
					#Keep trying until finds a location in existing region.
					cc = 0
					while  id_list[b] != map_img[rc] and cc<50:
						cc = cc+1
						#If located outside
						if trap_list[b] == 1 or id_list[b] != 0:
							r = np.random.normal(size=[2])*scale_in
						elif id_list[b] ==0 or trap_list[b] ==0:
							r = np.random.normal(size=[2])*scale_out

						#Update
						track_arr[b][:,i] = track_arr[b][:,i-1] + r
						y = int(round(track_arr[b][0,i]))
						x = int(round(track_arr[b][1,i]))
						rc = y*width + x
					#If we try 50x then isn't going anywhere, best just keep in present location
					if cc ==50:
						pass
						part_break += 1
						
						
			#If our particle has moved there is a chance it has crossed the simulation boundary.    
			if (float(y)-offset)**2 + (x-offset)**2>= R2:
				
					#We calculate the opposing point on the simulation and reposition.
					phi = np.arctan2(float(y)-float(offset), float(x)-float(offset));
					x1 = int(round(((x-offset) -  2.0*(R) * np.cos(phi)+offset)));
					y1 = int(round(((y-offset) -  2.0*(R) * np.sin(phi)+offset)));
					#Update
					track_arr[b][0,i] = float(y1)
					track_arr[b][1,i] = float(x1)
					rc = y1*width + x1
					id_list[b] = map_img[rc]
			else:
				id_list[b] = map_img[rc]

	print 'Particle bounded 50 times on ',part_break,' occassions'
	return track_arr
