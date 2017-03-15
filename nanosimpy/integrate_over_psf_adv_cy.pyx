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

def integrate_over_psf_adv(psf_py,track_arr_py,num_of_mol_py,psy_py,psx_py,bin_points_py):
	cdef dict psf = psf_py
	cdef dict track_arr = track_arr_py
	cdef int num_of_mol = num_of_mol_py
	cdef float psy = psy_py
	cdef float psx = psx_py
	cdef np.ndarray[DTYPE_i, ndim=1] bin_points = bin_points_py


	psf['trace'] ={}
	sys.stdout.write('\n')
	for ki in range(0, psf['number_FWHMs']):
		sys.stdout.write('\r')
		sys.stdout.write("processing FWHM %d" % psf['FWHMs'][ki])
		sys.stdout.flush()
		trace = 0
		
		for b in range(0,num_of_mol):
			
			trace += psf['V'][ki][np.round(np.sqrt((track_arr[b][1][bin_points]-psx)**2+(track_arr[b][0][bin_points]-psy)**2),0).astype(np.int32)]
			
		psf['trace'][ki] = trace
	return psf

