"""  Wraps _Energy.c for use in python.  """

# Declare the prototype of the C function we are interested in calling.
cdef extern from "_Energy.h":
	double _Energy(int v_i, int J_i, double massC, double massO)

# Import the Python- and C-level symbols of numpy
# import numpy as np
# cimport numpy as np

# Numpy must be initialized to avoid segfaults
# np.import_array()

# create the wrapper code
def co_energy(v_i, J_i, massC, massO):
	# return _Energy(int v_i, int J_i, double massC, double massO)
	return _Energy(v_i, J_i, massC, massO)
