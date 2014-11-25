from distutils.core import setup, Extension
import numpy as np
from Cython.Distutils import build_ext

setup(
	cmdclass={'build_ext': build_ext},
	ext_modules=[Extension("co_energy",
				 sources=["co_energy.pyx", "_Energy.c"],
				 include_dirs=[np.get_include()])],
	)