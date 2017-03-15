from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
	include_dirs=[numpy.get_include()],
    ext_modules = cythonize("brownian_domain_trap_cy.pyx")
)
setup(
	include_dirs=[numpy.get_include()],
    ext_modules = cythonize("brownian_trap_cy.pyx")
)