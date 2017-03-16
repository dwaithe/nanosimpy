from setuptools import setup
import numpy
from distutils.extension import Extension
from Cython.Build import cythonize

setup(name='nanosimpy',
      version='0.0.1',
      description='The nano simulation python libary',
      url='https://github.com/dwaithe/nanosimpy',
      author='Dominic Waithe',
      author_email='dominic_waithe@hotmail.com',
      license='MIT',
      packages=['nanosimpy'],
      
      
      include_dirs=[numpy.get_include()],
      ext_modules =  cythonize("nanosimpy/*.pyx"),
      zip_safe=False)