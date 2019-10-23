from distutils.core import setup, Extension
import numpy
import os

name = "interp1"      # interpolation1 module
version = "1.0"

setup(name=name, version=version,
      ext_modules=[Extension(name='_interp1', 
             sources=["interp.i", "src/interp.c"],
             include_dirs=['src', numpy.get_include(), '.'])
    ])
