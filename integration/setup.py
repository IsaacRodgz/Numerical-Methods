from distutils.core import setup, Extension
import numpy
import os

name = "interp2"      # interpolation1 module
version = "1.0"

setup(name=name, version=version,
      ext_modules=[Extension(name='_interp2',
             sources=["interp.i", "src/interp2.c"],
             include_dirs=['src', numpy.get_include(), '.'])
    ])
