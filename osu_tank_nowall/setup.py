from setuptools import setup
from Cython.Build import cythonize
from Cython.Distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

setup(
    cmdclass = {'build_ext':build_ext},
    ext_modules = [Extension("mangrove",
                             sources=['mangrove_allbay_hd.pyx'],
                             include_dirs = [numpy.get_include()])]
)
