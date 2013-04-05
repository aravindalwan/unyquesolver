import distribute_setup
distribute_setup.use_setuptools()

from setuptools import setup, Extension
import sys, os
import pkg_resources
import numpy

here = os.path.abspath(os.path.dirname(__file__))
README = open(os.path.join(here, 'README.rst')).read()
NEWS = open(os.path.join(here, 'NEWS.txt')).read()

version = '0.1'

install_requires = [
    # List your project dependencies here.
    # For more details, see:
    # http://packages.python.org/distribute/setuptools.html#declaring-dependencies
    ]

numpy_inc = numpy.get_include()
pyublas_inc = pkg_resources.resource_filename('pyublas', 'include')

setup(name='unyquesolver',
      version=version,
      description='Uncertainty Quantification Environment for modeling' + \
          ' uncertainties in engineering systems',
      long_description=README + '\n\n' + NEWS,
      classifiers=[
        # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        ],
      keywords='uncertainty estimation FEM',
      author='Aravind Alwan',
      author_email='aalwan2@illinois.edu',
      url='',
      license='GPLv3',
      packages=['unyquesolver'],
      ext_modules = [
        Extension('unyquesolver._internals', [
                'unyquesolver/internals/util/ublas.cpp',
                'unyquesolver/internals/fem/function.cpp',
                'unyquesolver/internals/fem/fem.cpp',
                'unyquesolver/internals/fem/nonelast.cpp',
                'unyquesolver/internals/fem/therm.cpp',
                'unyquesolver/internals/fem/elec.cpp',
                'unyquesolver/internals/fem/eles.cpp',
                'unyquesolver/internals/fem/fluid.cpp',
                'unyquesolver/internals/fem/solver.cpp',
                'unyquesolver/internals/wrapper.cpp',
                ],
                  include_dirs = [
                numpy_inc,
                pyublas_inc,
                'unyquesolver/internals/util',
                'unyquesolver/internals/fem',
                ],
                  libraries=[
                ],
                  depends = [
                'unyquesolver/internals/util/ublas.hpp',
                'unyquesolver/internals/fem/function.hpp',
                'unyquesolver/internals/fem/fem.hpp',
                'unyquesolver/internals/fem/nonelast.hpp',
                'unyquesolver/internals/fem/therm.hpp',
                'unyquesolver/internals/fem/elec.hpp',
                'unyquesolver/internals/fem/eles.cpp',
                'unyquesolver/internals/fem/fluid.hpp',
                'unyquesolver/internals/fem/solver.hpp',
                ],
                  ),
        ],
      zip_safe=False,
      install_requires=install_requires,
      entry_points={
        'console_scripts':
            []
        }
      )
