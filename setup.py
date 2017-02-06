#!/usr/bin/env python

import sys
import os.path
from Cython.Build import cythonize

python_version_min = '3.4'

try:
    from setuptools import setup, Extension
    from setuptools.command.build_py import build_py
except ImportError:
    sys.stderr.write("Could not import 'setuptools', falling back to 'distutils'.\n")
    from distutils.core import setup, Extension
    from distutils.command.build_py import build_py

if sys.version_info < tuple(map(int, python_version_min.split('.'))):
    sys.stderr.write("Error in setup script for HTSeq:\n")
    sys.stderr.write("Sorry, this version of HTSeq is for Python "+str(python_version_min)+" or higher.\n")
    sys.exit(1)

try:
    import numpy
except ImportError:
    sys.stderr.write("Setup script for HTSeq: Failed to import 'numpy'.\n")
    sys.stderr.write("Please install numpy and then try again to install HTSeq.\n")
    sys.exit(1)

numpy_include_dir = os.path.join(os.path.dirname(numpy.__file__), 'core', 'include')


# Update version from VERSION file into module
with open('VERSION') as fversion:
    version = fversion.readline().rstrip()
with open('HTSeq/_version.py', 'wt') as fversion:
    fversion.write('__version__ = "'+version+'"')


class Build_with_preprocess(build_py):
    def swig_and_cython(self):
        import os
        from subprocess import call

        def c(x): return call(x, shell=True)

        # CYTHON
        print('Cythonizing...')
        cython = os.getenv('CYTHON', 'cython')
        c(cython+' -3 src/HTSeq/_HTSeq.pyx -o src/_HTSeq.c')
        print('done')

        # SWIG
        print('SWIGging...')
        swig = os.getenv('SWIG', 'swig')
        c(swig+' -Wall -c++ -python -py3 src/StepVector.i')
        print('done, correcting and moving...')
        pyswigged = 'src/StepVector.py'
        c("2to3 --no-diffs --write --nobackups "+pyswigged)
        c("sed -i 's/    import builtins as __builtin__/    import builtins/' "+pyswigged)
        c("sed -i 's/\.next/.__next__/' "+pyswigged)
        c('mv '+pyswigged+' HTSeq/StepVector.py')
        print('done')

    def run(self):
        self.swig_and_cython()
        super().run()


setup(name='HTSeq',
      version=version,
      author='Simon Anders',
      author_email='sanders@fs.tum.de',
      maintainer='Fabio Zanini',
      maintainer_email='fabio.zanini@stanford.edu',
      url='http://www-huber.embl.de/users/anders/HTSeq/',
      description="A framework to process and analyze data from " +
                  "high-throughput sequencing (HTS) assays",
      classifiers=[
         'Development Status :: 5 - Production/Stable',
         'Topic :: Scientific/Engineering :: Bio-Informatics',
         'Intended Audience :: Developers',
         'Intended Audience :: Science/Research',
         'License :: OSI Approved :: GNU General Public License (GPL)',
         'Operating System :: POSIX',
         'Programming Language :: Python'
      ],
      requires=['numpy', 'python (>='+python_version_min+')'],
      ext_modules=[
         Extension(
             'HTSeq._HTSeq',
             ['src/_HTSeq.c'],
             include_dirs=[numpy_include_dir],
             extra_compile_args=['-w']),
         Extension(
             'HTSeq._StepVector',
             ['src/StepVector_wrap.cxx'],
             extra_compile_args=['-w']),
      ],
      py_modules=[
         'HTSeq._HTSeq_internal',
         'HTSeq.StepVector',
         'HTSeq._version',
         'HTSeq.scripts.qa',
         'HTSeq.scripts.count'
      ],
      scripts=[
         'scripts/htseq-qa',
         'scripts/htseq-count',
      ],
      cmdclass={
          'build_py': Build_with_preprocess,
          }
      )
