#!/usr/bin/env python

from distutils.core import setup, Extension

import os.path
try:
   import numpy
except ImportError:
   sys.stderr.write( "Setup script for HTSeq: Failed to import 'numpy'.\n" )
   sys.stderr.write( "Please install numpy and then try again to install HTSeq.\n" )
   
numpy_include_dir = os.path.join( os.path.dirname( numpy.__file__ ), 'core', 'include' )
 

setup( name = 'HTSeq',
       version = file("VERSION").readline().rstrip(),
       description = 'Processing high-throughput sequencing data',
       author = 'Simon Anders',
       author_email = 'sanders@fs.tum.de',
       url = 'http://htseq.sf.net',
       
       py_modules = [ 
          'HTSeq._HTSeq_internal', 
          'HTSeq.StepVector' 
       ],
       ext_modules = [ 
          Extension( 'HTSeq._HTSeq', 
             ['src/_HTSeq.c'], include_dirs=[numpy_include_dir], extra_compile_args=['-w'] ),
          Extension( 'HTSeq._StepVector', 
             ['src/StepVector_wrap.cxx'], extra_objects=['src/step_vector.h'], extra_compile_args=['-w'] ),
       ]
     )


