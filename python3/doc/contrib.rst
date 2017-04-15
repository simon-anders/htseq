.. _contrib:

**********************
Notes for Contributors
**********************

If you intend to contribute to the development of HTSeq, these notes will
help you to get started.

Source code
-----------

The source code is on Github_. To check out the repository, use
  
::
  
  git clone https://github.com/simon-anders/htseq.git 
   
.. _Github: https://github.com/simon-anders/htseq

Languages
---------

HTSeq is mostly written in Python and is compatible with both Python 2.7 and
Python 3.4 and above. However, the codebases for Python 2/3 are separate and
development happens mainly on the Python 3 branch.

A good part of HTSeq is actually not written in Python but in 
Cython_. In case you don't know it yet: Cython, a fork from Pyrex, is a
kind of Python compiler. You annotate Python code with additional type
informations (the lines starting with ``cdef`` in the source code). Cython
will then transform the Cython source file (with extension ``pyx``) into
a C file, which calls the appropriate funnctions of Python's C API. Without
type annotation, this looks and feels the same as normal Python and is not 
really faster, either. With type annotation, significant performance gains 
are possible, especially in inner loops.

A small part, namely the StepVector class, is written in C++ and exported with
SWIG. (SWIG_, the "Simple Wrapper and Interface Generator" is a very useful
tool to generate C/C++ code to wrap an existing C/C++ library such that it
becomes accessible as a native library within a number of scripting languages.)
I am not so happy with this any more (the abstraction panelty of the object-oriented 
SWIG wrapping turned out to be a bit high) and ultimatively want to rewrite this
part.

.. _Cython: http://www.cython.org/
.. _SWIG: http://www.swig.org/


Build process
-------------

HTSeq follows the standard python packaging guidelines and relies on a
``setup.py`` script that is simultaneously compatible withPython 2 and 3. To
build the code, run::

  python setup.py build

and to install::

  python setup.py install

If you are not modifying the low-level C/C++/Cython interfaces, you can do
without Cython and SWIG. This is how users normally install HTSeq using
``pip``. If you do modify those files, the ``setup.py`` has a preprocessing
step that calls Cython and/or SWIG if these programs are found. You set
the ``SWIG`` and ``CYTHON`` environment variables to point to your executables
if you have special requirements.
    
To test during development, HTSeq relies on Continuous Integration (CI), at
the moment Travis CI is set up.

To build the documentation, Sphinx_ was used. Just go into the appropriate
``doc`` folder and call::

  make html

to regenerate the documentation. Docs are stored on readthedocs_.

.. _Sphinx: http://www.sphinx-doc.org/
.. _readthedocs: https://readthedocs.org/

Distributing
------------

To wrap up a package, call::

  python setup.py sdist
 
This makes a directory ``dists`` and in there, a tarball with all the source
files (Python and C/C++). If you are a maintainer of HTSeq, you can upload
this file onto PyPI on the testing server. Then, you should run the Tracis CI
tests that try to install HTSeq directly from PyPI (without the source code).
If all goes well, you can upload the tar file onto the live PyPI server.

Files
-----

The package contains source files for Python 2 and 3 in separate folders.
Within each of those folders, the following files are found:

``HTSeq/__init__.py``:
   The outer face of HTSeq. This file defines the name space of HTSeq and contains
   the definition of all classes without performance-critical methods. The file
   imports ``_HTSeq`` in its own namespace, so that, for the user, it does not matter
   whether an object is defined here or in ``_HTSeq.pyx``.
   
``src/HTSeq/_HTSeq.pyx``:
   The core of HTSeq. All classes with perfomance-critical methods are defined here.
   For most of it, this file looks as a normal Python file. Only where performance
   is critical, type annotation has been added. See the Cython manual for details.
   
``src/HTSeq/_HTSeq.pxd``:
   The "header file" for ``_HTSeq.pyx``. It contains the type annotation for all the
   fields of the classes defined in ``_HTSeq.pyx``. If a user would want to write her
   own Cython code, she could use Cython's ``cimport`` directive to import this header
   file and so make Cython aware of the typed definitions of fields and methods in
   ``_HTSeq.pyx``, which may improve performance because it allows Cython to kick out
   all unnecessary type checking.
   
``HTSeq/_HTSeq_internal.py``:
   There are a few limitation to the standard Python code allowed in Cython files;
   most importantly, the ``yield`` statement is not yet supported. Hence, ``_HTSeq.pyx``
   imports this file, and whenever a method in ``_HTSeq.pyx`` needs a ``yield``, 
   it calls a function which is put in here.
   
``src/step_vector.h``:
   The C++ ``step_vector`` class template. As this is a pure template, there is no 
   ``step_vector.cc`` file with definitions. If you want to use a ``step_vector`` in
   a C++ project, this is all you need.
   
``src/StepVector.i``:
   An input file to SWIG, which produces the Python wrapper around ``step_vector.h``, i.e.,
   the ``StepVector`` module containing the ``StepVector`` class. Note that this file contains
   not only SWIG directives but also Python and come C++ code. 
   
``src/AutoPyObjPtr.i``: 
   A very small SWIG library that allows SWIG-wrapped C++ container classes to
   store Python objects in a way that Python's garbage collector is happy with.

``HTSeq/scripts/count.py`` and ``HTSeq/scripts/qa.py``:
   The source code for the stand-alone scripts ``htseq-count`` and ``htseq-qa``.
   They reside in the sub-package ``HTSeq.scripts``, allowing to call the scripts
   with, e.g., ``python -m HTSeq.scripts.qa``.

``scripts/htseq-count`` and ``scripts/htseq-qa``:
   Short stubs to call the scripts from the command line simply as, e.g., ``htseq-qa``.

``doc/``:
   this documentation, in Sphinx reStructuredText format, and a Makefile to drive
   Sphinx. 

``test/test.py``
  Performs all the deoctests in the documentation, using the example data in the
  ``example_data`` directory.

Furthermore, there are these files to support development:

``setup.py``:
  A typical setuptools setup.py file.
  
Finally, there are these files

``VERSION``:
  a one-line text-fil with the version number. It is read by ``setup.py``, used
  by ``build_it`` to generate the one-line Python file ``HTSeq/_version.py`` and
  also used when building the documentation.
  
``MANIFEST.in``:
  Brings some files to the attention of ``setup.py sdist`` which would otherwise not
  be included
  
``LICENCE``:
  The GPL, v3
  
``README.md``:
  Points the user to the web site.      
  
and these directories

``example_files/``:   
   a few example files to be use by the doctests in the documentation.
