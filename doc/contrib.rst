.. _contrib:

**********************
Notes for Contributers
**********************

If you intend to contribute to the development of HTSeq, these notes will
help you to get started.

Source code
-----------

The source code is on an Subversion repository, hosted on SourceForge.

To check out the repository, use
  
::
  
  svn co https://htseq.svn.sourceforge.net/svnroot/htseq/trunk htseq 
   
To browse the repository, see here_.
  
.. _here: http://htseq.svn.sourceforge.net/viewvc/htseq/

Languages
---------

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

I do not want to burden the user with having to install 
SWIG or Cython. Both these tools work by generating C/C++ code which then can
be compiled without the need of any files from SWIG or Cython. Hence, I've
divided the build process into two steps:

* Step 1: Generate the C/C++ files from the SWIG and Cython source files.
  This is done by the calling ``make`` in the ``src`` directory. Note that
  the ``Makefile`` there contains only calls to ``cython`` and ``swig`` but
  not to the C compiler. (Note: I am using Cython 0.11. Compiling with
  Cython 0.12 does not work at the moment, but I will update at some point.)
  
* Step 2: The C files are compiled and copied together with the Python source
  files into a ``build`` directory. This is done by calling ``python setup.py build``
  in the root directory. It creates (as usual for a setup.py script) a new
  directory ``build`` and in it a subdirectory for the machine architecture,
  which then contains the package directory. 
  
To test during development, set the ``PYTHONPATH`` to point to the maschine-specific
directory in the ``build`` directory, so that Python can find the ``HTSeq`` directory
that ``setup.py build`` puts there. Whenever you make a change, call the shell
script ``build_it``, which contains just two lines: the first calls ``make`` in ``src``,
the second calls ``setup.py build``.

Distributing
------------

To wrap up a package, call ``build_it`` (or at least ``make`` in ``src``) 
and then ``setup.py sdist``. This makes a directory ``dists`` and in there,
a tarball with all the source files (Python and C/C++) and all kinds of other stuff
(unfortunately including the ``example_files`` directory, that I hence always delete manually
before running ``setup.py sdist`` to keep the package lean). The tarball contains, when unpacked
the ``setup.py`` script, which allows installing with ``setup.py install``.

I am using setuptools_ (and not distutils_) so that I can make Python eggs with
``setup.py bdist_egg``. For Windows binaries, I use ``setup.py bdist_wininst --compiler=mingw32``
on my virtual Windows box.

.. _setuptools: http://peak.telecommunity.com/DevCenter/setuptools
.. _distutils: http://docs.python.org/library/distutils.html

Files
-----

The package contains the following source files:

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


Furthermore, there are these files to support development:

``src/Makefile``:
  Generates C/C++ files from SWIG and Cython source files but does no C/C++ compiling.
  
``setup.py``:
  A typical setuptools setup.py file.
  
``build_it``:
  A three-line shell script that
  * generates a file ``HTSeq/_version.py`` from the file ``VERSION``.
  * calls ``make`` in ``src`` to process ``src/Makefile``
  * runs ``setup.py build`` (see above)
  
``clean``:
  Another two-line shell script to first call ``make clean`` in ``src`` and then
  delete whatever ``setup.py`` may have written.
  
``test.py``
  Performs all the deoctests in the documentation, using the example data in the
  ``example_data`` directory.
  
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
  
``README``:
  Points the user to the web site.      
  
and these directories

``doc/``:
   this documentation, in Sphinx reStructuredText format, and a Makefile to drive
   Sphinx. 
   
``example_files/``:   
   a few example files to be use by the doctests in the documentation.
