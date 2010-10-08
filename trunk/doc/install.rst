.. _install:

****************************
Prequisites and installation
****************************

Prequisites
===========

To use HTSeq, you need at least version 2.5 of Python_ (Python 3 does not work yet), 
together with NumPy_,
a commonly used Python package for numerical calculations. Mac and Linux users 
will find that this is often pre-installed. (To check whether you have a working
NumPy installation, start Python and simply type ``import numpy``. If you do not
get an error, NumPy is available. Users of Debian and Ubuntu Linux should install
the package ``python-numpy``.) 

.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _`Enthought Python Distribution`: http://www.enthought.com/products/epd.php


Download
========

HTSeq is available from the `Python Package Index (PyPI)`_. 

.. _`Python Package Index (PyPI)`: http://pypi.python.org/

Download HTSeq from the

  `HTSeq package download page on PyPI`_,

.. _`HTSeq package download page on PyPI`: http://pypi.python.org/pypi/HTSeq
 
where you will find various version of the package. Either pick a *binary*
archive suitable for your operation system and Python version, or download the
*source* tarball, which should work with all supported operating systems and Python
versions.


Installing a binary package
===========================

For Windows, the binary package contains an automatic installer, Simply execute it.

If you have a binary package (not containing a subdirectory called ``src``) for Linux or MacOS,
unpack it, open a terminal and type, in the unpacked directory,
::

   python setup.py install
  
If you do not have write permission for Python's ``site-package``
directory, you can install it locally in your home directory as follows: 

   *Python 2.6* users simply type

   ::

      python setup.py install --user

   to get it installed in a place in your home directory (typically ``~/.local``) 
   where Python will find it.

   For a local installation with *Python 2.5*, you have to specify a 
   place where to install it to:

   ::

      python setup.py install --home <some_path>
      
   and then tell Python (every time before using HTSeq) where to look for it by setting
   the environment variable ``PYTHONPATH``:

   ::

     export PYTHONPATH=$PYTHONPATH:<some_path>/lib/python   


To test your installation, simply start Python (**not** in the directory
with ``setup.py``) and type ``import HTSeq``. No error 
message should appear.

If you have problems with the installation, please do not hesitate to contact me
(anders *at* embl *dot* de).


Installing a source package
===========================

If you have a source package (containing a subdirectory called ``src``), you need to 
first build the package before installing it. For this, make sure you have the
GNU tool chain and the Python header files installed:

* On a Mac, you need to install XCode (available from the Apple web site or from
  your second Mac OS X install CD).

* On Ubuntu or Debian Linux, install the ``build-essential`` and the 
  ``python-dev`` package. For other Linux distributions, similar packages are available.

* On MS Windows, MinGW_ is a commonly used build environment. Using it may be
  a bit tricky, so use a binary package if possible. If non is available, see below
  for detailed instructions.

.. _MinGW: http://www.mingw.org/

Unpack the source package, open a terminal, ``cd`` 
into the unpacked directory, and type
::

   python setup.py build
  
Afterwards, proceed as with the binary package.


Installation of a source package on MS Windows
..............................................

If the current version is not available as binary package for Windows, either 
ask me to build a Windows package, or proceed as follows:

- In addition to Python and numpy, also install MinGW_.

.. _MinGW: http://www.mingw.org/

- Download the source package and unpack it (e.g., with 7-zip) onto your desktop.

- Start a command line terminal and type the following commands (substituting paths
  where they are different on your system):

  ::  
  
     cd \Documents and Settings\anders\Desktop\HTSeq-0.4.4
     PATH=%PATH%;C:\Python26;C:\MinGW\bin
     python setup.py build --compiler=mingw32
     python setup.py install
    
