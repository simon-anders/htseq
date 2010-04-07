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
get an error, NumPy is available.) 

.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _`Enthought Python Distribution`: http://www.enthought.com/products/epd.php


Download
========

HTSeq is available from the `Python Package Index (PyPI)`_. 

.. _`Python Package Index (PyPI)`: http://pypi.python.org/

If you have ``easy_install`` (part of the Python `setuptools`_) you can simply
type

.. _`setuptools`: http://pypi.python.org/pypi/setuptools

::

   easy_install HTSeq
   
and ``easy_install`` will download and install everything automatically.   

Alternatively, you can download HTSeq from the

  `HTSeq package download page on PyPI`_,

.. _`HTSeq package download page on PyPI`: http://pypi.python.org/pypi/HTSeq
 
where you will find various version of the package. Either pick a *binary*
archive suitable for your operation system and Python version, or download the
*source* tarball, which should work with all supported operating systems and Python
versions.


Installing a binary package
===========================

If you have a binary package (not containing a subdirectory called ``src``),
unpack it, open a terminal and type, in the unpacked directory,
::

   python setup.py install
   
If you do not have write permission for Python's ``site-package``
directory, you can install it locally in your home directory by typing

::

   python setup.py install --user


To test your installation, simply start Python (*not* in the directory
with ``setup.py``) and type ``import HTSeq``. No error 
message should appear.

If you have problems with the installation, please do not hesitate to contact me
(anders *at* embl *dot* de).


Installing a source package
===========================

If you have a source package (containing a subdirectory called ``src``), you need to 
first build the package before installing it. For this, make sure you have the
GNU tool chain installed:

* On a Mac, you need to install XCode (available from the Apple web site).

* On Ubuntu or Debian Linux, install the ``build-essential`` package. For other
  Linux distributions, similar packages are available.

* On MS Windows, MinGW_ is a commonly used build environment. Using it may be
  a bit tricky, so use a binary package if possible.

.. _MinGW: http://www.mingw.org/

Unpack the source package, open a terminal, ``cd`` 
into the unpacked directory, and type
::

   python setup.py build
  
Afterwards, proceed as with the binary package.


