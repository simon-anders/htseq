.. _install:

Prequisites and installation
============================

To use HTSeq, you need at least version 2.5 of Python_ (Python 3 does not work yet), 
together with NumPy_,
a commonly used Python package for numerical calculations. Mac and Linux users 
will find that this is often pre-installed. (To check whether you have a working
NumPy installation, start Python and simply type ``import numpy``. If you do not
get an error, NumPy is available.) Especially Windows users might find the
`Enthought Python Distribution`_ convenient, which contains
Python, NyumPy and other add-ons and is available for all common platforms.

.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _`Enthought Python Distribution`: http://www.enthought.com/products/epd.php

To install HTSeq, get a package tarball from [no web page yet, ask me by e-mail 
for one] and unpack it.

If you have a binary package (not containing a subdirectory called ``src``), open
a terminal and type, in the unpacked directory,
::

   python setup.py install
   
This will not work if you do not have write permission for Python's ``site-package``
directory. In that case, just find the directory ``HTSeq`` inside the unpacked
tarball and add the full path of the directory that contains the ``HTSeq`` directory
to Python's search path. This is done by either adding it to the ``PYTHONPATH``
environment variable or by putting, at the beginning of all your scripts or interactive
sessions::

   import sys
   sys.path.append ('/path/to/the/directory')

If you have a source package (conatining a subdirectory called ``src``), you need to 
first build the package before installing it. For this, make sure you have the GNU tool chain installed. 
(On a Mac, you need to install XCode (available from the Apple web site). On Ubuntu 
Linux, install the ``build-essential`` package. ) Then, open a terminal, ``cd`` 
into the unpacked directory, and type
::

   python setup.py build
  
Afterwards, proceed as with the binary package.

To test your installation, simply start Python and type ``import HTSeq``. No error 
message should appear.


