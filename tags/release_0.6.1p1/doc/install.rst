.. _install:

****************************
Prequisites and installation
****************************

HTSeq is available from the `Python Package Index (PyPI)`_. Download HTSeq from the

.. _`Python Package Index (PyPI)`: http://pypi.python.org/

  `HTSeq package download page on PyPI`_,

.. _`HTSeq package download page on PyPI`: http://pypi.python.org/pypi/HTSeq
 
where you will find various version of the package. 

To use HTSeq, you need at least version 2.5 of Python_ (Python 3 does not work yet), 
together with NumPy_,
a commonly used Python package for numerical calculations,
and matplotlib_, a plotting library. 

.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://matplotlib.org/

HTSeq can be installed like any other Python package. For users unfamiliar with this,
more detailled instructions are given below.


Installation on Linux
=====================

Make sure that you the standard GNU build environment installed, as well as Python together with its development headers and numpy and matplotlib. Users of Ubuntu Linux simply type::

   sudo apt-get install build-essential python2.7-dev python-numpy python-matplotlib

This command installs the required packages as listed (or simply does nothing if they are already installed.

For users of RedHat or RedHat-derived distros (Fedora, CentOS), the equivalent command
seems to be (untested)::

   sudo yum groupinstall "Development Tools"
   sudo yum install python-devel numpy python-matplotlib

To install HTSeq itself, download the *source* package from the `HTSeq PyPI page`_, unpack the tarball,
go into the directory with the unpacked files and type there

.. _`HTSeq PyPI page`: http://pypi.python.org/pypi/HTSeq

::

   python setup.py install --user

to install HTSeq for the user currently logged in. To make HTSeq available to all users, use instead::

   python setup.py build
   sudo python setup.py install

To test the installation, change to another director than the build directory, start Python
(by typing ``python``) and then try whether typing ``import HTSeq`` causes an error meesage.

Installation on MacOS X
=======================

Mac users should install NumPy as explained here_ in the NumPy/SciPy documentation. Note that you need
to install Xcode to be able to compile NumPy. Due to the
mess that Apple recently made out of Xcode, the whole process may be a slight bit more cumbersome than necessary, especially if you work with OSX Lion, so read the instructions carefully.

.. _here: http://www.scipy.org/Installing_SciPy/Mac_OS_X

If you want to produce plots or use htseq-qa, you will also need matplotlib. (For htseq-count, it
is not required.) There seems to be a binary package (a "Python egg") available on the matplotlib
SourceForge page.

To install HTSeq itself, download the *source* package from the `HTSeq PyPI page`_, unpack the tarball,
go into the directory with the unpacked files and type there:

.. _`HTSeq PyPI page`: http://pypi.python.org/pypi/HTSeq

::

   python setup.py build

to compile HTSeq. If you get an error regarding the availability of a C compiler, you may need to
set environment variables to point Python to the . The NumPy/SciPy installation instructions above cover this topic well and
apply here, too, so simply do the same as you did to install NumPy.

Once building has been successful, use::

   python setup.py --user

to install HTSeq for the current users. To make HTSeq available to all users, use instead::

   python setup.py build
   sudo python setup.py install

To test the installation, change to another director than the build directory, start Python
(by typing ``python``) and then try whether typing ``import HTSeq`` causes an error meesage.


MS Windows
==========

If you have not yet installed Python, do so first. You can find an automatic installer 
for Windows on
the `Python download page`_. Make sure to use Python 2.7, not Python 3.3. 

.. _`Python download page`: http://www.python.org/getit/

Then install the newest version of NumPy. Look on `NumPy's PyPI page`_ for the automatic installer.

.. _`NumPy's PyPI page`: https://pypi.python.org/pypi/numpy

If you want to produce plots or use htseq-qa, you will also need matplotlib. (For htseq-count, it
is not required.) Follow the installation instructions on their web page.

To install HTSeq itself, simply download the Windows installer from the `HTSeq download page`_
and run it.

.. _`HTSeq download page`: http://pypi.python.org/pypi/HTSeq

To test your installation, start Python and then try whether typing ``import HTSeq`` 
causes an error meesage.

If you get the error message "ImportError: DLL load failed", you are most likely
missing the file MSVCR110.DLL on your system, which you can get by downloading and
installing the file "vcredist_x86.exe" from `this page`_.

.. _`this page`: http://www.microsoft.com/en-us/download/details.aspx?id=30679

