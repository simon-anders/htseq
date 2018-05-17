|Build Status| |Documentation Status|

HTSeq
=====

HTSeq is a Python library to facilitate processing and analysis of data
from high-throughput sequencing (HTS) experiments.

Requirements
~~~~~~~~~~~~

To use ``HTSeq`` you will need:

-  ``Python 2.7``\ or ``Python >= 3.4`` (tested up to 3.6)
-  ``numpy``
-  ``pysam >= 0.9.0``

To run the ``htseq-qa`` script, you will also need:

-  ``matplotlib >=1.4``

To **build** the package from source, you will **also** need:

-  ``Cython``
-  ``SWIG >=3.0.8``

The latter packages are not required if you have already built ``HTSeq``
and are transferring the binaries onto another machine with a compatible
environment (architechture, shared libraries). If you are not sure,
chances are you need them.

Both **Linux** and **OSX** are supported and binaries are provided for virtually
all Linux versions and for some OSX versions (the latter only for Python 2.7
and 3.6). A source package which should not require ``Cython`` nor ``SWIG``
is provided for all other cases.

**Windows is not officially supported** as we don't have access to a Continuous
Integration Windows machine that supports ``pysam``. However, if you have built
``HTSeq`` for Windows, please open an issue and we'll try and include it in the
release.

Installation
~~~~~~~~~~~~

PIP
^^^

To install directly from PyPI:

.. raw:: html

   <div class="highlight highlight-source-shell">

::

    pip install HTSeq

.. raw:: html

   </div>

If this fails, please install all dependencies first:

.. raw:: html

   <div class="highlight highlight-source-shell">

::

    pip install 'matplotlib>=1.4'
    pip install Cython
    pip install 'pysam>=0.9'
    pip install HTSeq

.. raw:: html

   </div>

**NOTE**: ``pysam==0.9.0`` has a bug so that ``pip Cython`` is
**required** at installation. ``pysam>=0.10.0`` should build without
Cython.

Using setup.py (distutils/setuptools)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the dependencies with your favourite tool (``pip``, ``conda``,
etc.).

To install ``HTSeq`` itself, run:

.. raw:: html

   <div class="highlight highlight-source-shell">

::

    python setup.py build install

.. raw:: html

   </div>

Documentation
~~~~~~~~~~~~~

Please see:

http://htseq.readthedocs.io

.. |Build Status| image:: https://camo.githubusercontent.com/12452733a10aadd3dfd477d0497f2f4a32935be3/68747470733a2f2f7472617669732d63692e6f72672f73696d6f6e2d616e646572732f68747365712e7376673f6272616e63683d6d6173746572
   :target: https://travis-ci.org/simon-anders/htseq
.. |Documentation Status| image:: https://camo.githubusercontent.com/d3d354c898588bb4b62f559a3a30fa6b6364dfc3/68747470733a2f2f72656164746865646f63732e6f72672f70726f6a656374732f68747365712f62616467652f3f76657273696f6e3d6d6173746572
   :target: http://htseq.readthedocs.io
