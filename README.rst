|Build Status| |Documentation Status|

**2020-04-11: The new official repository is:** https://github.com/htseq/htseq.

HTSeq
=====

HTSeq is a Python library to facilitate processing and analysis of data
from high-throughput sequencing (HTS) experiments. A popular use of ``HTSeq``
is ``htseq-count``, a tool to quantify gene expression in RNA-Seq and similar
experiments.

Requirements
~~~~~~~~~~~~

To use ``HTSeq`` you will need:

-  ``Python 2.7``\ or ``Python >= 3.5`` (tested up to 3.8)
-  ``numpy``
-  ``pysam >= 0.9.0``

To run the ``htseq-qa`` script, you will also need:

-  ``matplotlib >=1.4``

Both **Linux** and **OSX** are supported and binaries are provided on for many
Linux and OSX versions. A source package which should not require ``Cython``
nor ``SWIG`` is provided for all other cases. To **build** the package completely
from source, you will **also** need:

-  ``Cython >=0.29.5``
-  ``SWIG >=3.0.8``

which are required for performance reasons.

**Windows is not officially supported** as we don't have access to a Continuous
Integration Windows machine that supports ``pysam``. Please do **not** open an
issue asking to support Windows installers: we do not know how to do that and 
do not have the bandwidth to learn. However, if you are interested in giving it
a try yourself, we are happy to provide as much support as we can.

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

To install a specific version (e.g. version 0.11.0):

.. raw:: html

   <div class="highlight highlight-source-shell">

::

    pip install 'HTSeq==0.11.0'

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

.. |Build Status| image:: https://travis-ci.org/htseq/htseq.svg?branch=master
    :target: https://travis-ci.org/htseq/htseq
.. |Documentation Status| image:: https://camo.githubusercontent.com/d3d354c898588bb4b62f559a3a30fa6b6364dfc3/68747470733a2f2f72656164746865646f63732e6f72672f70726f6a656374732f68747365712f62616467652f3f76657273696f6e3d6d6173746572
   :target: http://htseq.readthedocs.io

Authors
~~~~~~~~~~~~~

Since 2016: Fabio Zanini @ http://fabilab.org.
2020-2015: Simon Anders, Wolfgang Huber
