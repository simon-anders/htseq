[![Build Status](https://travis-ci.org/simon-anders/htseq.svg?branch=master)](https://travis-ci.org/simon-anders/htseq)

# HTSeq
HTSeq is a Python library to facilitate processing and analysis of data from high-throughput sequencing (HTS) experiments. 

**NOTE**: htseq is currently undergoing a major update. The `master` (python2.7) and `python3` (python>=3.4) branches on github require `pysam>=0.9`, but the PyPI version refers to the last stable release, which is only compatible with `pysam<=0.8.0`.We will update the PyPI version within the next few weeks (hopefully!).

### Requirements
To use `HTSeq` you will need:
- `Python ==2.7` (see the `python3` branch for Python 3 support)
- `numpy`
- `pysam >= 0.9.0`

To run the `htseq-qa` script, you will also need:
- `matplotlib >=1.4`

To **build** the package from source, you will **also** need:
- `Cython`
- `SWIG >=3.0.8`

The latter packages are not required if you have already built `HTSeq` and are transferring the binaries onto another machine with a compatible environment (architechture, shared libraries). If you are not sure, chances are you need them.

### Installation
#### PIP
Install all runtime and build dependencies. Because setuptools has historic issues with cython and optional dependencies, you need to install them before anything.
```bash
pip install Cython
pip install 'matplotlib>=1.4'
```

To install `HTSeq` itself, you can download the source from github and run:
```bash
pip install .
```
or install directly from PyPI (once the new version becomes available):
```bash
pip install HTSeq
```

#### Using setup.py (distutils/setuptools)
Install the dependencies with your favourite tool (`pip`, `conda`, etc.).

To install `HTSeq` itself, run:
```bash
python setup.py build install
```

### Documentation
Please see:

   http://www-huber.embl.de/users/anders/HTSeq/

