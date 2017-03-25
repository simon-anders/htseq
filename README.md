[![Build Status](https://travis-ci.org/simon-anders/htseq.svg?branch=master)](https://travis-ci.org/simon-anders/htseq)

# HTSeq
HTSeq is a Python library to facilitate processing and analysis of data from high-throughput sequencing (HTS) experiments. 

### Requirements
To use `HTSeq` you will need:
- `Python 2.7`or `Python >= 3.4` (tested up to 3.6)
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
To install directly from PyPI:
```bash
pip install HTSeq
```
If this fails, please install all dependencies first:
```bash
pip install 'matplotlib>=1.4'
pip install Cython
pip install 'pysam>=0.9'
pip install HTSeq
```
**NOTE**: `pysam==0.9.0` has a bug so that `pip Cython` is __required__ at installation. `pysam>=0.10.0` should build without Cython.

#### Using setup.py (distutils/setuptools)
Install the dependencies with your favourite tool (`pip`, `conda`, etc.).

To install `HTSeq` itself, run:
```bash
python setup.py build install
```

### Documentation
Please see:

   http://www-huber.embl.de/users/anders/HTSeq/

