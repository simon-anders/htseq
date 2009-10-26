CYTHON=cython
SWIG=swig

NUMPYDIR=`python -c 'import os.path; import numpy; print os.path.dirname( numpy.__file__ )'`
CFLAGS=`python-config --cflags` -w -shared -fPIC -I $(NUMPYDIR)/core/include

all: _HTSeq.so _StepVector.so StepVector.py

clean:
	rm -f _HTSeq.so _HTSeq.c *.pyc *.gch *_wrap.cxx StepVector.py _*.so \
	   *.pyc *.pyo
   
_HTSeq.so: _HTSeq.c
	$(CC) $(CFLAGS) $^ -o $@

_HTSeq.c: _HTSeq.pyx _HTSeq.pxd
	$(CYTHON) _HTSeq.pyx
   
StepVector_wrap.cxx StepVector.py: StepVector.i AutoPyObjPtr.i
	$(SWIG) -Wall -c++ -python StepVector.i

_StepVector.so: StepVector_wrap.cxx step_vector.h
	$(CXX) $(CFLAGS) StepVector_wrap.cxx -o _StepVector.so
