CYTHON=cython
SWIG=swig

all: StepVector_wrap.cxx ../HTSeq/StepVector.py _HTSeq.c

clean:
	rm -f _HTSeq.so _HTSeq.c *.pyc *.gch *_wrap.cxx StepVector.py _*.so \
	   *.pyc *.pyo
   
_HTSeq.c: _HTSeq.pyx _HTSeq.pxd
	$(CYTHON) _HTSeq.pyx
   
StepVector_wrap.cxx ../HTSeq/StepVector.py: StepVector.i AutoPyObjPtr.i
	$(SWIG) -Wall -c++ -python StepVector.i 
	mv StepVector.py ../HTSeq

