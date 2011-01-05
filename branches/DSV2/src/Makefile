CYTHON=cython
SWIG=swig

all: _HTSeq.c DynamicStepVector_wrap.cxx ../HTSeq/DynamicStepVector.py

clean:
	rm -f _HTSeq.c DynamicStepVector_wrap.cxx ../HTSeq/DynamicStepVector.py
   
_HTSeq.c: _HTSeq.pyx _HTSeq.pxd
	$(CYTHON) _HTSeq.pyx
   
DynamicStepVector_wrap.cxx ../HTSeq/DynamicStepVector.py: DynamicStepVector.i AutoPyObjPtr.i DynamicStepVector.hpp
	$(SWIG) -Wall -c++ -python DynamicStepVector.i 
	mv DynamicStepVector.py ../HTSeq
