CYTHON=cython
SWIG=swig

all: StepVector_wrap.cxx ../HTSeq/StepVector.py _HTSeq.c DynamicStepVector_wrap.cxx \
		../HTSeq/DynamicStepVector.py

clean:
	rm -f StepVector_wrap.cxx ../HTSeq/StepVector.py _HTSeq.c DynamicStepVector_wrap.cxx ../HTSeq/DynamicStepVector.py
   
_HTSeq.c: _HTSeq.pyx _HTSeq.pxd
	$(CYTHON) _HTSeq.pyx
	
StepVector_wrap.cxx ../HTSeq/StepVector.py: StepVector.i AutoPyObjPtr.i step_vector.h
	$(SWIG) -Wall -c++ -python StepVector.i 
	mv StepVector.py ../HTSeq

DynamicStepVector_wrap.cxx ../HTSeq/DynamicStepVector.py: DynamicStepVector.i AutoPyObjPtr.i DynamicStepVector.hpp
	$(SWIG) -Wall -c++ -python DynamicStepVector.i 
	mv DynamicStepVector.py ../HTSeq
