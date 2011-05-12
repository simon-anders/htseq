CYTHON=cython
SWIG=swig

all: StepVector_wrap.cxx ../HTSeq/StepVector.py _HTSeq.c

clean:
	rm -f StepVector_wrap.cxx ../HTSeq/StepVector.py _HTSeq.c
   
_HTSeq.c: HTSeq/_HTSeq.pyx HTSeq/_HTSeq.pxd
	$(CYTHON) HTSeq/_HTSeq.pyx -o _HTSeq.c
   
StepVector_wrap.cxx ../HTSeq/StepVector.py: StepVector.i AutoPyObjPtr.i step_vector.h
	$(SWIG) -Wall -c++ -python StepVector.i 
	mv StepVector.py ../HTSeq

