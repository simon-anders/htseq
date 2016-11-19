CYTHON=cython
SWIG=${SWIG=swig}


all: StepVector_wrap.cxx ../HTSeq/StepVector.py _HTSeq.c

clean:
	rm -f StepVector_wrap.cxx ../HTSeq/StepVector.py _HTSeq.c
   
_HTSeq.c: HTSeq/_HTSeq.pyx HTSeq/_HTSeq.pxd
	$(CYTHON) -3 HTSeq/_HTSeq.pyx -o _HTSeq.c
   
StepVector_wrap.cxx ../HTSeq/StepVector.py: StepVector.i AutoPyObjPtr.i step_vector.h
	$(SWIG) -Wall -c++ -python -py3 StepVector.i 
	
	# Swig doesn't work perfectly with python3 yet
	2to3 -w StepVector.py
	# 2to3 takes care of __next__ and sys.maxsize, but it
	# messes up conditional formatting of import builtins,
	# so we fix it by hand
	sed -i 's/    import builtins as __builtin__/    import builtins/' StepVector.py

	#2to3 does not always know about instances' "next" method
	sed -i 's/\.next/.__next__/' StepVector.py
	
	mv StepVector.py ../HTSeq

