CXXFLAGS=-I${HOME}/local/eigen-eigen-bdd17ee3b1b3/ -I${HOME}/local/src/Boost.NumPy/

fowf.o: fowf.cpp 
	g++ -c fowf.cpp -m32
fowf:
	ghc -lstdc++ --make fowf.hs fowf.o


fowf.so:fowf.cpp fowfw.cpp fowf.hpp
	g++ -I`python -c 'from distutils.sysconfig import *; print get_python_inc()'` -DPIC -bundle -fPIC -o fowf.so fowf.cpp fowfw.cpp ${CXXFLAGS} -lboost_python  -framework Python

utest: fowf.so utest.py
	python utest.py 



