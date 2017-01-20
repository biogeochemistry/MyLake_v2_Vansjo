mex -v CC=gcc LD=g++ COPTIMFLAGS='-O3 -DNDEBUG -ffast-math'  -I/Users/MarkelovIgor/eigen pH.cpp

mex -v CC=gcc LD=g++ COPTIMFLAGS='-O3 -DNDEBUG -ffast-math' -I/usr/local/include/ -I/Users/MarkelovIgor/eigen -L/usr/local/lib/ -liphreeqc pH_phreeqc.cpp