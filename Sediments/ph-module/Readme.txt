Instead of compiling MEX in command line within Matlab, an alternative way is to use CMake tools. 
Using CMake could be advantageous for building large MEX project or building MEX with lots of external dependencies. 
This example is a simple demonstration of how to easily compile Matlab MEX using CMake. 

To compile the test MEX under Linux, 
first set MATLAB_ROOT environment variable to your installed matlab path,
such as 'export MATLAB_ROOT=/usr/local/MATLAB/R2012b',
then, simply do

mkdir build
cd build
cmake ../src
make
make install

To compile the test MEX under Windows,
first set MATLAB_ROOT environment variable to your installed matlab path,
then, use cmake or cmake-gui to generate building project according to installed compiler (e.g. MSVC),
then, build the generated project using this compiler.

The test MEX source code is located under /src/mex/mexAdd. The compiled test MEX 'mexAdd' will 
be installed into /bin by default. C=mexAdd(A,B) basically do element-by-element addition for 1D or 
2D matrix A and B, return matrix C. Nothing fancy.

To add new MEX source code, for example mexXXX.cpp, simply do
1. add a new folder 'mexXXX' under /src/mex
2. add one line 'add_subdirectory(mexXXX)' to CMakeLists.txt under /src/mex
3. copy CMakeLists.txt under /src/mex/mexAdd to folder /src/mex/mexXXX
4. change first line to set(CPP_FILE mexXXX) in copied CMakeLists.txt
5. follow compiling steps as described above

Another two examples are also available for
building MEX with CUDA support using CMake
http://www.mathworks.com/matlabcentral/fileexchange/45505-cudamexcmake
building MEX with OpenMP support using CMake
http://www.mathworks.com/matlabcentral/fileexchange/45501-openmpmexcmake