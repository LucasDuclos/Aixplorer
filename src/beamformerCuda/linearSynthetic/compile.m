clear all
delete syntheticlinear.o
delete beamformerSyntheticLinear.mexglx
! mkdir /tmp/nvcc_compiler/
! ln -sf `which gcc` /tmp/nvcc_compiler/gcc
! ln -sf `which g++` /tmp/nvcc_compiler/g++
! nvcc -c --output-file=syntheticlinear.o --compiler-bindir=/tmp/nvcc_compiler/ beamformer.cu -Xcompiler -fPIC
mex  CC=gcc CXX=g++ CXXFLAGS="-fPIC -pthread" beamformerSyntheticLinear.cpp beamformer_kernel.cpp syntheticlinear.o -L/usr/local/cuda/lib -lcuda -lcudart -I/usr/local/cuda/include
delete syntheticlinear.o