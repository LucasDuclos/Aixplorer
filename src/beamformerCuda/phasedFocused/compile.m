clear all
delete focusedPhasedBeamformer.o
delete beamformerPhasedFocused.mexglx
! mkdir /tmp/nvcc_compiler/
! ln -sf `which gcc-4.4` /tmp/nvcc_compiler/gcc
! ln -sf `which g++-4.4` /tmp/nvcc_compiler/g++
! nvcc -c --output-file=focusedPhasedBeamformer.o --compiler-bindir=/tmp/nvcc_compiler/ beamformer.cu -Xcompiler -fPIC
mex  CC=gcc-4.4 CXX=g++-4.4 CXXFLAGS="-fPIC -pthread" beamformerPhasedFocused.cpp beamformer_kernel.cpp focusedPhasedBeamformer.o -L/usr/local/cuda/lib -lcuda -lcudart -I/usr/local/cuda/include
delete focusedPhasedBeamformer.o