%sclear all
delete beamformerCurvedSynthetic.mex*
! mkdir /tmp/nvcc_compiler/
! ln -sf `which gcc-4.4` /tmp/nvcc_compiler/gcc
! ln -sf `which g++-4.4` /tmp/nvcc_compiler/g++
! nvcc -c --output-file=syntheticCurvedBeamformer.o --use_fast_math -O3 --compiler-bindir=/tmp/nvcc_compiler/ beamformer.cu -Xcompiler -fPIC -DCUDA_CODE
! nvcc -c --output-file=syntheticCurvedBeamformerCPU.o --compiler-bindir=/tmp/nvcc_compiler/ beamformer_kernel.cu -Xcompiler -fPIC
mex  -output beamformerCurvedSynthetic CC=gcc-4.4 CXX=g++-4.4 CXXFLAGS="-fPIC -pthread" beamformerGPU.cpp syntheticCurvedBeamformerCPU.o syntheticCurvedBeamformer.o -L/usr/local/cuda/lib -lcuda -lcudart -I/usr/local/cuda/include
delete syntheticCurvedBeamformer*.o