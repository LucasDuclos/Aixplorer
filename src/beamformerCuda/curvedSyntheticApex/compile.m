! nvcc -c --output-file=temp.o beamformer.cu -Xcompiler -fPIC
%! nvcc -c --output-file=syntheticcurved.o  --compiler-bindir="C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin" -I="C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\include" beamformer.cu
mex  beamformerSyntheticCurvedApex.cpp temp.o -lcuda -lcudart
delete temp.o