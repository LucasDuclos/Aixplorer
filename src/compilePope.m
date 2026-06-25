try
    currentPath = pwd;
    cd ../lib/popeRubi;

    if isunix
        delete('*.a');
        delete('*.o');

        %ssibuntu version
%         mex  CC=gcc CXX=g++ -c CXXFLAGS="-g -fPIC -pthread -msse -msse2 -msse3" ../../src/popeRubi/reconstruct.cpp
%         mex  CC=gcc CXX=g++ -c CXXFLAGS="-g -fPIC -pthread -msse -msse2 -msse3" ../../src/popeRubi/complut.cpp

        %other distrib version
        mex  -c CXXFLAGS="-ggdb -fPIC -pthread -msse -msse2" ../../src/popeRubi/reconstruct.cpp
        mex  -c CXXFLAGS="-ggdb -fPIC -pthread -msse -msse2" ../../src/popeRubi/complut.cpp

        mex  CC=gcc CXX=g++ -c CXXFLAGS="-g -fPIC -pthread -msse -msse2 -msse3" ../../src/popeRubi/CompLutLinear.cpp 
        mex  CC=gcc CXX=g++ -c CXXFLAGS="-g -fPIC -pthread -msse -msse2 -msse3" ../../src/popeRubi/CompLutPhased.cpp 
        mex  CC=gcc CXX=g++ -c CXXFLAGS="-g -fPIC -pthread -msse -msse2 -msse3" ../../src/popeRubi/CompLutLinearPhotoacoustics.cpp
        mex  CC=gcc CXX=g++ -c CXXFLAGS="-g -fPIC -pthread -msse -msse2 -msse3" ../../src/popeRubi/CompLutCurvedSteering.cpp 
        mex  CC=gcc CXX=g++ -c CXXFLAGS="-g -fPIC -pthread -msse -msse2 -msse3" ../../src/popeRubi/CompLut_FlatSynthetic_Curved.cpp

        mex  CC=gcc CXX=g++ ../../src/popeRubi/reconstructRubi.cpp reconstruct.o
        mex  CC=gcc CXX=g++ ../../src/popeRubi/complutRubi.cpp complut.o CompLutLinear.o CompLutPhased.o CompLutLinearPhotoacoustics.o CompLutCurvedSteering.o CompLut_FlatSynthetic_Curved.o

        delete('*.a');
        delete('*.o');

%         mex  ../popeRubi/reconstructRubi.cpp reconstruct.o
%         mex  ../popeRubi/complutRubi.cpp complut.o

    elseif ispc
        delete('*.obj');
        delete('*.mexw64*');
        delete('*.mexw32*');
        mex  -c ../../src/popeRubi/reconstruct.cpp
        mex  -c ../../src/popeRubi/complut.cpp
        mex  -c ../../src/popeRubi/CompLutLinear.cpp
        mex  -c ../../src/popeRubi/CompLutPhased.cpp
        mex  -c ../../src/popeRubi/CompLutLinearPhotoacoustics.cpp
        mex  -c ../../src/popeRubi/CompLutCurvedSteering.cpp
        mex  -c ../../src/popeRubi/CompLut_FlatSynthetic_Curved.cpp

        mex  CXXFLAGS="-fpermissive" ../../src/popeRubi/reconstructRubi.cpp reconstruct.obj
        mex  ../../src/popeRubi/complutRubi.cpp complut.obj CompLutLinear.obj CompLutPhased.obj CompLutLinearPhotoacoustics.obj CompLutCurvedSteering.obj CompLut_FlatSynthetic_Curved.obj

        delete('*.obj');
%         delete('*.mexw64*');
%         delete('*.mexw32*');

    else
        error ('Platform Issue. This platform is not known');
    end

    addpath(pwd)

    disp ('script finished with success!')
    cd (currentPath);

catch ME
    disp ('Error during remote compilation!')
    disp (ME.message)
    disp (ME.identifier)
    cd (currentPath);
end
