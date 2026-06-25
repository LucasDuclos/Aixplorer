
function compileRemote (varargin)
% COMPILEREMOTE  Compile a remote client for matlab
%   Work under both linux64/32(gcc) 
%       and windows XP32 & VISTA or 7 64/32 (VS2008 or Intel)
%
%
% Usage :
%--------
%
% * LINUX   : compileRemote
%
% * WINDOWS : compileRemote([libLocation], [boost compiler], [boost version])
%    with :
%       liblocation :  - 'relative'(default) the libraries to link are in
%                           ..\..\dependancies
%                      - 'programfiles' the libraries to link
%                           are in C:\Program Files\...
%
%       boost compiler : - vs2008(default) normally used with matlab. With
%                           64 bits Windows version use intel to compile
%                           the pope (MMX + VS C++ compiler = caca)
%                         - intel : To use with 64bits windows to get pope.
%                           Please install the compiler before use...
%
%       boost version  : - 40 (default) or 44
%
% DEPENDANCIES. Remote needs :
%       - xerces_c 2.8  header + binaries
%       - boost 1.40    header + binaries
%       - asio 1.4.1    binaries
%
%Copyright 2010 Supersonic Imagine
%$Revision: 1.00$  $Date: 2010/09/28$
%M-File function.    
    try 
        currentPath = pwd;
        cd remoteBuild;

        if isunix
            delete('*.a');
            delete('*.o');
            delete('*.m');
            copyfile('../src/*.m', '.')
        
            disp ('Building...');
            mex  -c -g ../libClientCommonRemote/message/messageCommon.cpp
            mex  -c -g ../libClientCommonRemote/message/CRemoteMessage.cpp -lxerces-c
            mex   -c -g ../libClientCommonRemote/message/xercesCommon.cpp
            mex   -c -g ../libClientCommonRemote/message/CCell.cpp
            mex   -c -g ../libClientCommonRemote/message/CCellContainer.cpp
            mex   -c -g ../libClientCommonRemote/message/CRemoteSequence.cpp -lxerces-c
            system('ar rcs libMessageRemote.a xercesCommon.o messageCommon.o CRemoteMessage.o CCell.o CCellContainer.o CRemoteSequence.o');

            mex   -c -g ../libClientCommonRemote/common.cpp
            mex   -c -g ../libClientCommonRemote/commonError.cpp
            mex   -c -g ../libClientCommonRemote/commonAsio.cpp
            mex   -c -g ../src/matlabParser.cpp

            disp ('Linking...');
            mex   -g ../src/remoteGetSweData.cpp libMessageRemote.a common.o commonError.o matlabParser.o commonAsio.o -lxerces-c
            mex   -g ../src/remoteSendMessage.cpp libMessageRemote.a common.o commonError.o matlabParser.o commonAsio.o -lxerces-c
            mex   -g ../src/remoteSendSequence.cpp libMessageRemote.a common.o commonError.o matlabParser.o commonAsio.o -lxerces-c
            mex   -g ../src/remoteGetRFData.cpp libMessageRemote.a common.o commonError.o matlabParser.o commonAsio.o -lxerces-c
            mex   -g ../src/remoteGetShmRFData.cpp libMessageRemote.a common.o commonError.o matlabParser.o commonAsio.o -lxerces-c
            mex   -g ../src/remoteGetBFData.cpp libMessageRemote.a common.o commonError.o matlabParser.o commonAsio.o -lxerces-c
            mex   -g ../src/remoteGetFFData.cpp libMessageRemote.a common.o commonError.o matlabParser.o commonAsio.o -lxerces-c
            mex   -g ../src/remoteGetHelixData.cpp libMessageRemote.a common.o commonError.o matlabParser.o commonAsio.o -lxerces-c


        elseif ispc %windows mode

            lib.target = 'relative'; %default value
            lib.boostVersion = '40';
            lib.boostCompiler = 'vs2008';    % intel or vs2008. Intel mandatory for Windows X64

            for k = 1:length(varargin)
                switch varargin{k}
                    case 'relative'
                        lib.target = 'relative';
                    case 'intel'
                        lib.boostCompiler = 'intel';
                    case '44'
                        lib.boostVersion = '44';
                end
            end

            switch lib.boostCompiler
                case 'intel'
                    lib.boostPostfix = ['-iw-mt-1_' lib.boostVersion];
                case 'vs2008'
                    lib.boostPostfix = ['-vc90-mt-1_' lib.boostVersion];
            end
            message = ['Launch windows build with dependancies in "' lib.target '" and for a boost build "' lib.boostCompiler '" version ' lib.boostVersion];
            disp (message);

            delete('*.obj');
            delete('*.pdb');
            delete('*.m');
            copyfile('../src/*.m', '.');

            % use mex with this syntax to force the script to resolve the
            % names in 8.3 format. In other case, the spaces (like program files)
            % should divide one argument in many arguments
            
            % -D_WIN32_WINNT=0x0501 => to avoid warnings
            disp ('Building...');
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/message/messageCommon.cpp')
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/message/CRemoteMessage.cpp', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib))
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/message/xercesCommon.cpp', getWinArg('xerces', 'include', lib))
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/message/CCell.cpp')
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/message/CCellContainer.cpp')
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/message/CRemoteSequence.cpp', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib))


            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/common.cpp')
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/commonError.cpp')
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../libClientCommonRemote/commonAsio.cpp', getWinArg('asio', 'include', lib), getWinArg('boost', 'include', lib), getWinArg('boost', 'lib', lib))
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-c', '-g', '../src/matlabParser.cpp')

            disp ('Linking...');
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-g', '../src/remoteSendMessage.cpp', 'messageCommon.obj','CRemoteMessage.obj','xercesCommon.obj','CCell.obj','CCellContainer.obj','CremoteSequence.obj','common.obj','commonError.obj','commonAsio.obj','matlabParser.obj', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib), getWinArg('asio', 'include', lib), getWinArg('boost', 'include', lib), getWinArg('boost', 'lib', lib), ['-llibboost_date_time' lib.boostPostfix],['-llibboost_regex' lib.boostPostfix])
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-g', '../src/remoteSendSequence.cpp', 'messageCommon.obj','CRemoteMessage.obj','xercesCommon.obj','CCell.obj','CCellContainer.obj','CremoteSequence.obj','common.obj','commonError.obj','commonAsio.obj','matlabParser.obj', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib), getWinArg('asio', 'include', lib), getWinArg('boost', 'include', lib), getWinArg('boost', 'lib', lib), ['-llibboost_date_time' lib.boostPostfix],['-llibboost_regex' lib.boostPostfix])
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-g', '../src/remoteGetRFData.cpp', 'messageCommon.obj','CRemoteMessage.obj','xercesCommon.obj','CCell.obj','CCellContainer.obj','CremoteSequence.obj','common.obj','commonError.obj','commonAsio.obj','matlabParser.obj', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib), getWinArg('asio', 'include', lib), getWinArg('boost', 'include', lib), getWinArg('boost', 'lib', lib), ['-llibboost_date_time' lib.boostPostfix],['-llibboost_regex' lib.boostPostfix])
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-g', '../src/remoteGetBFData.cpp', 'messageCommon.obj','CRemoteMessage.obj','xercesCommon.obj','CCell.obj','CCellContainer.obj','CremoteSequence.obj','common.obj','commonError.obj','commonAsio.obj','matlabParser.obj', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib), getWinArg('asio', 'include', lib), getWinArg('boost', 'include', lib), getWinArg('boost', 'lib', lib), ['-llibboost_date_time' lib.boostPostfix],['-llibboost_regex' lib.boostPostfix])
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-g', '../src/remoteGetFFData.cpp', 'messageCommon.obj','CRemoteMessage.obj','xercesCommon.obj','CCell.obj','CCellContainer.obj','CremoteSequence.obj','common.obj','commonError.obj','commonAsio.obj','matlabParser.obj', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib), getWinArg('asio', 'include', lib), getWinArg('boost', 'include', lib), getWinArg('boost', 'lib', lib), ['-llibboost_date_time' lib.boostPostfix],['-llibboost_regex' lib.boostPostfix])
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-g', '../src/remoteGetSweData.cpp', 'messageCommon.obj','CRemoteMessage.obj','xercesCommon.obj','CCell.obj','CCellContainer.obj','CremoteSequence.obj','common.obj','commonError.obj','commonAsio.obj','matlabParser.obj', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib), getWinArg('asio', 'include', lib), getWinArg('boost', 'include', lib), getWinArg('boost', 'lib', lib), ['-llibboost_date_time' lib.boostPostfix],['-llibboost_regex' lib.boostPostfix])
            mex  ('-D_WIN32_WINNT=0x0501', '-DXML_LIBRARY', '-g', '../src/remoteGetHelixData.cpp', 'messageCommon.obj','CRemoteMessage.obj','xercesCommon.obj','CCell.obj','CCellContainer.obj','CremoteSequence.obj','common.obj','commonError.obj','commonAsio.obj','matlabParser.obj', getWinArg('xerces', 'lib', lib), '-lxerces-c_static_2_MD', getWinArg('xerces', 'include', lib), getWinArg('asio', 'include', lib), getWinArg('boost', 'include', lib), getWinArg('boost', 'lib', lib), ['-llibboost_date_time' lib.boostPostfix],['-llibboost_regex' lib.boostPostfix])
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
end

% return a string with the type (-L or -I) and the correct path to the librairie
% librairieName: 
%       - xerces
%       - boost
%       - asio
%
% fileType:
%       - lib
%       - include

function pathArg = getWinArg (librairieName, fileType, lib)

	switch computer
        case 'PCWIN64'
            programFiles = getenv('PROGRAMFILES(X86)');
        otherwise 
            programFiles = getenv('PROGRAMFILES');
    end
    
    switch lib.target
        case 'programfiles'
            lPath = strtrim(programFiles );
        case 'relative'
            lPath = strtrim('..\..\dependancies' );
        otherwise
            error ('Windows method path not defined');
    end
    
    switch librairieName
        case 'xerces'
            
            switch computer
                case 'PCWIN64'
                     pathArg = [lPath '\xerces-c_2_8_0-x86_64-windows-vc_8_0'];
                otherwise 
                    pathArg = [lPath '\xerces-c_2_8_0-x86-windows-vc_8_0'];
            end
            switch fileType
                case 'lib'
                   pathArg = ['-L' pathArg  '\lib'];
                case 'include'
                   pathArg = ['-I' pathArg  '\include'];
                otherwise
                    error (['unknown type asked: ' fileType]);
            end
        case 'boost'
           pathArg = [lPath '\boost\boost_1_' lib.boostVersion];
            switch fileType
                case 'lib'
                    switch computer
                        case 'PCWIN64'
                             pathArg = ['-L' pathArg '\lib64'];
                        otherwise 
                             pathArg = ['-L' pathArg '\lib32'];
                    end
                case 'include'
                   pathArg = ['-I' pathArg];
                otherwise
                    error (['unknown type asked: ' fileType]);
            end
        case 'asio'
           pathArg = [lPath '\asio-1.4.1'];
            switch fileType
                case 'include'
                   pathArg = ['-I' pathArg '\include'];
                otherwise
                    error (['unknown type asked: ' fileType]);
            end
        otherwise
            error (['unknown librarie asked : ' librairieName]);
    end
end
