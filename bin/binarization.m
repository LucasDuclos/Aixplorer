% Function to binarize remote
function binarization()

% ============================================================================ %
% ============================================================================ %

% Initialization
BinPath = pwd;
Slash   = strfind(BinPath, '\');
if ( ~isempty(Slash) )
    BinPath(Slash) = '/';
end
SrcPath = BinPath;
SrcPath(end-2:end) = 'src';

% ============================================================================ %
% ============================================================================ %

% List packages and classes
FileList = dir(SrcPath); cd(SrcPath);
Packages = [];

% Loop over all packages
for k = 1 : length(FileList)
    FileList(k).name

    % Build packages list
    if strcmpi(FileList(k).name(1), '+') ...
            && ~strcmpi(FileList(k).name, '+scripts') ...
            && ~strcmpi(FileList(k).name, '+tools') ... TODO: put +tools/* in @ classes
            && ~strcmpi(FileList(k).name, '+test')
        Packages(end+1).SrcPath = [SrcPath '/' FileList(k).name '/'];
        Packages(end).BinPath   = [BinPath '/' FileList(k).name '/'];
        Packages(end).Info      = meta.package.fromName(FileList(k).name(2:end));
    end
    
end

% ============================================================================ %

% Binarize all methods and classes
while ( ~isempty(Packages) )
% for k = 1 : length(Packages)
    
    disp(['Binarize package: ' Packages(1).SrcPath]);

    % Binarize classes and class methods
    if ( ~isempty(Packages(1).Info.Classes) )
        for l = 1 : length(Packages(1).Info.Classes)
            
            % Extract class information
            Idx       = strfind(Packages(1).Info.Classes{l}.Name, '.');
            ClassName = Packages(1).Info.Classes{l}.Name(Idx(end)+1:end);
            
            % Binarize class
            pcode([Packages(1).SrcPath '@' ClassName], '-inplace');
            mkdir( [ Packages(1).BinPath '@' ClassName '/' ] );
            movefile( ...
                [Packages(1).SrcPath '@' ClassName '/*.p'], ...
                [Packages(1).BinPath '@' ClassName '/']);

        end
    end
    
% ============================================================================ %

    % Binarize functions
    if ( ~isempty(Packages(1).Info.Functions) )
        for l = 1 : length(Packages(1).Info.Functions)
            
            FuncName = ...
                [Packages(1).SrcPath Packages(1).Info.Functions{l}.Name '.m'];
            if ( exist(FuncName, 'file') == 2 )
                pcode(FuncName, '-inplace');
                
                movefile( ...
                    [Packages(1).SrcPath Packages(1).Info.Functions{l}.Name '.p'], ...
                    [Packages(1).BinPath Packages(1).Info.Functions{l}.Name '.p']);
            end
            
        end
    end

% ============================================================================ %

    % Additional packages
    if ( ~isempty(Packages(1).Info.Packages) )
        for l = 1 : length(Packages(1).Info.Packages)
            
            % Extract package information
            Idx      = strfind(Packages(1).Info.Packages{l}.Name, '.');
            PackName = Packages(1).Info.Packages{l}.Name(Idx(end)+1:end);
            
            % Add package
            Packages(end+1).SrcPath = [Packages(1).SrcPath '+' PackName '/'];
            Packages(end).BinPath   = [Packages(1).BinPath '+' PackName '/'];
            Packages(end).Info      = ...
                meta.package.fromName(Packages(1).Info.Packages{l}.Name);
        end
    end

% ============================================================================ %

    % Binarize config files
    if ( strcmpi(Packages(1).Info.Name, 'system.config') )
        
        % Backup hardware and probe files
        PackPath = Packages(1).SrcPath;
        Idx      = strfind(PackPath, '/');
        HwName   = [PackPath(1:Idx(end-1)) '@hardware/hardware.'];
        PrName   = [PackPath(1:Idx(end-1)) '@probe/probe.'];
%         copyfile([HwName 'm'], [HwName 'bak']);
%         copyfile([PrName 'm'], [PrName 'bak']);
        
        % Binarize hardware and probe files
        FileList = dir(PackPath);
        for k = 1 : length(FileList)
            if ( (length(FileList(k).name) > 1) ...
                    && strcmpi(FileList(k).name(end-1:end), '.m') )
                
                % Binarize hardware file
                if ( strcmpi(FileList(k).name(1:3), 'hw_') )
                    
                    FileName = FileList(k).name(1:end-1);
                    copyfile([PackPath FileName 'm'], [HwName 'm']);
                    pcode([HwName 'm'], '-inplace');
                    movefile([HwName 'p'], [PackPath FileName 'p'])
                    
                % Binarize probe file
                elseif ( strcmpi(FileList(k).name(1:3), 'pr_') )
                    
                    FileName = FileList(k).name(1:end-1);
                    copyfile([PackPath FileName 'm'], [PrName 'm']);
                    % pcode([PrName 'm'], '-inplace');
                    % movefile([PrName 'p'], [PackPath FileName 'p'])
                    mkdir( Packages(1).BinPath );
                    copyfile([PackPath FileName 'm'], Packages(1).BinPath);
                    
                end
                
            end
        end
        
        % Restore hardware and probe files
%         movefile([HwName 'bak'], [HwName 'm']);
%         movefile([PrName 'bak'], [PrName 'm']);
        
        % Move binarized hardware files
        movefile([PackPath 'hw_*.p'], Packages(1).BinPath);
        % Move binarized probe files
%         movefile([PackPath 'pr_*.p'], Packages(1).BinPath);
        
    % Nothing to binarize
    elseif ( isempty(Packages(1).Info.Classes) ...
            && isempty(Packages(1).Info.Functions) ...
            && isempty(Packages(1).Info.Packages) )
        
        Packages(1).Info
        disp('Nothing to binarize...');
        
    end

% ============================================================================ %

    % Update the package list
    Packages = Packages(2:end);

end

% ============================================================================ %

% Copy all mex-files of +remote/@tw
Listing = dir([SrcPath '/+remote/@tw/*.mex*']);
for k = 1 : length(Listing)
    copyfile([SrcPath '/+remote/@tw/' Listing(k).name], ...
        [BinPath '/+remote/@tw/' Listing(k).name]);
end

% Copy all mex-files of +processing
Listing = dir([SrcPath '/+processing/*.mex*']);
for k = 1 : length(Listing)
    copyfile([SrcPath '/+processing/' Listing(k).name], ...
        [BinPath '/+processing/' Listing(k).name]);
end

% Copy all scripts files
mkdir( [ BinPath '/+scripts/' ] )
Listing = dir( [ SrcPath '/+scripts/*.m' ] );
for k = 1 : length(Listing)
    copyfile( [ SrcPath '/+scripts/' Listing(k).name ], ...
        [ BinPath '/+scripts/' Listing(k).name ] );
end

% Copy all constants files & delete constant p files
Listing = dir([SrcPath '/+common/@constants/*.m']);
for k = 1 : length(Listing)
    copyfile([SrcPath '/+common/@constants/' Listing(k).name], ...
        [BinPath '/+common/@constants/' Listing(k).name]);
end
Listing = dir([BinPath '/+common/@constants/*.p']);
for k = 1 : length(Listing)
    delete( [ BinPath '/+common/@constants/' Listing(k).name ] );
end

% Make Aixplorer the default hardware
copyfile( [ BinPath '/+system/+config/hw_aixplorer.p' ], [ BinPath '/+system/@hardware/hardware.p' ] )

% Make SL15-4 the default probe
delete(  [ BinPath '/+system/@probe/probe.p' ] )
copyfile( [ BinPath '/+system/+config/pr_SL15-4.m' ], [ BinPath '/+system/@probe/probe.m' ] )

% Copy help file
pcode( [ SrcPath '/legHAL_help.m' ], '-inplace');
movefile( [ SrcPath '/legHAL_help.p' ], [ BinPath '/legHAL_help.p' ] )


% ============================================================================ %

% Finalize binarization
cd(BinPath);

% ============================================================================ %
% ============================================================================ %



end
