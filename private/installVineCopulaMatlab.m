function check = installVineCopulaMatlab(libdir,IncludeDirs,LibDirs,varargin)
%INSTALLVINECPP Setting up the VineCPP toolbox
% Purpose
%        The function adds the folder VineCPP with its subfolders to the
%        current search path of MATLAB. One can optionally apply this
%        change in the search path permanently.
%        Furthermore, the C++ files are compiled.
%
%
%
% Author: Malte Kurz

Path = mfilename('fullpath');
Path = Path(1:end-length(mfilename)-8);
cd(Path);

PathRoot = matlabroot;

check = 1;

%% Add the path of the toolbox to the MATLAB search path
addpath(Path)

%% Compiling C++-Code
CPP_files = {'VineCopulaMatlab.cpp'};
try
    OldCd = cd;
    cd([Path '/private'])
    if isunix %% Linux system
        %% Linking the C++ shared library VineCopulaCPP for usage from MATLAB
        PathRoot = matlabroot;
        system(['ln -sf ' libdir '/libVineCopulaCPP.so.1.0 ' PathRoot '/bin/glnxa64/libVineCopulaCPP.so']);
        
        for i=1:size(CPP_files,2);
            eval(['mex ' CPP_files{i} ' CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" -lVineCopulaCPP -lnlopt -lgfortran'])
        end
    elseif ispc % Windows system
        %% Linking the C++ shared library VineCopulaCPP for usage from MATLAB
        system(['cp "' libdir '/libVineCopulaCPP.dll" "."']);
        if nargin==5
            setenv('MINGWROOT',varargin{1})
            eval(['mex -setup:' varargin{2}])
        end
        I = findstr('-I/',IncludeDirs);
        I = [I length(IncludeDirs)+1];
        IncludeDirsTilde = '';
        for i=1:length(I)-1
            IncludeDirsTilde = [IncludeDirsTilde IncludeDirs(I(i):I(i)+1) IncludeDirs(I(i)+3) ':' IncludeDirs(I(i)+4:I(i+1)-1)];
        end
        I = findstr('-L/',LibDirs);
        I = [I length(LibDirs)+1];
        LibDirsTilde = '';
        for i=1:length(I)-1
            LibDirsTilde = [LibDirsTilde LibDirs(I(i):I(i)+1) LibDirs(I(i)+3) ':' LibDirs(I(i)+4:I(i+1)-1)];
        end
        IncludeDirs = IncludeDirsTilde;
        LibDirs = LibDirsTilde;
        for i=1:size(CPP_files,2);
            eval(['mex ' CPP_files{i} ' ' IncludeDirs ' ' LibDirs])
        end
    end
    cd(OldCd)
catch ME
    warning('VineCopulaMatlabToolbox:CppCompilation',['There was a problem compiling the C++ code. The error was: \n' ME.message])
    check = 0;
end

%% Setting a primary seed
%SetSeedVineCPP;

warning('off','VineCPP:CopulaDataOnBounds')

if check == 1
    disp('The installation of the VineCopulaMatlab toolbox is completed.');
else
    disp('The installation of the VineCopulaMatlab toolbox failed.');
    pause
end


end


