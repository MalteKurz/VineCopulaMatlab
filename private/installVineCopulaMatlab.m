function check = installVineCopulaMatlab
%INSTALLVINECPP Setting up the VineCPP toolbox
% Purpose
%        The function adds the folder VineCPP with its subfolders to the
%        current search path of MATLAB. One can optionally apply this
%        change in the search path permanently.
%        Furthermore, the C++ files are compiled.
%
%
% Usage
%               installVineCPP
%
%
%
% Author: Malte Kurz

Path = mfilename('fullpath');
Path = Path(1:end-length(mfilename)-8);
cd(Path);

check = 1;

%% Add the path of the toolbox to the MATLAB search path
addpath(Path)

%% Linking the C++ shared library VineCopulaCPP for usage from MATLAB
PathRoot = matlabroot;
system(['cd src; mv libVineCPP.so.1.0 ' PathRoot '/bin/glnxa64']);
system(['ln -sf ' PathRoot '/bin/glnxa64/libVineCPP.so.1.0 ' PathRoot '/bin/glnxa64/libVineCPP.so.1']);
system(['ln -sf ' PathRoot '/bin/glnxa64/libVineCPP.so.1.0 ' PathRoot '/bin/glnxa64/libVineCPP.so']);

%% Compiling C++-Code
CPP_files = {'PCAIC.cpp',...
    'PCCDF.cpp',...
    'PCFit.cpp',...
    'PCHfun.cpp',...
    'PCVfun.cpp',...
    'PCIndepTest.cpp',...
    'PCInvHfun.cpp',...
    'PCInvVfun.cpp',...
    'PCNegLL.cpp',...
    'PCPDF.cpp',...
    'PCRand.cpp',...
    'PCSelect.cpp',...
    'VineFit.cpp',...
    'VineFitSeq.cpp',...
    'VineGetPseudoObs.cpp',...
    'VineNegLL.cpp',...
    'VineRand.cpp',...
    'VineStructureSelect.cpp',...
    'CvMTestStatCPP.cpp',...
    'RandUniform.cpp',...
    'RandNormal.cpp',...
    'SetSeed.cpp',...
    'FastKendallTau.cpp'...
    };
try
    OldCd = cd;
    cd([Path '/private'])
    for i=1:size(CPP_files,2);
        eval(['mex ' CPP_files{i} ' -I' Path '/src  CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -lVineCPP -lnlopt -lgfortran'])
    end
    cd(OldCd)
catch ME
    warning('VineCPPToolbox:CppCompilation',['There was a problem compiling the C++ code. The error was: \n' ME.message])
    check = 0;
end

%% Setting a primary seed
%SetSeedVineCPP;

warning('off','VineCPP:CopulaDataOnBounds')

if check == 1
    disp('The installation of the VineCopulaMatlab toolbox is completed.');
else
    error('The installation of the VineCopulaMatlab toolbox failed.');
end


end


