function check = installVineCPP
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
Path = Path(1:end-length(mfilename)-1);

install = questdlg(['Please make sure that the script ''installVineCPP.m'' is still in the main folder of the toolbox that should be installed. The toolbox VineCPP will be installed in the path ' Path ', which includes adding this path to the (non-permanent) search path of your MATLAB installation.'],'Installation of the VineCPP toolbox','Continue','Cancel','Cancel');

if strcmp(install,'Continue')
    check = 1;
    
    %% Add the path of the toolbox to the MATLAB search path
    addpath(Path)
    
    
    %% Ask for permanent path setting
    button = questdlg('Do you want to add the VineCPP toolbox permanently to your search path?','Installation of the VineCPP toolbox','Yes','No','No');
    
    if strcmp(button,'Yes')
        try
            savepath
        catch ME
            warning('VineCPPToolbox:Installation',['There was a problem saving the path permanently. The error was: \n' ME.message])
            check = 0;
        end
    end
    
    %% Writing absolute search paths into the file PathToBoundsAndSeed.hpp
    fid = fopen([Path '/private/PathToBoundsAndSeed.hpp'], 'w');
    fprintf(fid, ['#define PathBounds "' Path '/private/bounds.txt" \n#define PathSeed "' Path '/private/Seed.dat"'],char);
    fclose(fid);
    
    button = questdlg('Should the C++ code be compiled with openmp to enable parallel computing?','Installation of the VineCPP toolbox','Yes','No','No');
    
    if strcmp(button,'Yes')
        %% Compiling C++-Code
        CPP_files = {'PCAIC.cpp VineCPP_helper.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp PairCopulaAIC.cpp -lnlopt',...
            'PCCDF.cpp VineCPP_helper.cpp PairCopulaCDF.cpp mvtdstpack.f -lstdc++',...
            'PCFit.cpp VineCPP_helper.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp -lnlopt',...
            'PCHfun.cpp VineCPP_helper.cpp PairCopulaHfun.cpp',...
            'PCVfun.cpp VineCPP_helper.cpp PairCopulaHfun.cpp',...
            'PCIndepTest.cpp VineCPP_helper.cpp PairCopulaIndepTest.cpp SDtau.cpp',...
            'PCInvHfun.cpp VineCPP_helper.cpp PairCopulaInvHfun.cpp PairCopulaHfun.cpp',...
            'PCInvVfun.cpp VineCPP_helper.cpp PairCopulaInvHfun.cpp PairCopulaHfun.cpp',...
            'PCNegLL.cpp VineCPP_helper.cpp PairCopulaNegLL.cpp',...
            'PCPDF.cpp VineCPP_helper.cpp PairCopulaPDF.cpp',...
            'PCSelect.cpp VineCPP_helper.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp PairCopulaAIC.cpp PairCopulaSelect.cpp PairCopulaIndepTest.cpp SDtau.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -lnlopt',...
            'VineFit.cpp VineCPP_helper.cpp PairCopulaHfun.cpp VineCopulaFit.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp VineCopulaNegLL.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -lnlopt',...
            'VineFitSeq.cpp VineCPP_helper.cpp PairCopulaHfun.cpp VineCopulaFit.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp VineCopulaNegLL.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -lnlopt',...
            'VineGetPseudoObs.cpp VineCPP_helper.cpp PairCopulaHfun.cpp VineCopulaGetPseudoObs.cpp',...
            'VineNegLL.cpp VineCPP_helper.cpp PairCopulaHfun.cpp VineCopulaNegLL.cpp PairCopulaNegLL.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"',...
            'VineRand.cpp VineCPP_helper.cpp PairCopulaHfun.cpp PairCopulaInvHfun.cpp VineCopulaRand.cpp',...
            'VineStructureSelect.cpp VineCopulaStructureSelect.cpp VineCPP_helper.cpp PairCopulaHfun.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp PairCopulaSelect.cpp PairCopulaAIC.cpp PairCopulaIndepTest.cpp SDtau.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -lnlopt',...
            'CvMTestStatCPP.cpp','RandUniform.cpp','RandNormal.cpp','SetSeed.cpp'};
        try
            OldCd = cd;
            cd([Path '/private'])
            for i=1:size(CPP_files,2);
                eval(['mex ' CPP_files{i}])
            end
            cd(OldCd)
        catch ME
            warning('VineCPPToolbox:CppCompilation',['There was a problem compiling the C++ code. The error was: \n' ME.message])
            check = 0;
        end
    else
        %% Compiling C++-Code
        CPP_files = {'PCAIC.cpp VineCPP_helper.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp PairCopulaAIC.cpp -lnlopt',...
            'PCCDF.cpp VineCPP_helper.cpp PairCopulaCDF.cpp mvtdstpack.f -lstdc++',...
            'PCFit.cpp VineCPP_helper.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp -lnlopt',...
            'PCHfun.cpp VineCPP_helper.cpp PairCopulaHfun.cpp',...
            'PCVfun.cpp VineCPP_helper.cpp PairCopulaHfun.cpp',...
            'PCIndepTest.cpp VineCPP_helper.cpp PairCopulaIndepTest.cpp SDtau.cpp',...
            'PCInvHfun.cpp VineCPP_helper.cpp PairCopulaInvHfun.cpp PairCopulaHfun.cpp',...
            'PCInvVfun.cpp VineCPP_helper.cpp PairCopulaInvHfun.cpp PairCopulaHfun.cpp',...
            'PCNegLL.cpp VineCPP_helper.cpp PairCopulaNegLL.cpp',...
            'PCPDF.cpp VineCPP_helper.cpp PairCopulaPDF.cpp',...
            'PCSelect.cpp VineCPP_helper.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp PairCopulaAIC.cpp PairCopulaSelect.cpp PairCopulaIndepTest.cpp SDtau.cpp -lnlopt',...
            'VineFit.cpp VineCPP_helper.cpp PairCopulaHfun.cpp VineCopulaFit.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp VineCopulaNegLL.cpp -lnlopt',...
            'VineFitSeq.cpp VineCPP_helper.cpp PairCopulaHfun.cpp VineCopulaFit.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp VineCopulaNegLL.cpp -lnlopt',...
            'VineGetPseudoObs.cpp VineCPP_helper.cpp PairCopulaHfun.cpp VineCopulaGetPseudoObs.cpp',...
            'VineNegLL.cpp VineCPP_helper.cpp PairCopulaHfun.cpp VineCopulaNegLL.cpp PairCopulaNegLL.cpp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"',...
            'VineRand.cpp VineCPP_helper.cpp PairCopulaHfun.cpp PairCopulaInvHfun.cpp VineCopulaRand.cpp',...
            'VineStructureSelect.cpp VineCopulaStructureSelect.cpp VineCPP_helper.cpp PairCopulaHfun.cpp PairCopulaFit.cpp PairCopulaNegLL.cpp PairCopulaSelect.cpp PairCopulaAIC.cpp PairCopulaIndepTest.cpp SDtau.cpp -lnlopt',...
            'CvMTestStatCPP.cpp','RandUniform.cpp','RandNormal.cpp','SetSeed.cpp'};
        try
            OldCd = cd;
            cd([Path '/private'])
            for i=1:size(CPP_files,2);
                eval(['mex ' CPP_files{i}])
            end
            cd(OldCd)
        catch ME
            warning('VineCPPToolbox:CppCompilation',['There was a problem compiling the C++ code. The error was: \n' ME.message])
            check = 0;
        end
    end
        
    %% Setting a primary seed
    SetSeedVineCPP;
    
    warning('off','VineCPP:CopulaDataOnBounds')
    
    if check == 1
        msgbox('The installation of the VineCPP toolbox is completed.','Installation of the VineCPP toolbox','custom',imread([Path '/private/Symbol.tif']));
    else
        msgbox('The installation of the VineCPP toolbox failed.', 'Error','error');
    end
    
    
end

end


