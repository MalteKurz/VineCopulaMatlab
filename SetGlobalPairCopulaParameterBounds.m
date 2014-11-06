function SetGlobalPairCopulaParameterBounds
%SETGLOBALPAIRCOPULAPARAMETERBOUNDS Setting the global (i.e., for the whole VineCPP toolbox) pair-copula parameter bounds, which are used for estimation and consistency checks
% Purpose
%        The GUI can be used to set the global (i.e., for the whole VineCPP
%        toolbox) pair-copula parameter bounds, which are used for
%        estimation and consistency checks. Possible copula families (the
%        default values for the parameter bounds are given in parenthesis,
%        where the ordering is as follows (lb1, ub1; lb2, ub2; lb3, ub3)):
%           0   Indep
%           1   AMH (-1, 1)
%           2   AsymFGM (0, 1)
%           3   BB1 (0, 6; 1, 6)
%           4   BB6 (1, 6; 1, 6)
%           5   BB7 (1, 6; 0.001, 6)
%           6   BB8 (1, 6; 0.001, 1)
%           7   Clayton (0, Inf)
%           8   FGM (-1, 1)
%           9   Frank (-30, 30)
%           10  Gaussian (-0.999, 0.999)
%           11  Gumbel (1, Inf)
%           12  IteratedFGM (-1, 1; -1, 1)
%           13  Joe (1, Inf)
%           14  PartialFrank (0, 30)
%           15  Plackett (0.001, Inf)
%           16  Tawn1 (1.001, 20; 0.001, 0.999)
%           17  Tawn2 (1.001, 20; 0.001, 0.999)
%           18  Tawn (1.0001, 20; 0.001, 0.999; 0.001, 0.999)
%           19  t (-0.999, 0.999; 1, 30)
%
%
% Usage
%               SetGlobalPairCopulaParameterBounds
%
%
%
% Author: Malte Kurz

Path = mfilename('fullpath');
Path = Path(1:end-length(mfilename)-1);

bounds = dlmread([Path '/private/bounds.txt']);
bounds_extreme = dlmread([Path '/private/bounds_extreme.txt']);
bounds_default = dlmread([Path '/private/bounds_default.txt']);
bounds(1,:) = [];
bounds_extreme(1,:) = [];
bounds_default(1,:) = [];

bounds_old = bounds;

% Create the GUI and hide it until the construction is finished.
f = figure('Visible','off','Position',[300,300,300,800]);

% Construct the components:

% Pair-copula families
htextFamily = uicontrol('Style','text','String','Select a pair-copula family',...
    'Position',[20,760,260,25],...
    'BackgroundColor',[.505,.639,.505],...
    'FontSize',12);
hpopupFamily = uicontrol('Style','popupmenu',...
    'String',{'AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'},...
    'Position',[20,720,260,25],...
    'Callback',{@PairCopulaFamily_Callback},...
    'FontSize',12);

% Bounds for the first copula parameters
htextLb1 = uicontrol('Style','text','String','Lower bound of the first parameter',...
    'Position',[20,650,260,25],...
    'BackgroundColor',[.505,.639,.505],...
    'FontSize',12);
heditLb1 = uicontrol('Style','edit',...
    'Position',[20,610,260,25],...
    'Callback',{@PairCopulaLb1_Callback},...
    'FontSize',12,...
    'String',num2str(bounds(1,1)));
htextUb1 = uicontrol('Style','text','String','Upper bound of the first parameter',...
    'Position',[20,570,260,25],...
    'BackgroundColor',[.505,.639,.505],...
    'FontSize',12);
heditUb1 = uicontrol('Style','edit',...
    'Position',[20,530,260,25],...
    'Callback',{@PairCopulaUb1_Callback},...
    'FontSize',12,...
    'String',num2str(bounds(1,2)));

% Bounds for the second copula parameters
htextLb2 = uicontrol('Style','text','String','Lower bound of the second parameter',...
    'Position',[20,480,260,25],...
    'Visible','off',...
    'BackgroundColor',[.505,.639,.505],...
    'FontSize',12);
heditLb2 = uicontrol('Style','edit',...
    'String',{''},...
    'Position',[20,440,260,25],...
    'Callback',{@PairCopulaLb3_Callback},...
    'Visible','off',...
    'FontSize',12);
htextUb2 = uicontrol('Style','text','String','Upper bound of the second parameter',...
    'Position',[20,400,260,25],...
    'Visible','off',...
    'BackgroundColor',[.505,.639,.505],...
    'FontSize',12);
heditUb2 = uicontrol('Style','edit',...
    'String',{''},...
    'Position',[20,360,260,25],...
    'Callback',{@PairCopulaUb3_Callback},...
    'Visible','off',...
    'FontSize',12);

% Bounds for the third copula parameters
htextLb3 = uicontrol('Style','text','String','Lower bound of the third parameter',...
    'Position',[20,310,260,25],...
    'Visible','off',...
    'BackgroundColor',[.505,.639,.505],...
    'FontSize',12);
heditLb3 = uicontrol('Style','edit',...
    'String',{''},...
    'Position',[20,270,260,25],...
    'Callback',{@PairCopulaLb3_Callback},...
    'Visible','off',...
    'FontSize',12);
htextUb3 = uicontrol('Style','text','String','Upper bound of the third parameter',...
    'Position',[20,230,260,25],...
    'Visible','off',...
    'BackgroundColor',[.505,.639,.505],...
    'FontSize',12);
heditUb3 = uicontrol('Style','edit',...
    'String',{''},...
    'Position',[20,190,260,25],...
    'Callback',{@PairCopulaUb3_Callback},...
    'Visible','off',...
    'FontSize',12);

% Buttons
hbuttonDefault = uicontrol('Style','pushbutton','String','Default',...
    'Position',[25,60,70,30],...
    'Callback',{@Default_Callback});
hbuttonSave = uicontrol('Style','pushbutton','String','Save',...
    'Position',[115,60,70,30],...
    'Callback',{@Save_Callback});
hbuttonClose = uicontrol('Style','pushbutton','String','Close',...
    'Position',[205,60,70,30],...
    'Callback',{@Close_Callback});
align([htextFamily,hpopupFamily,htextLb1,heditLb1,htextUb1,heditUb1,htextLb2,heditLb2,htextUb2,heditUb2,htextLb3,heditLb3,htextUb3,heditUb3],'Center','None');

set([htextFamily,hpopupFamily,htextLb1,heditLb1,htextUb1,heditUb1,htextLb2,heditLb2,htextUb2,heditUb2,htextLb3,heditLb3,htextUb3,heditUb3,hbuttonDefault,hbuttonSave,hbuttonClose],...
    'Units','normalized');
% Name
set(f,'Name','Setting the global parameter bounds for the pair-copulas')
% Move the GUI to the center of the screen.
movegui(f,'center')
% Make the GUI visible.
set(f,'Visible','on');

% Set the callback, when the GUI is closed by the user.
set(f,'CloseRequestFcn',@Close_Callback)

% Callback for the Pair-copula families
    function PairCopulaFamily_Callback(source,eventdata)
        
        % Obtain the new data.
        str = get(source, 'String');
        val = get(source,  'Value');
        
        % Apply changes (especially set the visibilities of the different
        % windows and get the current values for the bounds of the chosen
        % copula).
        switch str{val};
            case {'AMH','AsymFGM','Clayton','FGM','Frank','Gaussian','Gumbel','Joe','PartialFrank','Plackett'}
                set(heditLb1,'String',num2str(bounds(val,1)))
                set(heditUb1,'String',num2str(bounds(val,2)))
                set(htextLb2,'Visible','off')
                set(heditLb2,'Visible','off')
                set(htextUb2,'Visible','off')
                set(heditUb2,'Visible','off')
                set(htextLb3,'Visible','off')
                set(heditLb3,'Visible','off')
                set(htextUb3,'Visible','off')
                set(heditUb3,'Visible','off')
            case {'BB1','BB6','BB7','BB8','IteratedFGM','Tawn1','Tawn2','t'}
                set(heditLb1,'String',num2str(bounds(val,1)))
                set(heditUb1,'String',num2str(bounds(val,2)))
                set(htextLb2,'Visible','on')
                set(heditLb2,'Visible','on','String',num2str(bounds(val,3)))
                set(htextUb2,'Visible','on')
                set(heditUb2,'Visible','on','String',num2str(bounds(val,4)))
                set(htextLb3,'Visible','off')
                set(heditLb3,'Visible','off')
                set(htextUb3,'Visible','off')
                set(heditUb3,'Visible','off')
            case {'Tawn'}
                set(heditLb1,'String',num2str(bounds(val,1)))
                set(heditUb1,'String',num2str(bounds(val,2)))
                set(htextLb2,'Visible','on')
                set(heditLb2,'Visible','on','String',num2str(bounds(val,3)))
                set(htextUb2,'Visible','on')
                set(heditUb2,'Visible','on','String',num2str(bounds(val,4)))
                set(htextLb3,'Visible','on')
                set(heditLb3,'Visible','on','String',num2str(bounds(val,5)))
                set(htextUb3,'Visible','on')
                set(heditUb3,'Visible','on','String',num2str(bounds(val,6)))
        end
    end

% Callback for the lower bounds for the first copula parameters
    function PairCopulaLb1_Callback(source,eventdata)
        
        % Obtain the new data and the current selections.
        val = get(hpopupFamily,'Value');
        BoundInput = get(source,  'String');
        
        % Apply changes.
        NewBound = str2double(BoundInput);
        if NewBound >= bounds_extreme(val,1) && NewBound < bounds_extreme(val,2) && NewBound < bounds(val,2)
            bounds(val,1) = NewBound;
        else
            h = errordlg(['Wrong input. The parameter bounds have to form a non-empty subset of the intervall (' num2str(bounds_extreme(val,1)) ' , ' num2str(bounds_extreme(val,2)) ').']);
            a = get(h,'Children');
            set(get(a(2),'Children'),'FontSize',13)
            set(h,'Children',a)
            set(heditLb1,'String',num2str(bounds(val,1)))
        end
    end

% Callback for the upper bounds for the first copula parameters
    function PairCopulaUb1_Callback(source,eventdata)
        
        % Obtain the new data and the current selections.
        val = get(hpopupFamily,'Value');
        BoundInput = get(source,  'String');
        
        % Apply changes.
        NewBound = str2double(BoundInput);
        if NewBound > bounds_extreme(val,1) && NewBound <= bounds_extreme(val,2) && NewBound > bounds(val,1)
            bounds(val,2) = NewBound;
        else
            h = errordlg(['Wrong input. The parameter bounds have to form a non-empty subset of the intervall (' num2str(bounds_extreme(val,1)) ' , ' num2str(bounds_extreme(val,2)) ').']);
            a = get(h,'Children');
            set(get(a(2),'Children'),'FontSize',13)
            set(h,'Children',a)
            set(heditUb1,'String',num2str(bounds(val,2)))
        end
    end

% Callback for the lower bounds for the second copula parameters
    function PairCopulaLb2_Callback(source,eventdata)
        
        % Obtain the new data and the current selections.
        val = get(hpopupFamily,'Value');
        BoundInput = get(source,  'String');
        
        % Apply changes.
        NewBound = str2double(BoundInput);
        if NewBound >= bounds_extreme(val,3) && NewBound < bounds_extreme(val,4) && NewBound < bounds(val,4)
            bounds(val,3) = NewBound;
        else
            h = errordlg(['Wrong input. The parameter bounds have to form a non-empty subset of the intervall (' num2str(bounds_extreme(val,3)) ' , ' num2str(bounds_extreme(val,4)) ').']);
            a = get(h,'Children');
            set(get(a(2),'Children'),'FontSize',13)
            set(h,'Children',a)
            set(heditLb2,'String',num2str(bounds(val,3)))
        end
    end

% Callback for the upper bounds for the second copula parameters
    function PairCopulaUb2_Callback(source,eventdata)
        
        % Obtain the new data and the current selections.
        val = get(hpopupFamily,'Value');
        BoundInput = get(source,  'String');
        
        % Apply changes.
        NewBound = str2double(BoundInput);
        if NewBound > bounds_extreme(val,3) && NewBound <= bounds_extreme(val,4) && NewBound > bounds(val,3)
            bounds(val,4) = NewBound;
        else
            h = errordlg(['Wrong input. The parameter bounds have to form a non-empty subset of the intervall (' num2str(bounds_extreme(val,3)) ' , ' num2str(bounds_extreme(val,4)) ').']);
            a = get(h,'Children');
            set(get(a(2),'Children'),'FontSize',13)
            set(h,'Children',a)
            set(heditUb2,'String',num2str(bounds(val,4)))
        end
    end

% Callback for the lower bounds for the third copula parameters
    function PairCopulaLb3_Callback(source,eventdata)
        
        % Obtain the new data and the current selections.
        val = get(hpopupFamily,'Value');
        BoundInput = get(source,  'String');
        
        % Apply changes.
        NewBound = str2double(BoundInput);
        if NewBound >= bounds_extreme(val,5) && NewBound < bounds_extreme(val,6) && NewBound < bounds(val,6)
            bounds(val,5) = NewBound;
        else
            h = errordlg(['Wrong input. The parameter bounds have to form a non-empty subset of the intervall (' num2str(bounds_extreme(val,5)) ' , ' num2str(bounds_extreme(val,6)) ').']);
            a = get(h,'Children');
            set(get(a(2),'Children'),'FontSize',13)
            set(h,'Children',a)
            set(heditLb3,'String',num2str(bounds(val,5)))
        end
    end

% Callback for the upper bounds for the third copula parameters
    function PairCopulaUb3_Callback(source,eventdata)
        
        % Obtain the new data and the current selections.
        val = get(hpopupFamily,'Value');
        BoundInput = get(source,  'String');
        
        % Apply changes.
        NewBound = str2double(BoundInput);
        if NewBound > bounds_extreme(val,5) && NewBound <= bounds_extreme(val,6) && NewBound > bounds(val,5)
            bounds(val,6) = NewBound;
        else
            h = errordlg(['Wrong input. The parameter bounds have to form a non-empty subset of the intervall (' num2str(bounds_extreme(val,5)) ' , ' num2str(bounds_extreme(val,6)) ').']);
            a = get(h,'Children');
            set(get(a(2),'Children'),'FontSize',13)
            set(h,'Children',a)
            set(heditUb3,'String',num2str(bounds(val,6)))
        end
    end

% Callback for the default values button
    function Default_Callback(source,eventdata)
        fs=get(0,'DefaultUicontrolFontSize');
        set(0,'DefaultUicontrolFontSize',13);
        button = questdlg('Which parameter bounds should be set back to default values?','Default parameter bounds','Default for the current pair-copula','Default for all pair-copulas','Cancel','Default for the current pair-copula');
        set(0,'DefaultUicontrolFontSize',fs);
        if strcmp(button,'Default for the current pair-copula')
            button = 1;
        elseif strcmp(button,'Default for all pair-copulas')
            button = 2;
            bounds = bounds_default;
        else
            button = 3;
        end
        
        % Obtain the current selections.
        str = get(hpopupFamily, 'String');
        val = get(hpopupFamily,'Value');
        
        % Apply changes.
        switch str{val};
            case {'AMH','AsymFGM','Clayton','FGM','Frank','Gaussian','Gumbel','Joe','PartialFrank','Plackett'}
                if button == 1
                    bounds(val,1) = bounds_default(val,1);
                    bounds(val,2) = bounds_default(val,2);
                end
                set(heditLb1,'String',num2str(bounds(val,1)))
                set(heditUb1,'String',num2str(bounds(val,2)))
            case {'BB1','BB6','BB7','BB8','IteratedFGM','Tawn1','Tawn2','t'}
                if button == 1
                    bounds(val,1) = bounds_default(val,1);
                    bounds(val,2) = bounds_default(val,2);
                    bounds(val,3) = bounds_default(val,3);
                    bounds(val,4) = bounds_default(val,4);
                end
                set(heditLb1,'String',num2str(bounds(val,1)))
                set(heditUb1,'String',num2str(bounds(val,2)))
                set(heditLb2,'Visible','on','String',num2str(bounds(val,3)))
                set(heditUb2,'Visible','on','String',num2str(bounds(val,4)))
            case {'Tawn'}
                if button == 1
                    bounds(val,1) = bounds_default(val,1);
                    bounds(val,2) = bounds_default(val,2);
                    bounds(val,3) = bounds_default(val,3);
                    bounds(val,4) = bounds_default(val,4);
                    bounds(val,5) = bounds_default(val,5);
                    bounds(val,6) = bounds_default(val,6);
                end
                set(heditLb1,'String',num2str(bounds(val,1)))
                set(heditUb1,'String',num2str(bounds(val,2)))
                set(heditLb2,'Visible','on','String',num2str(bounds(val,3)))
                set(heditUb2,'Visible','on','String',num2str(bounds(val,4)))
                set(heditLb3,'Visible','on','String',num2str(bounds(val,5)))
                set(heditUb3,'Visible','on','String',num2str(bounds(val,6)))
        end
    end

% Callback for the save button.
    function Save_Callback(source,eventdata)
        fs=get(0,'DefaultUicontrolFontSize');
        set(0,'DefaultUicontrolFontSize',13);
        button = questdlg('Save the current parameter bounds as global bounds for the whole TestOnSimplifiedVineCopulas toolbox?','Parameter bounds','Yes','No','No');
        set(0,'DefaultUicontrolFontSize',fs);
        if strcmp(button,'Yes')
            bounds = [zeros(1,6);bounds];
            dlmwrite([Path '/private/bounds.txt'],bounds,'\t')
            bounds(1,:) = [];
            bounds_old = bounds;
        end
    end

% Callback for the close button.
    function Close_Callback(source,eventdata)
        if sum(sum(bounds~=bounds_old)) > 0
            fs=get(0,'DefaultUicontrolFontSize');
            set(0,'DefaultUicontrolFontSize',13);
            button = questdlg('Save the changes in the parameter bounds as global bounds for the whole TestOnSimplifiedVineCopulas toolbox?','Parameter bounds','Yes','No','No');
            set(0,'DefaultUicontrolFontSize',fs);
            if strcmp(button,'Yes')
            bounds = [zeros(1,6);bounds];
            dlmwrite([Path '/private/bounds.txt'],bounds,'\t')
            end
        end
        delete(gcf)
    end

end
