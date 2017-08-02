% G_module is used to generate the code associated with the gap constructor G.

% The standard application of the Minimalist Principle for a gap must be modified if the gap
% is followed immediatley by a monochromator (see G_module_mono support function below for details).
% To deal with this dependency, the G_module has been forked with the original module written
% by Mads renamed as G_module_standard and the monochromator exception written by Leland renamed 
% as G_module_mono.

% Modified by Leland Harriger, NCNR Oct. 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% START OF MAIN FUNCTION G_MODULE() %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [McStasStr] = G_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)

% Check if module immediately after the gap is a monochromator
if (index ~= 1 && ~strcmp(globalinfo.modules{globalinfo.modulelist(index-1)}, 'M'))
    % Next module is not a monochromator, do standard Minimalist calculation.
    [McStasStr] = G_module_standard(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults);
else
    % Next module is a monochromator, do monochromator Minimalist calculation.
    [McStasStr] = G_module_mono(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% END OF MAIN FUNCTION G_MODULE() %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [McStasStr] = G_module_standard(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
%M_module adds text to McStas bot to make space between two guide elements
   % add the relevant component to the trace section
   % add the relevant input to the input list
   % add the relevant calculations to the initialize section
   % add the relevant variable names to the declare list
   
   
   
num=num2str(index);
numM1=num2str(index-1);

if (index==last || index ==1) 
    if (index==last)
        disp('ERROR, do not place G at the start of input-string, change closest_element instead')
        disp('ERROR, the first G module will be ignored')
    end
    if (index==1)
        disp('ERROR, do not place G at the end of input-string, change sample_dist instead')
        disp('ERROR, the last G module will be ignored')
    end
else % add a free space to the instrument file

defaultnames=fieldnames(defaults);
% Input not concerning the geometry
% WARNING, the indexes are used in the part extracting data from options
% linp=length(McStasStr.input);
% R0Index=linp+1;
% McStasStr.input{R0Index}=['R0' num];
% if ismember(defaultnames,'R0')
% McStasStr.inputvalue(R0Index)=defaults.R0;
% else    
% McStasStr.inputvalue(R0Index)=0.99;
% end
% 
% QcIndex=linp+2;
% McStasStr.input{QcIndex}=['Qc' num];
% if ismember(defaultnames,'Qc')
% McStasStr.inputvalue(QcIndex)=defaults.Qc;    
% else
% McStasStr.inputvalue(QcIndex)=0.0217;
% end
% 
% alphaIndex=linp+3;
% McStasStr.input{alphaIndex}=['alpha' num];
% if ismember(defaultnames,'alpha')
% McStasStr.inputvalue(alphaIndex)=defaults.alpha;
% else
% McStasStr.inputvalue(alphaIndex)=6.07;
% end
% 
% mIndex=linp+4;
% McStasStr.input{mIndex}=['m' num];
% if ismember(defaultnames,'m')
% McStasStr.inputvalue(mIndex)=defaults.m;
% else
% McStasStr.inputvalue(mIndex)=3;
% end
% 
% WIndex=linp+5;
% McStasStr.input{WIndex}=['W' num];
% if ismember(defaultnames,'W')
% McStasStr.inputvalue(WIndex)=defaults.W;
% else
% McStasStr.inputvalue(WIndex)=0.003;
% end

optimlimits=-ones(1,4);
locked=-ones(1,2); % width,height
% Extracting options from the option cell
for ops=1:length(options)
    for i=1:length(options{ops})
        if strcmp(options{ops}(i),'='); 
            equals=i;
        end
    end
    
    % check for options on reflectivity curve
    if strcmp(options{ops}(1:equals-1),'minStartWidth')
       if optimlimits(1)>0
           optimlimits(1)=max([optimlimits(1) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(1)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'maxStartWidth')
       if optimlimits(2)>0
           optimilimts(2)=min([optimlimits(2) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(2)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'minStartHeight')
       if optimlimits(3)>0
           optimilimts(3)=max([optimlimits(3) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(3)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'maxStartHeight')
       if optimlimits(4)>0
           optimilimts(4)=min([optimlimits(4) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(4)=str2num(options{ops}(equals+1:end));   
       end
    elseif strcmp(options{ops}(1:equals-1),'StartWidth')
       if locked(1)>0
           locked(1)=min([locked(1) str2num(options{ops}(equals+1:end))]);
       else
           locked(1)=str2num(options{ops}(equals+1:end));   
       end
    elseif strcmp(options{ops}(1:equals-1),'StartHeight')
       if locked(2)>0
           locked(2)=min([locked(2) str2num(options{ops}(equals+1:end))]);
       else
           locked(2)=str2num(options{ops}(equals+1:end));   
       end
       % These avoid errors assosicated with options which are treated by
       % the main loop. changed is for when options are removed. 
    elseif (strcmp(options{ops}(1:equals-1),'start') || strcmp(options{ops}(1:equals-1),'length') || strcmp(options{ops}(1:equals-1),'changed') ||...
            strcmp(options{ops}(1:equals-1),'maxstart') || strcmp(options{ops}(1:equals-1),'minlength') || strcmp(options{ops}(1:equals-1),'maxlength') ||...
            strcmp(options{ops}(1:equals-1),'minstart') || strcmp(options{ops}(1:equals-1),'minstart'))
    else
       disp(['ERROR, the input ' options{ops} ' is not valid in module G (gap)'])
    end
    
end
    
    
    
    
% trace string   
l{1}='// Free space module of the length in the following arm';
l{2}=['COMPONENT EndOfelement_' num '= Arm()'];
l{3}=['AT (0,0,length' num ' ) RELATIVE PREVIOUS'];

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];
McStasStr.trace_seg=[McStasStr.trace_seg '\n\n' tracestring];
%McStasStr.trace=[tracestring '\n\n' McStasStr.trace ];
clear l;


if McStasStr.minimalist>-0.5
% start
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['startx' num];
    McStasStr.declare{ldec+2}=['starty' num];

    % Here a certain solution is chossen. Same var_div_req, but larger end
    % of the former guide element which enlarges the phase-space
    % requirement. But it could be done with larger var_div_req and then
    % not iluminating the entire guide start of the latter element.
    
    
    l{1} = ['startx' num ' = endx' num ' + 2*length' num '*tan(var_divreq_x*DEG2RAD);'];
    l{2} = ['starty' num ' = endy' num ' + 2*length' num '*tan(var_divreq_y*DEG2RAD);'];   
    
    % update var_divreq_x/y to current value
    % NON TRIVIAL CHANGE NEEDED!
    %l{1} = ['var_divreq_x = atan(endx' num '*tan(var_divreq_x*DEG2RAD)/startx' num ')*RAD2DEG;'];
    %l{2} = ['var_divreq_y = atan(endy' num '*tan(var_divreq_y*DEG2RAD)/starty' num ')*RAD2DEG;'];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    
else
    
if optimlimits(1)>0
    xlimdata(1)=optimlimits(1);
else
    xlimdata(1)=0.005;
end
if optimlimits(2)>0
    xlimdata(2)=optimlimits(2);
else
    xlimdata(2)=0.15;
end
xlimdata(3)=0.03;
if xlimdata(3)>xlimdata(2) || xlimdata(3)<xlimdata(1)
    xlimdata(3)=mean(xlimdata(1:2));
end    

if optimlimits(3)>0
    ylimdata(1)=optimlimits(3);
else
    ylimdata(1)=0.005;
end
if optimlimits(4)>0
    ylimdata(2)=optimlimits(4);
else
    ylimdata(2)=0.15;
end
ylimdata(3)=0.03;
if ylimdata(3)>ylimdata(2) || ylimdata(3)<ylimdata(1)
    ylimdata(3)=mean(ylimdata(1:2));
end    

% Tese next functions write the part of the instrument file for this module
% concerning start, end and length parameters.
% start
McStasStr=guide_writer_start(McStasStr,index,last,xlimdata,ylimdata,locked,globalinfo);

end

% % end
% 
%    % calculate the end so that it will illuminate the next module
%    % the size of the opening of the next module is known, and the
%    % divergencfe requirement could be kept in a single variable and
%    % overwriten at each step where nessecary.
%    
%    % what can change the divergency requirement?
%    % Everything other than a guide with start=end
%    
%    % When these modules are used, update the div_requirement accordingly
%    % but it should not be divreq_x/y, as these are needed elsewhere.
%    % let's call it var_divreq_x and var_divreq_y 
%    
%     ldec=length(McStasStr.declare);
%     McStasStr.declare{ldec+1}=['endx' num];
%     McStasStr.declare{ldec+2}=['endy' num];
%     
%     l{1}=['endx' num ' = startx' num2str(index-1) ';'];
%     l{2}=['endy' num ' = starty' num2str(index-1) ';'];
%     
%     % How do i calculate the the end of the current element?
%     
%     % case 1: locked to the former element end num = start num - 1
%     % case 2: M before, propagate in free space, M can do that!
%     
%     % This will lead to kinks and M modules have a start and end
%     
%     initstring='';
%     for i=1:length(l)
%         initstring=[initstring l{i} '\n'];
%     end
%     clear l;
%     McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    


% end
McStasStr=guide_writer_end(McStasStr,index,last);

% length 
McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,[0.01 0.25 0.02]);
end % end of overall logic, runs if no error.


end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% START OF SUPPORT FUNCTION G_MODULE_MONO() %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This is a fork of the gap module written by Leland that picks up the special exception when the gap is followed by a
% monochromator. In this situation, the Minimalist Principle must be modified. There are four different possibilities.
% Namely:
% 1)    If the monochromator parameter BeamDir = 'transmit', then the standard Minimalist Principle applies except now it
%       needs to calculate the acceptance for the entrance of the component immediately after the monochromator. Thus, the
%       equations needs to be modified from:
%       ['startx' num ' = endx' num ' + 2*length' num '*tan(var_divreq_x*DEG2RAD);'];
%       ['starty' num ' = endy' num ' + 2*length' num '*tan(var_divreq_y*DEG2RAD);'];
%       to
%       ['startx' num ' = endx' numM1 ' + 2*(length' num ' + length' numM1 ')*tan(var_divreq_x*DEG2RAD);'];
%       ['starty' num ' = endy' numM1 ' + 2*(length' num ' + length' numM1 ')*tan(var_divreq_y*DEG2RAD);'];
%
% 2)    In the BeamDir = 'transmit' case, if the user does not want the standard Minimialist Principle for a
%       gap calculated, then it can be turned off in the standard way by setting options.minimalist = -1. In this
%       case startx and starty of the gap are optimized.
%
% 3)    If the monochromator parameter BeamDir = 'reflect' then the Minimalist Principle needs to be modified to only accepted
%       neutron with a divergence that can be accepted by the mosaic of monochromator (mos). Thus, the equations need to be modified from:
%       ['startx' num ' = endx' num ' + 2*length' num '*tan(var_divreq_x*DEG2RAD);'];
%       ['starty' num ' = endy' num ' + 2*length' num '*tan(var_divreq_y*DEG2RAD);'];
%       to
%       ['startx' num ' = endx' num ' + 2*length' num '*tan(hmos_rad' numM1 ');'];
%       ['starty' num ' = endy' num ' + 2*length' num '*tan(vmos_rad' numM1 ');'];
%
% 4)    In the BeamDir = 'reflect' case, if the user does not want the Monochromator Minimialist Principle for a
%       gap calculated, then it can be turned off using the monochromator setting MinimalistMono = -1. In this
%       case startx and starty of the gap are optimized.
%
% For the most part G_module_mono and G_module_standard are exactly the same. The small modified block resides between the comments:
% 'Beginning/End of Changes to G_module'
%
function [McStasStr] = G_module_mono(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
%M_module adds text to McStas bot to make space between two guide elements
   % add the relevant component to the trace section
   % add the relevant input to the input list
   % add the relevant calculations to the initialize section
   % add the relevant variable names to the declare list
   
   
   
num=num2str(index);
numM1=num2str(index-1);

if (index==last || index ==1) 
    if (index==last)
        disp('ERROR, do not place G at the start of input-string, change closest_element instead')
        disp('ERROR, the first G module will be ignored')
    end
    if (index==1)
        disp('ERROR, do not place G at the end of input-string, change sample_dist instead')
        disp('ERROR, the last G module will be ignored')
    end
else % add a free space to the instrument file

defaultnames=fieldnames(defaults);
% Input not concerning the geometry
% WARNING, the indexes are used in the part extracting data from options
% linp=length(McStasStr.input);
% R0Index=linp+1;
% McStasStr.input{R0Index}=['R0' num];
% if ismember(defaultnames,'R0')
% McStasStr.inputvalue(R0Index)=defaults.R0;
% else    
% McStasStr.inputvalue(R0Index)=0.99;
% end
% 
% QcIndex=linp+2;
% McStasStr.input{QcIndex}=['Qc' num];
% if ismember(defaultnames,'Qc')
% McStasStr.inputvalue(QcIndex)=defaults.Qc;    
% else
% McStasStr.inputvalue(QcIndex)=0.0217;
% end
% 
% alphaIndex=linp+3;
% McStasStr.input{alphaIndex}=['alpha' num];
% if ismember(defaultnames,'alpha')
% McStasStr.inputvalue(alphaIndex)=defaults.alpha;
% else
% McStasStr.inputvalue(alphaIndex)=6.07;
% end
% 
% mIndex=linp+4;
% McStasStr.input{mIndex}=['m' num];
% if ismember(defaultnames,'m')
% McStasStr.inputvalue(mIndex)=defaults.m;
% else
% McStasStr.inputvalue(mIndex)=3;
% end
% 
% WIndex=linp+5;
% McStasStr.input{WIndex}=['W' num];
% if ismember(defaultnames,'W')
% McStasStr.inputvalue(WIndex)=defaults.W;
% else
% McStasStr.inputvalue(WIndex)=0.003;
% end

optimlimits=-ones(1,4);
locked=-ones(1,2); % width,height
% Extracting options from the option cell
for ops=1:length(options)
    for i=1:length(options{ops})
        if strcmp(options{ops}(i),'='); 
            equals=i;
        end
    end
    
    % check for options on reflectivity curve
    if strcmp(options{ops}(1:equals-1),'minStartWidth')
       if optimlimits(1)>0
           optimlimits(1)=max([optimlimits(1) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(1)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'maxStartWidth')
       if optimlimits(2)>0
           optimilimts(2)=min([optimlimits(2) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(2)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'minStartHeight')
       if optimlimits(3)>0
           optimilimts(3)=max([optimlimits(3) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(3)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'maxStartHeight')
       if optimlimits(4)>0
           optimilimts(4)=min([optimlimits(4) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(4)=str2num(options{ops}(equals+1:end));   
       end
    elseif strcmp(options{ops}(1:equals-1),'StartWidth')
       if locked(1)>0
           locked(1)=min([locked(1) str2num(options{ops}(equals+1:end))]);
       else
           locked(1)=str2num(options{ops}(equals+1:end));   
       end
    elseif strcmp(options{ops}(1:equals-1),'StartHeight')
       if locked(2)>0
           locked(2)=min([locked(2) str2num(options{ops}(equals+1:end))]);
       else
           locked(2)=str2num(options{ops}(equals+1:end));   
       end
       % These avoid errors assosicated with options which are treated by
       % the main loop. changed is for when options are removed. 
    elseif (strcmp(options{ops}(1:equals-1),'start') || strcmp(options{ops}(1:equals-1),'length') || strcmp(options{ops}(1:equals-1),'changed') ||...
            strcmp(options{ops}(1:equals-1),'maxstart') || strcmp(options{ops}(1:equals-1),'minlength') || strcmp(options{ops}(1:equals-1),'maxlength') ||...
            strcmp(options{ops}(1:equals-1),'minstart') || strcmp(options{ops}(1:equals-1),'minstart'))
    else
       disp(['ERROR, the input ' options{ops} ' is not valid in module G (gap)'])
    end
    
end
    
    
    
    
% trace string   
l{1}='// Free space module of the length in the following arm';
l{2}=['COMPONENT EndOfelement_' num '= Arm()'];
l{3}=['AT (0,0,length' num ' ) RELATIVE PREVIOUS'];

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];
McStasStr.trace_seg=[McStasStr.trace_seg '\n\n' tracestring];
%McStasStr.trace=[tracestring '\n\n' McStasStr.trace ];
clear l;

%%%% Beginning of Changes to G_module %%%%

Mono_options = Parse_options(globalinfo.options{last-index+2});
% start
ldec=length(McStasStr.declare);
McStasStr.declare{ldec+1}=['startx' num];
McStasStr.declare{ldec+2}=['starty' num];


if (McStasStr.minimalist>-0.5  &&  strcmp(Mono_options.BeamDir, 'transmit'))


    % Here a certain solution is chossen. Same var_div_req, but larger end
    % of the former guide element which enlarges the phase-space
    % requirement. But it could be done with larger var_div_req and then
    % not iluminating the entire guide start of the latter element.
    
    l{1} = ['startx' num ' = endx' numM1 ' + 2*(length' num ' + length' numM1 ')*tan(var_divreq_x*DEG2RAD);'];
    l{2} = ['starty' num ' = endy' numM1 ' + 2*(length' num ' + length' numM1 ')*tan(var_divreq_y*DEG2RAD);'];
    
    % update var_divreq_x/y to current value
    % NON TRIVIAL CHANGE NEEDED!
    %l{1} = ['var_divreq_x = atan(endx' num '*tan(var_divreq_x*DEG2RAD)/startx' num ')*RAD2DEG;'];
    %l{2} = ['var_divreq_y = atan(endy' num '*tan(var_divreq_y*DEG2RAD)/starty' num ')*RAD2DEG;'];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
    
elseif (Mono_options.MinimalistMono == 1  &&  strcmp(Mono_options.BeamDir, 'reflect'))
    l{1} = ['startx' num ' = endx' num ' + 2*length' num '*tan(hmos_rad' numM1 ');'];
    l{2} = ['starty' num ' = endy' num ' + 2*length' num '*tan(vmos_rad' numM1 ');'];
        initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
    
else
    
    
if optimlimits(1)>0
    xlimdata(1)=optimlimits(1);
else
    xlimdata(1)=0.005;
end
if optimlimits(2)>0
    xlimdata(2)=optimlimits(2);
else
    xlimdata(2)=Mono_options.length;
end
xlimdata(3)= mean(xlimdata(1:2));
    

if optimlimits(3)>0
    ylimdata(1)=optimlimits(3);
else
    ylimdata(1)=0.005;
end
if optimlimits(4)>0
    ylimdata(2)=optimlimits(4);
else
    ylimdata(2)=Mono_options.BladeH;
end
ylimdata(3)=mean(ylimdata(1:2));


%%%% End of Changes to G_module %%%%%



% Tese next functions write the part of the instrument file for this module
% concerning start, end and length parameters.
% start
McStasStr=guide_writer_start(McStasStr,index,last,xlimdata,ylimdata,locked,globalinfo);

end

% % end
% 
%    % calculate the end so that it will illuminate the next module
%    % the size of the opening of the next module is known, and the
%    % divergencfe requirement could be kept in a single variable and
%    % overwriten at each step where nessecary.
%    
%    % what can change the divergency requirement?
%    % Everything other than a guide with start=end
%    
%    % When these modules are used, update the div_requirement accordingly
%    % but it should not be divreq_x/y, as these are needed elsewhere.
%    % let's call it var_divreq_x and var_divreq_y 
%    
%     ldec=length(McStasStr.declare);
%     McStasStr.declare{ldec+1}=['endx' num];
%     McStasStr.declare{ldec+2}=['endy' num];
%     
%     l{1}=['endx' num ' = startx' num2str(index-1) ';'];
%     l{2}=['endy' num ' = starty' num2str(index-1) ';'];
%     
%     % How do i calculate the the end of the current element?
%     
%     % case 1: locked to the former element end num = start num - 1
%     % case 2: M before, propagate in free space, M can do that!
%     
%     % This will lead to kinks and M modules have a start and end
%     
%     initstring='';
%     for i=1:length(l)
%         initstring=[initstring l{i} '\n'];
%     end
%     clear l;
%     McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    


% end
McStasStr=guide_writer_end(McStasStr,index,last);

% length 
McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,[0.01 0.25 0.02]);
end % end of overall logic, runs if no error.


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% END OF SUPPORT FUNCTION G_MODULE_MON() %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% START OF SUPPORT FUNCTION PARSE_OPTIONS() %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StructOut = Parse_options(options)
for ops=1:length(options)
    command = options{ops};
    try % Assume command is numeric
        eval(['StructArg.' command ';']);
    catch % If not numeric then store command as a string;
        for n = 1:length(command)
            if strcmp(command(n),'=');
                equals=n;
            end
        end
        try
            eval(['StructArg.' command(1:equals-1) '=''' command(equals+1:end) ''';']);
        end
    end
end
StructOut = StructArg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% END OF SUPPORT FUNCTION PARSE_OPTIONS() %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

