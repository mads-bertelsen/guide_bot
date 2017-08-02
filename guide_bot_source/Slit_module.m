function [McStasStr] = Slit_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
%Smodule adds text to McStas bot to make a straight guide
   % Adds the component to trace 
   % Handles the start parameters
   % Handles the end parameters
   % Handles the length parameter
   % Handles the endPoint parameter
   % Handles the options given to the function
   
num=num2str(index);
numM1=num2str(index-1);

defaultnames=fieldnames(defaults);
% linp=length(McStasStr.input);
% R0Index=linp+1;
% McStasStr.input{R0Index}=['R0' num];
% if ismember(defaultnames,'R0')
% McStasStr.inputvalue(R0Index)=defaults.R0;
% else    
% McStasStr.inputvalue(R0Index)=0.99;
% end
horizontal_displacement = 0;
vertical_displacement = 0;


optimlimits=-ones(1,4);
locked=-ones(1,2); % width,height
% Extracting options from the option cell
for ops=1:length(options)
    for i=1:length(options{ops})
        if strcmp(options{ops}(i),'='); 
            equals=i;
        end
    end
    
    % This component works with some exceptions in the option handling in
    % the main mcstas_bot file.
    % Basicly, StartWidth options are passed to the following module,
    % which will effectivly make the slit have that option because it has
    % equal start and end dimensions.
    % There is a special case, if there is no module after which can accept
    % the StartWidth like options.
    % I have yet to decide what to do in this case, so far it does not
    % work.
    
    if strcmp(options{ops}(1:equals-1),'R0')
       McStasStr.inputvalue(R0Index)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'minStartWidth')
%        if optimlimits(1)>0
%            optimlimits(1)=max([optimlimits(1) str2num(options{ops}(equals+1:end))]);
%        else
%            optimlimits(1)=str2num(options{ops}(equals+1:end));
%        end
    disp('The minStartWidth options is ignored on slit when it is the last module;')
    elseif strcmp(options{ops}(1:equals-1),'maxStartWidth')
%        if optimlimits(2)>0
%            optimilimts(2)=min([optimlimits(2) str2num(options{ops}(equals+1:end))]);
%        else
%            optimlimits(2)=str2num(options{ops}(equals+1:end));
%        end
    disp('The maxStartWidth options is ignored on slit when it is the last module;')
    elseif strcmp(options{ops}(1:equals-1),'minStartHeight')
%        if optimlimits(3)>0
%            optimilimts(3)=max([optimlimits(3) str2num(options{ops}(equals+1:end))]);
%        else
%            optimlimits(3)=str2num(options{ops}(equals+1:end));
%        end
    disp('The minStartHeight options is ignored on slit when it is the last module;')
    elseif strcmp(options{ops}(1:equals-1),'maxStartHeight')
%        if optimlimits(4)>0
%            optimilimts(4)=min([optimlimits(4) str2num(options{ops}(equals+1:end))]);
%        else
%            optimlimits(4)=str2num(options{ops}(equals+1:end));   
%        end
    disp('The maxStartHeight options is ignored on slit when it is the last module;')
    elseif strcmp(options{ops}(1:equals-1),'StartWidth')
%        if locked(1)>0
%            locked(1)=min([locked(1) str2num(options{ops}(equals+1:end))]);
%        else
%            locked(1)=str2num(options{ops}(equals+1:end));   
%        end
    disp('The StartWidth options is ignored on slit when it is the last module;')
    elseif strcmp(options{ops}(1:equals-1),'StartHeight')
%        if locked(2)>0
%            locked(2)=min([locked(2) str2num(options{ops}(equals+1:end))]);
%        else
%            locked(2)=str2num(options{ops}(equals+1:end));   
%        end
    disp('The StartHeight options is ignored on slit when it is the last module;')
    elseif strcmp(options{ops}(1:equals-1),'H_displacement')
        horizontal_displacement = str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'V_displacement')
        vertical_displacement = str2num(options{ops}(equals+1:end));
    elseif (strcmp(options{ops}(1:equals-1),'start') || strcmp(options{ops}(1:equals-1),'minstart') || strcmp(options{ops}(1:equals-1),'changed') ||...
            strcmp(options{ops}(1:equals-1),'maxstart'))
    elseif (strcmp(options{ops}(1:equals-1),'length') || strcmp(options{ops}(1:equals-1),'maxlength') || strcmp(options{ops}(1:equals-1),'minlength'))
        disp('ERROR,The Slit module have no length, use a gap (G) after or before instead. Option ignored.')
    else
       disp(['ERROR, the input ' options{ops} ' is not valid in module Slit'])
    end
    
end

% trace string   
l{1}=['COMPONENT slit_' num '= Slit('];
l{2}=['xwidth=startx' num ', yheight=starty' num ')'];
l{3}=['AT (' num2str(horizontal_displacement) ',' num2str(vertical_displacement) ',u) RELATIVE PREVIOUS'];    
l{4}='ROTATED (0,0,0) RELATIVE PREVIOUS';
l{5}='';
l{6}=['COMPONENT EndOfelement_' num '= Arm()'];
l{7}=['AT (-' num2str(horizontal_displacement) ',-' num2str(vertical_displacement) ',u) RELATIVE PREVIOUS'];    

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];
McStasStr.trace_seg=[McStasStr.trace_seg '\n\n' tracestring];
clear l;

% trace string

% Limits
if optimlimits(1)>0
    xlimdata(1)=optimlimits(1);
else
    xlimdata(1)=0.005;
end
if optimlimits(2)>0
    xlimdata(2)=optimlimits(2);
else
    xlimdata(2)=0.2;
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
    ylimdata(2)=0.2;
end
ylimdata(3)=0.03;
if ylimdata(3)>ylimdata(2) || ylimdata(3)<ylimdata(1)
    ylimdata(3)=mean(ylimdata(1:2));
end    

% Tese next functions write the part of the instrument file for this module
% concerning start, end and length parameters.
% start
%McStasStr=guide_writer_start(McStasStr,index,last,xlimdata,ylimdata,locked,globalinfo);
% The last two vectors is the [min,max,guess] of (x,y) size

% end
%McStasStr=guide_writer_end(McStasStr,index,last);

% length
% This module does not have a length!
% McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,[0.1 0.95]);

ldec=length(McStasStr.declare);
McStasStr.declare{ldec+1}=['endPoint' num];
McStasStr.declare{ldec+2}=['length' num];
l{1} = ['length' num ' = 2*u;'];
if (index==1)
    l{end+1}=['endPoint' num ' = Mod_sample - sample_dist;'];
else
    l{end+1}=['endPoint' num ' = endPoint' numM1 ' - length' numM1 ' - u;'];
end

initstring='';
for i=1:length(l)
    initstring=[initstring l{i} '\n'];
end
clear l;
McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    


if (index==last) && McStasStr.minimalist==1 %start can be calculated
    
    % declare start and end
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['startx' num];
    McStasStr.declare{ldec+2}=['starty' num];
    McStasStr.declare{ldec+3}=['endx' num];
    McStasStr.declare{ldec+4}=['endy' num];
    
    % var_divreq known
    % calculate dimensions of the guide opening so that it is illuminated
    % in phase-space
    
    l{1}=['if (0<mod_x*mod_x - endx' num '*tan(var_divreq_x*DEG2RAD)*8*guide_start) {'];
    l{1+end}=['startx' num ' = mod_x/2 - sqrt(mod_x*mod_x - endx' num '*tan(var_divreq_x*DEG2RAD)*8*guide_start)/2;'];
    l{1+end}=['}'];
    l{1+end}=['else {'];
    l{1+end}=['startx' num ' = mod_x/2;'];
    l{1+end}=['}'];
    l{1+end}=['if (0<mod_y*mod_y - endy' num '*tan(var_divreq_y*DEG2RAD)*8*guide_start) {'];
    l{1+end}=['starty' num ' = mod_y/2 - sqrt(mod_y*mod_y - endy' num '*tan(var_divreq_y*DEG2RAD)*8*guide_start)/2;']; 
    l{1+end}=['}'];
    l{1+end}=['else {'];
    l{1+end}=['starty' num ' = mod_y/2;'];
    l{1+end}=['}'];
    l{1+end}=['endx' num ' = startx' num ';'];
    l{1+end}=['endy' num ' = starty' num ';'];
    
    % Add this calculation to the mcstas initialize section
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
       

elseif (index==1) % first module after the sample
    
    % declare start  
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['startx' num];
    McStasStr.declare{ldec+2}=['starty' num];
    McStasStr.declare{ldec+3}=['endx' num];
    McStasStr.declare{ldec+4}=['endy' num];
    
    % Calculate end (towards sample)
    l{1}=['endx' num ' = sizeX + 2*sample_dist*tan(divreq_x*DEG2RAD);'];
    l{2}=['endy' num ' = sizeY + 2*sample_dist*tan(divreq_y*DEG2RAD);'];
    l{3}=['startx' num ' = endx' num ';'];
    l{4}=['starty' num ' = endy' num ';'];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];   
    
else
    
    % declare start and end
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['startx' num];
    McStasStr.declare{ldec+2}=['starty' num];
    McStasStr.declare{ldec+3}=['endx' num];
    McStasStr.declare{ldec+4}=['endy' num];
    
    
    % Connect the slit to previous module
    l{1}=['endx' num ' = startx' num2str(index-1) ';'];
    l{2}=['endy' num ' = starty' num2str(index-1) ';'];
    l{3}=['startx' num ' = endx' num ';'];
    l{4}=['starty' num ' = endy' num ';'];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];   
    
end


end

