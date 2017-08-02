function [McStasStr] = C_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
%Smodule adds text to McStas bot to make a straight guide
   % Adds the component to trace 
   % Handles the start parameters
   % Handles the end parameters
   % Handles the length parameter
   % Handles the endPoint parameter
   % Handles the options given to the function

if (index==1 && index==last)
    disp('ERROR, C (curved guide) can not be the only component!')
elseif (index==1)
    disp('ERROR, C (curved guide) can not be the last component!')
elseif (index==last)
    disp('ERROR, C (curved guide) can not be the first component!')
else
% Run the module as all assumptions are met.

% default options   
for ops=1:length(options)    
    keyword='rot=';
    klength=length(keyword);
    % check for possible option rot
    if length(options{ops}) > klength
        if strcmp(options{ops}(1:klength),keyword)
            rot=str2num(options{ops}(klength+1:end));
            rot=abs(rot);
            overwrite=1;
        end
    end
end
    
num=num2str(index);
numM1=num2str(index-1);

% initialize
rotlogic = 0; % default is calculating the angle via RT
% Default is a curved guide, but can be overwritten by options
channels = 1;
mirror_thickness = 0.0001;

defaultnames=fieldnames(defaults);
% Input not concerning the geometry
% WARNING, the indexes are used in the part extracting data from options
linp=length(McStasStr.input);
R0Index=linp+1;
McStasStr.input{R0Index}=['R0' num];
if max(ismember(defaultnames,'R0'))
McStasStr.inputvalue(R0Index)=defaults.R0;
else    
McStasStr.inputvalue(R0Index)=0.99;
end

QcIndex=linp+2;
McStasStr.input{QcIndex}=['Qc' num];
if max(ismember(defaultnames,'Qc'))
McStasStr.inputvalue(QcIndex)=defaults.Qc;    
else
McStasStr.inputvalue(QcIndex)=0.0217;
end

alphaIndex=linp+3;
McStasStr.input{alphaIndex}=['alpha' num];
if max(ismember(defaultnames,'alpha'))
McStasStr.inputvalue(alphaIndex)=defaults.alpha;
else
McStasStr.inputvalue(alphaIndex)=6.07;
end

mIndex=linp+4;
McStasStr.input{mIndex}=['m' num];
if max(ismember(defaultnames,'m'))
McStasStr.inputvalue(mIndex)=defaults.m;
else
McStasStr.inputvalue(mIndex)=3;
end

WIndex=linp+5;
McStasStr.input{WIndex}=['W' num];
if max(ismember(defaultnames,'W'))
McStasStr.inputvalue(WIndex)=defaults.W;
else
McStasStr.inputvalue(WIndex)=0.003;
end

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
    if strcmp(options{ops}(1:equals-1),'R0')
       McStasStr.inputvalue(R0Index)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'Qc')
       McStasStr.inputvalue(QcIndex)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'alpha')
       McStasStr.inputvalue(alphaIndex)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'m')
       McStasStr.inputvalue(mIndex)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'W')
       McStasStr.inputvalue(WIndex)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'rot')
       rotlogic = 1;
       rot = str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'bender')
       channels = str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'mirror_thickness')
       mirror_thickness = str2num(options{ops}(equals+1:end));   
    elseif strcmp(options{ops}(1:equals-1),'minStartWidth')
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
    elseif (strcmp(options{ops}(1:equals-1),'start') || strcmp(options{ops}(1:equals-1),'length') || strcmp(options{ops}(1:equals-1),'changed') ||...
            strcmp(options{ops}(1:equals-1),'maxstart') || strcmp(options{ops}(1:equals-1),'minstart') || strcmp(options{ops}(1:equals-1),'maxlength') ||...
            strcmp(options{ops}(1:equals-1),'minlength') || strcmp(options{ops}(1:equals-1),'rots') || strcmp(options{ops}(1:equals-1),'rot') ||...
            strcmp(options{ops}(1:equals-1),'los_start') || strcmp(options{ops}(1:equals-1),'los_end') || strcmp(options{ops}(1:equals-1),'los_divide') ||...
            strcmp(options{ops}(1:equals-1),'rotd'))
    else
       disp(['ERROR, the input ' options{ops} ' is not valid in module C (Curved guide)'])
    end
    
end

% trace string
% globalinfo.rotsign should be used to select the direction of the curve.
if globalinfo.rotlogic(index)>0
l{1}=['COMPONENT Curved_guide_' num '= Bender('];
l{2}=['w=startx' num ', h=starty' num ', l=length' num ','];
l{3}=['r=curve_radius' num ', k = channels' num ', d = ' num2str(mirror_thickness) ', ' ];
l{4}=['R0a=R0' num ', Qca=Qc' num ', alphaa=alpha' num ', ma=m' num ', Wa=W' num ','];
l{5}=['R0i=R0' num ', Qci=Qc' num ', alphai=alpha' num ', mi=m' num ', Wi=W' num ','];
l{6}=['R0s=R0' num ', Qcs=Qc' num ', alphas=alpha' num ', ms=m' num ', Ws=W' num ')'];
l{7}='AT (0,0, u) RELATIVE PREVIOUS';    
if globalinfo.rotsign(index)>0    
l{8}='ROTATED (0,0,0) RELATIVE PREVIOUS'; % Bend -x direction
else
l{8}='ROTATED (0,0,180) RELATIVE PREVIOUS'; % Bend x direction
end
l{9}='';
l{10}=['COMPONENT EndOfelement_' num '= Arm()'];
l{11}=['AT (0,0,length' num ' + 2*u) RELATIVE PREVIOUS'];
if globalinfo.rotsign(index)>0    
l{end+1}='ROTATED (0,0,0) RELATIVE PREVIOUS'; % Bend -x direction
else
l{end+1}='ROTATED (0,0,-180) RELATIVE PREVIOUS'; % Bend x direction
end
else
% This string should have a curved guide which curves in the vertical plane 
% It will have to be checked if this actually works
l{1}=['COMPONENT Curved_guide_' num '= Bender('];
l{2}=['w=starty' num ', h=startx' num ', l=length' num ',']; % since rotated, replace x and y.
l{3}=['r=curve_radius' num ', k = channels' num ', d = ' num2str(mirror_thickness) ', ' ];
l{4}=['R0a=R0' num ', Qca=Qc' num ', alphaa=alpha' num ', ma=m' num ', Wa=W' num ','];
l{5}=['R0i=R0' num ', Qci=Qc' num ', alphai=alpha' num ', mi=m' num ', Wi=W' num ','];
l{6}=['R0s=R0' num ', Qcs=Qc' num ', alphas=alpha' num ', ms=m' num ', Ws=W' num ')'];
l{7}='AT (0,0, u) RELATIVE PREVIOUS';    
if globalinfo.rotsign(index)>0    
l{8}='ROTATED (0,0,90) RELATIVE PREVIOUS'; % Bend -y direction
else
l{8}='ROTATED (0,0,-90) RELATIVE PREVIOUS'; % Bend y direction
end
l{9}='';
l{10}=['COMPONENT EndOfelement_' num '= Arm()'];
l{11}=['AT (0,0,length' num ' + 2*u) RELATIVE PREVIOUS'];
if globalinfo.rotsign(index)>0    
l{end+1}='ROTATED (0,0,-90) RELATIVE PREVIOUS'; % Bend -y direction
else
l{end+1}='ROTATED (0,0,90) RELATIVE PREVIOUS'; % Bend y direction
end
end
    
% The McStas compontn Bender does not actually bend the guide, but the
% transmission is as if it was bent which is done by changing the neutrons
% velocity instead of curving the guide. Because of this, the EndOfelement
% arm does not need to be rotated or translated.

% curved guide pseudo code:

% This crappy component does not support different opening and exit
% dimensions, and thus can not use the standard functions / logic.

% if index = 1 The end is calculated from requirements
% if index = last The start is calculated from requirements
% if index = last = 1 ERROR impossible with this crap 
% if index + 1 and index -1 have different dimensions, what then? ERROR
% if index have max min start, it should be passed on to index -1 start

% need to calculate curve_radius from rot
% 1: could do it with the component
% 2: could do it my self because the curve radius is very important

% r = 2*pi*L / (rot * DEG2RAD)

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];
clear l;

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


% The last two vectors is the [min,max,guess] of (x,y) size
% renamed due to function renamning
startxpar=xlimdata;
startypar=ylimdata;

% Custom written start and end writer, but standard length.
% The width is never optimized, but taken from start num - 1.
% If the index is first or last, the width is calculated.

% Decided that it is not needed to start or end with a curved guide.
% Thus the index == 1 and index == last are not nessecary.


% What to do if index == last and minimalist == 0 ? 
% if (index==last) && McStasStr.minimalist==1 %start can be calculated
%     
%     % declare start  
%     ldec=length(McStasStr.declare);
%     McStasStr.declare{ldec+1}=['startx' num];
%     McStasStr.declare{ldec+2}=['starty' num];
%     
%     % var_divreq known
%     % calculate dimensions of the guide opening so that it is illuminated
%     % in phase-space
%     
%     l{1}=['if (0<mod*mod - endx' num '*tan(var_divreq_x*DEG2RAD)*8*guide_start) {'];
%     l{1+end}=['startx' num ' = mod/2 - sqrt(mod*mod - endx' num '*tan(var_divreq_x*DEG2RAD)*8*guide_start)/2;'];
%     l{1+end}=['}'];
%     l{1+end}=['else {'];
%     l{1+end}=['startx' num ' = mod/2;'];
%     l{1+end}=['}'];
%     l{1+end}=['if (0<mod*mod - endy' num '*tan(var_divreq_y*DEG2RAD)*8*guide_start) {'];
%     l{1+end}=['starty' num ' = mod/2 - sqrt(mod*mod - endy' num '*tan(var_divreq_y*DEG2RAD)*8*guide_start)/2;']; 
%     l{1+end}=['}'];
%     l{1+end}=['else {'];
%     l{1+end}=['starty' num ' = mod/2;'];
%     l{1+end}=['}'];
%     
%     % Add this calculation to the mcstas initialize section
%     initstring='';
%     for i=1:length(l)
%         initstring=[initstring l{i} '\n'];
%     end
%     clear l;
%     McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    
if (index == last)
    disp('ERROR, module C (curved guide) can not start a guide');
elseif (index == 1)
    disp('ERROR, module C (curved guide) can not end a guide');
else
    % End should allready be set to the start of the guide element after
    % the curved guide.
    
    % can be calculated
    % declare start 
    % assign start the end of the guide
    
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['startx' num];
    McStasStr.declare{ldec+2}=['starty' num];
    
    l{1}=['startx' num ' = endx' num ';'];
    l{2}=['starty' num ' = endy' num ';'];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
end

if (index==1) % first module after the sample
    disp('ERROR, module C (curved guide) can not end a guide')
%     ldec=length(McStasStr.declare);
%     McStasStr.declare{ldec+1}=['endx' num];
%     McStasStr.declare{ldec+2}=['endy' num];
%     
%     l{1}=['endx' num ' = sizeX + 2*sample_dist*tan(divreq_x*DEG2RAD);'];
%     l{2}=['endy' num ' = sizeY + 2*sample_dist*tan(divreq_y*DEG2RAD);'];
%     
%     initstring='';
%     for i=1:length(l)
%         initstring=[initstring l{i} '\n'];
%     end
%     clear l;
%     McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
elseif (index == last)
    disp('ERROR, module C (curved guide) can not start a guide')
%         ldec=length(McStasStr.declare);
%         McStasStr.declare{ldec+1}=['endx' num];
%         McStasStr.declare{ldec+2}=['endy' num];
% 
%         l{1}=['endx' num ' = startx' num ';'];
%         l{2}=['endy' num ' = starty' num ';'];
% 
%         initstring='';
%         for i=1:length(l)
%             initstring=[initstring l{i} '\n'];
%         end
%         clear l;
%         McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    
else % index not equal to 1 or last.
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['endx' num];
    McStasStr.declare{ldec+2}=['endy' num];

    l{1}=['endx' num ' = startx' num2str(index-1) ';'];
    l{2}=['endy' num ' = starty' num2str(index-1) ';'];

    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    
end

linp=length(McStasStr.input);
McStasStr.input{linp+1}=['channels' num];
McStasStr.inputvalue(linp+1)=channels;

% Calculate the curve radius
ldec=length(McStasStr.declare);
McStasStr.declare{ldec+1}=['curve_radius' num];

if globalinfo.curveoverwrite(index) == 1
    % overwrite and put as input
    linp=length(McStasStr.input);
    McStasStr.input{linp+1}=['rot' num];
    McStasStr.inputvalue(linp+1)=rot;
    % r = 2*pi*L / (rot * DEG2RAD)
    l{1}=['curve_radius' num ' = length' num '/(rot' num '*DEG2RAD);'];
else
    current_los_logic = false(1,last);
    % experimental code
    relevant_ghost_list = ghost_globalinfo.active_los{ghost_globalinfo.los_section};
    if globalinfo.curve_logic(ghost_globalinfo.real_index(relevant_ghost_list(end)))
        %disp('first if entered         --- // ---')
        % check if it should be removed
        % Should be removed if only the end point is within current active los
        same_real_index_logic = ghost_globalinfo.real_index == ghost_globalinfo.real_index(relevant_ghost_list(end));
        same_real_index = find(same_real_index_logic);
        if same_real_index(end) == relevant_ghost_list(end)
            % remove from list, only the end point is used!
            relevant_ghost_list = relevant_ghost_list(1:end-1);
        end
    end
    % end experimental code. WORKS! Keep it for now.
    
    % current_los_logic(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section})) = 1;
    current_los_logic(ghost_globalinfo.real_index(relevant_ghost_list)) = 1;
    current_raytracers_logic = current_los_logic.*globalinfo.raytracers;
    current_raytracers_index = find(current_raytracers_logic);
    master_rot_index = current_raytracers_index(1);
    
    if index == master_rot_index
    % last call in a active raytracer section
    % calculate with increment method as declare
    % Need to start the while loop here
        McStasStr.declare{end+1}=['rot' num];
    if sum(ismember(McStasStr.declare,'var_divreq_x_protected'))<0.5 % not always needed
        McStasStr.declare{end+1}='var_divreq_x_protected';
        McStasStr.declare{end+1}='var_divreq_y_protected';
    end
   % save needed variables to protect
   % start while
   % do angle increment
    l{1}='';
    l{end+1}=['var_divreq_x_protected = var_divreq_x;'];
    l{end+1}=['var_divreq_y_protected = var_divreq_y;'];
    l{end+1}=[''];
    l{end+1}=['rot' num ' = 0;'];
    l{end+1}=['los_logic = 1;'];
    l{end+1}='while(los_logic==1) {';
    l{end+1}=['rot' num ' = rot' num ' + 0.0002;'];
    l{end+1}=['curve_radius' num ' = length' num '/(rot' num '*DEG2RAD);'];
    l{end+1}=[''];
    l{end+1}=['var_divreq_x = var_divreq_x_protected;'];
    l{end+1}=['var_divreq_y = var_divreq_y_protected;'];
    l{end+1}=[''];
    
    else % if not last los_blocker in raytracer section.
        
            linp=length(McStasStr.input);
            McStasStr.input{linp+1}=['rot_ratio' num];
            McStasStr.optimize(linp+1)=1;
            % An important application for gradient renormalization
            % These optimal ranges needs to be considered
            McStasStr.optimvals.min(linp+1)=0.2;
            McStasStr.optimvals.max(linp+1)=5;
            McStasStr.optimvals.guess(linp+1)=1;
            
            McStasStr.declare{end+1}=['rot' num];
            
            l{1}='';
            l{end+1}=['rot' num ' = rot' num2str(master_rot_index) '*rot_ratio' num ';'];
            l{end+1}=['curve_radius' num ' = length' num '/(rot' num '*DEG2RAD);'];
            l{end+1}=[''];
    end
end

initstring='';
for i=1:length(l)
    initstring=[initstring l{i} '\n'];
end
clear l;
McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    

% length
McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,[0.1 0.95]);


   
% 
% COMPONENT Bender = Bender(
%     w = startx1, h = starty1, r = 400, k = 1, d = 0.001, l = 50,
%     R0a = R0a, Qca = qca, alphaa = alphaa, ma = ma, Wa = Wa,
%     R0i = R0i, Qci = Qci, alphai = alphai, mi = mi, Wi = Wi,
%     R0s = R0s, Qcs = Qcs, alphas = alphas, ms = ms, Ws = Ws)
%   AT (0, 0, u) RELATIVE Previous
%   ROTATED (0, 0, 0) RELATIVE Previous

end % end if checking if C is the only module
end