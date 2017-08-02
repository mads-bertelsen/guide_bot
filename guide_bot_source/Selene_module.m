function [McStasStr] = Selene_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
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
% Input not concerning the geometry
% WARNING, the indexes are used in the part extracting data from options
linp=length(McStasStr.input);
R0Index=linp+1;
McStasStr.input{R0Index}=['R0' num];
if ismember(defaultnames,'R0')
McStasStr.inputvalue(R0Index)=defaults.R0;
else    
McStasStr.inputvalue(R0Index)=0.99;
end

QcIndex=linp+2;
McStasStr.input{QcIndex}=['Qc' num];
if ismember(defaultnames,'Qc')
McStasStr.inputvalue(QcIndex)=defaults.Qc;    
else
McStasStr.inputvalue(QcIndex)=0.0217;
end

alphaIndex=linp+3;
McStasStr.input{alphaIndex}=['alpha' num];
if ismember(defaultnames,'alpha')
McStasStr.inputvalue(alphaIndex)=defaults.alpha;
else
McStasStr.inputvalue(alphaIndex)=6.07;
end

mIndex=linp+4;
McStasStr.input{mIndex}=['m' num];
if ismember(defaultnames,'m')
McStasStr.inputvalue(mIndex)=defaults.m;
else
McStasStr.inputvalue(mIndex)=3;
end

WIndex=linp+5;
McStasStr.input{WIndex}=['W' num];
if ismember(defaultnames,'W')
McStasStr.inputvalue(WIndex)=defaults.W;
else
McStasStr.inputvalue(WIndex)=0.003;
end

max_m = -1;
force_slit = 0;
gravity = 1;
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
    elseif strcmp(options{ops}(1:equals-1),'max_m')
       max_m = str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'force_slit')
       if str2num(options{ops}(equals+1:end)) == 1
           force_slit = 1;
       else
           force_slit = 0;
       end
    elseif strcmp(options{ops}(1:equals-1),'gravity')
       if str2num(options{ops}(equals+1:end)) == 0
           gravity = 0;
       else
           gravity = 1;
       end
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
            strcmp(options{ops}(1:equals-1),'minlength'))
    else
       disp(['ERROR, the input ' options{ops} ' is not valid in module E (Elliptic guide)'])
    end
    
end

l{1}=[''];

if (index == last && globalinfo.selene_moderator_focus == 0 || force_slit == 1)
l{end+1}=['COMPONENT slit_collimation1 = Slit('];
l{end+1}=['    xmin = -startx' num '*0.5, xmax = startx' num '*0.5,'];
l{end+1}=['    ymin = -starty' num '*0.5, ymax = starty' num '*0.5)'];
l{end+1}=['  AT (0, 0, 0) RELATIVE PREVIOUS'];
l{end+1}=[''];
elseif index == last && globalinfo.selene_moderator_focus == 1
    if globalinfo.selene_moderator_focus_v == 1 && globalinfo.selene_moderator_focus_h == 1
    % Do nothing!
    elseif globalinfo.selene_moderator_focus_v == 1
    l{end+1}=['COMPONENT slit_collimation1 = Slit('];
    l{end+1}=['    xmin = -startx' num '*0.5, xmax = startx' num '*0.5,'];
    l{end+1}=['    ymin = -10, ymax = 10)'];
    l{end+1}=['  AT (0, 0, 0) RELATIVE PREVIOUS'];
    l{end+1}=[''];
    elseif globalinfo.selene_moderator_focus_h == 1
    l{end+1}=['COMPONENT slit_collimation1 = Slit('];
    l{end+1}=['    xmin = 10, xmax = 10,'];
    l{end+1}=['    ymin = -starty' num '*0.5, ymax = starty' num '*0.5)'];
    l{end+1}=['  AT (0, 0, 0) RELATIVE PREVIOUS'];
    l{end+1}=[''];    
    end
end

l{end+1} = ['COMPONENT arm_selene_' num ' = Arm()'];
l{end+1}=['AT (0, 0, u) RELATIVE PREVIOUS'];
l{end+1}=['  ROTATED (selene_epsilon_y' num ', selene_epsilon_x' num ', 0) RELATIVE PREVIOUS'];
l{end+1}=[''];
l{end+1}=['COMPONENT slit_before_guide11' num ' = Slit('];
l{end+1}=['    xmin = -selene_entry_x' num ', xmax = -selene_entry_x' num '/4.0,'];
l{end+1}=['    ymin = selene_entry_y' num '/4.0, ymax = selene_entry_y' num ')'];
l{end+1}=['  AT (0, 0, selene_distance' num '-0.003) RELATIVE arm_selene_' num ''];
l{end+1}=[''];
l{end+1}=['COMPONENT sele_guide_11' num ' = Elliptic_guide_gravity('];
l{end+1}=['    l = selene_length' num ', linxw = selene_distance' num ', linyh = selene_distance' num ','];
l{end+1}=['    loutxw = selene_distance' num ', loutyh = selene_distance' num ','];
l{end+1}=['    xwidth = 2*selene_entry_x' num ', yheight = 2*selene_entry_y' num ','];
l{end+1}=['    dimensionsAt = "entrance", mright = selene_coating_x' num ', mtop = selene_coating_y' num ','];
l{end+1}=['    mleft = 0, mbottom = 0,'];
if gravity == 1
l{end+1}=['    R0=R0' num ', Qc=Qc' num ', alpha=alpha' num ', W=W' num ',enableGravity=1)'];
else 
l{end+1}=['    R0=R0' num ', Qc=Qc' num ', alpha=alpha' num ', W=W' num ',enableGravity=0)'];
end
l{end+1}=['  AT (0,0,selene_distance' num ') RELATIVE arm_selene_' num ''];
l{end+1}=[''];
l{end+1}=['COMPONENT sele_guide_12' num ' = Elliptic_guide_gravity('];
l{end+1}=['    l = selene_length' num ', linxw = selene_distance' num ', linyh = selene_distance' num ','];
l{end+1}=['    loutxw = selene_distance' num ', loutyh = selene_distance' num ','];
l{end+1}=['    xwidth = 2*selene_entry_x' num ', yheight = 2*selene_entry_y' num ','];
l{end+1}=['    dimensionsAt = "entrance", mleft = selene_coating_x' num ', mbottom = selene_coating_y' num ','];
l{end+1}=['    mtop = 0, mright = 0,'];
if gravity == 1
l{end+1}=['    R0=R0' num ', Qc=Qc' num ', alpha=alpha' num ', W=W' num ',enableGravity=1)'];
else 
l{end+1}=['    R0=R0' num ', Qc=Qc' num ', alpha=alpha' num ', W=W' num ',enableGravity=0)'];
end
l{end+1}=['  AT (0,0,2*selene_c' num '+selene_distance' num ') RELATIVE arm_selene_' num ''];
l{end+1}=['  ROTATED (0,0,0) RELATIVE arm_selene_' num ''];
l{end+1}=[''];
l{end+1}=['COMPONENT EndOfelement_' num '= Arm()'];
l{end+1}=['  AT (0,0,4*selene_c' num ') RELATIVE arm_selene_' num ''];
l{end+1}=['  ROTATED (-selene_epsilon_y' num ',-selene_epsilon_x' num ',0) RELATIVE arm_selene_' num ''];
% remember to adjust for sample_dist
% the available sample space is: 2*selene_c-selene_distance-selene_length

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];

% Selene will use the same code for coating optimization as each mirror
% will need a constant m value.

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace_seg=[McStasStr.trace_seg '\n\n' tracestring];
clear l;

% Tese next functions write the part of the instrument file for this module
% concerning start, end and length parameters.

% code for calculating the focus points and smallaxis based on start and
% end determined by the minimalist concept and ifit optimizer

% Problem with the smallaxis. It can not be smaller than the average of the
% openings. factor from 1 to 5 ?

    linp=length(McStasStr.input);
    McStasStr.input{linp+1}=['selene_max_coat_x' num ];
    if max_m < 0
        McStasStr.inputvalue(linp+1)=6;
    else
        McStasStr.inputvalue(linp+1)=max_m;
    end
    
    McStasStr.input{linp+2}=['selene_max_coat_y' num ];
    if max_m < 0
        McStasStr.inputvalue(linp+2)=6;
    else
        McStasStr.inputvalue(linp+2)=max_m;
    end
    
    
    
    if ~isfield(McStasStr.input,'Selene_min_wavelength')
    linp=length(McStasStr.input);
    McStasStr.input{linp+1}='Selene_min_wavelength';
    McStasStr.inputvalue(linp+1)=demands.WaveLmin;
    end
    
    
    
    
    McStasStr.declare{end+1}=['divergence_x' num];
    McStasStr.declare{end+1}=['divergence_y' num];
    McStasStr.declare{end+1}=['total_length' num];
    McStasStr.declare{end+1}=['delta_theta_x' num];
    McStasStr.declare{end+1}=['delta_theta_y' num];
    McStasStr.declare{end+1}=['selene_xi' num];
    McStasStr.declare{end+1}=['selene_lambda_min' num];
    McStasStr.declare{end+1}=['selene_epsilon_x' num];
    McStasStr.declare{end+1}=['selene_epsilon_y' num];
    McStasStr.declare{end+1}=['selene_boa_x' num];
    McStasStr.declare{end+1}=['selene_boa_y' num];
    McStasStr.declare{end+1}=['selene_coating_x' num];
    McStasStr.declare{end+1}=['selene_coating_y' num];
    McStasStr.declare{end+1}=['selene_c' num];
    McStasStr.declare{end+1}=['selene_b_x' num];
    McStasStr.declare{end+1}=['selene_b_y' num];
    McStasStr.declare{end+1}=['selene_length' num];
    McStasStr.declare{end+1}=['selene_distance' num];
    McStasStr.declare{end+1}=['selene_entry_x' num];
    McStasStr.declare{end+1}=['selene_entry_y' num];
    if index == 1
        McStasStr.declare{end+1}=['original_sample_dist'];
    end
    if index == last
        McStasStr.declare{end+1}=['selene_mod_x'];
        McStasStr.declare{end+1}=['selene_mod_y'];
        McStasStr.declare{end+1}=['selene_max_mod_x'];
        McStasStr.declare{end+1}=['selene_max_mod_y'];
    end
        
    
    
    
    
    
  l{1} = '';
  l{end+1} = ['divergence_x' num ' = 2*var_divreq_x;'];
  l{end+1} = ['divergence_y' num ' = 2*var_divreq_y;'];
  l{end+1} = ['total_length' num ' = length' num ';']; 
  l{end+1} = [''];
  l{end+1} = ['delta_theta_x' num ' = divergence_x' num '*3.14159/180;'];
  l{end+1} = ['delta_theta_y' num ' = divergence_y' num '*3.14159/180;'];
  l{end+1} = ['selene_lambda_min' num ' = Selene_min_wavelength;'];
  
  l{end+1} = ['/* important selene input parameter. */'];  
    l{end+1} = ['/* the percentage of space with guide */'];  
    l{end+1} = ['selene_xi' num ' = 0.6;'];
  
  l{end+1} = ['selene_c' num '        = total_length' num ' / 4.0 ;'];
  l{end+1} = ['selene_length' num '   = 2.0 * selene_c' num ' * selene_xi' num ' ;'];
  l{end+1} = ['selene_distance' num ' = ( 1.0 - selene_xi' num ' ) * selene_c' num ' ;'];
  
  
  if globalinfo.selene_moderator_focus == 1 && index == last
      % Need to ensure that the guide starts at least 2 meters after the
      % guide
      
      l{end+1} = ['if (selene_distance1  < closest_element ) {'];    
      l{end+1} = ['selene_xi' num ' = - (closest_element - selene_c1)/selene_c1 ;'];  
      l{end+1} = ['}'];
      l{end+1} = [''];
  end    
  
  
  
  % why use such a low slope?
  %double my_alpha        =  2.74 ;    // slope of SM reflectivity
 
  % In this section delta theta must be in radians
  l{end+1} = ['/* derived quantities */'];
  l{end+1} = ['selene_epsilon_x' num ' = delta_theta_x' num ' / ( 2.0 * selene_xi' num ' ) ;'];
  l{end+1} = ['selene_epsilon_y' num '  = delta_theta_y' num ' / ( 2.0 * selene_xi' num ' ) ;'];
  l{end+1} = ['selene_boa_x' num '      = sqrt( pow(selene_epsilon_x' num ',2.0) - pow(delta_theta_x' num ',2.0)/4.0 ) ;'];
  l{end+1} = ['selene_boa_y' num '      = sqrt( pow(selene_epsilon_y' num ',2.0) - pow(delta_theta_y' num ',2.0)/4.0 ) ;'];
  % Original from Jochen, assumes m*qc is nicly reflecting, but in the
  % mcstas model that is not the case.
  %l{end+1} = ['selene_coating_x' num '  = 1.2 * 4.0*3.14159 * selene_epsilon_x' num ' / ( selene_lambda_min' num ' * 0.022 ) ;'];
  %l{end+1} = ['selene_coating_y' num '  = 1.2 * 4.0*3.14159 * selene_epsilon_y' num ' / ( selene_lambda_min' num ' * 0.022 ) ;'];
  % Model from Werner.
  l{end+1} = ['selene_coating_x' num '  = 1.17*1.2 * 4.0*3.14159 * selene_epsilon_x' num ' / ( selene_lambda_min' num ' * 0.022 ) ;'];
  l{end+1} = ['selene_coating_y' num '  = 1.17*1.2 * 4.0*3.14159 * selene_epsilon_y' num ' / ( selene_lambda_min' num ' * 0.022 ) ;'];
  l{end+1} = [''];
  l{end+1} = ['if (selene_coating_x' num ' > selene_max_coat_x' num ') selene_coating_x' num ' = selene_max_coat_x' num ';'];
  l{end+1} = ['if (selene_coating_y' num ' > selene_max_coat_y' num ') selene_coating_y' num ' = selene_max_coat_y' num ';'];
  l{end+1} = [''];
  l{end+1} = ['printf("Selene' num ' coating_x = %%lf \\n",selene_coating_x' num ');'];
  l{end+1} = ['printf("Selene' num ' coating_y = %%lf \\n",selene_coating_y' num ');'];
  l{end+1} = [''];
  l{end+1} = ['/* parameters for the guides */'];
  
  l{end+1} = ['selene_b_x' num '        = selene_c' num ' / sqrt( pow(selene_boa_x' num ',-2.0) - 1.0 ) ;'];
  l{end+1} = ['selene_b_y' num '        = selene_c' num ' / sqrt( pow(selene_boa_y' num ',-2.0) - 1.0 ) ;'];
  l{end+1} = ['selene_entry_x' num '    = selene_b_x' num ' * sqrt( 1.0 - pow(selene_xi' num ',2.0) ) ;'];
  l{end+1} = ['selene_entry_y' num '    = selene_b_y' num ' * sqrt( 1.0 - pow(selene_xi' num ',2.0) ) ;'];
  l{end+1} = [''];
  l{end+1} = ['/* rad2deg */'];
  l{end+1} = ['selene_epsilon_x' num ' *= 180.0/3.14159 ;'];
  l{end+1} = ['selene_epsilon_y' num ' *= 180.0/3.14159 ;'];
  l{end+1} = [''];
  if index == 1  % TEMP SOLUTION correct selene_xi dynamicaly to improve.
    current_sample_dist = McStasStr.inputvalue(globalinfo.sample_dist_index);
    l{end+1} = ['original_sample_dist = ' num2str(current_sample_dist) ';'];
    l{end+1} = ['if (original_sample_dist > selene_distance' num ') {'];
    l{end+1} = ['printf("Selene' num ' does not leave the requested sample space of %%lf m, but only %%lf m, reduce selene_xi.\\n",original_sample_dist,selene_distance1);'];
    l{end+1} = ['}'];
    l{end+1} = [''];
    
    McStasStr.inputvalue(globalinfo.sample_dist_index) = 0.0;
  end  
    % Handle moderator size to allow for better performance:
    % The lines below are in the moderator code in the trace section
    %l{end+1}='yheight = selene_mod_y, xwidth = selene_mod_x,';
    %l{end}=['dist = selene_source_focus_dist, focus_xw = selene_focus_width, focus_yh = selene_focus_height,'];
    
  if index == last && globalinfo.selene_moderator_focus == 0
  % Need to focus on a smaller part of the moderator for performance:
  % moderator_x_max = (selene_entry_x+startx)*selene_distance/(guide_start+selene_distance)
  % moderator_y_max = (selene_entry_y+starty)*selene_distance/(guide_start+selene_distance)
    l{end+1} = ['selene_max_mod_x = (selene_entry_x' num '*3/4+startx' num ')*(guide_start+selene_distance' num ')/selene_distance' num '-selene_entry_x' num '*3/4;'];
    l{end+1} = ['selene_max_mod_y = (selene_entry_y' num '*3/4+starty' num ')*(guide_start+selene_distance' num ')/selene_distance' num '-selene_entry_y' num '*3/4;'];
    l{end+1} = ['if (mod_x > selene_max_mod_x)'];
    l{end+1} = ['selene_mod_x = selene_max_mod_x;'];
    l{end+1} = ['else'];
    l{end+1} = ['selene_mod_x = mod_x;'];
    l{end+1} = [''];
    l{end+1} = ['if (mod_y > selene_max_mod_y)'];
    l{end+1} = ['selene_mod_y = selene_max_mod_y;'];
    l{end+1} = ['else'];
    l{end+1} = ['selene_mod_y = mod_y;'];
    l{end+1} = [''];
  end
  
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
    
    
    % END
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['endx' num];
    McStasStr.declare{ldec+2}=['endy' num];
    
    l{1}=['endx' num ' = startx' num ';'];
    l{2}=['endy' num ' = starty' num ';'];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    
    
    % START
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['startx' num];
    McStasStr.declare{ldec+2}=['starty' num];
    
    if index == 1
    l{1}=['startx' num ' = sizeX;'];
    l{2}=['starty' num ' = sizeY;'];
    else
    l{1}=['startx' num ' = startx' numM1 ';'];
    l{2}=['starty' num ' = starty' numM1 ';'];
    end
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
    

% length
McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,[0.1 0.95]);




end

