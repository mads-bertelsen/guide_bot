function [McStasStr] = K_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
%C_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
%Smodule adds text to McStas bot to make a straight guide
   % add the relevant component to the trace section
   % add the relevant input to the input list
   % add the relevant calculations to the initialize section
   % add the relevant variable names to the declare list
      
% default options   
for ops=1:length(options)    
    keyword='rot=';
    klength=length(keyword);
    % check for possible option rot
    if strcmp(options{ops}(1:klength),keyword)
       rot=str2num(options{ops}(klength+1:end));
       rot=abs(rot);
       overwrite=1;
    end
       
end
   
num=num2str(index);
numM1=num2str(index-1);


if (index==last || index ==1) 
    if (index==last)
        disp('ERROR, do not place K at the start of the input-string, change closest_element instead')
        disp('ERROR, the first K module will be ignored')
    end
    if (index==1)
        disp('ERROR, do not place K at the end of the input-string, change sample_dist instead')
        disp('ERROR, the last K module will be ignored')
    end
else % add a free space to the instrument file

% trace string   
if globalinfo.rotlogic(index)>0
    l{1}='// Horizontal kink module of the length in the following arm';
    l{2}=['COMPONENT CenterArm_' num '= Arm()'];
    if globalinfo.rotsign(index)<0
    l{3}=['AT (translation' num ',0,length' num ' ) RELATIVE PREVIOUS'];
    else
    l{3}=['AT (-translation' num ',0,length' num ' ) RELATIVE PREVIOUS'];
    end
    l{4}='';
    l{5}=['COMPONENT EndOfelement_' num '= Arm()'];
    l{6}='AT (0,0,0) RELATIVE PREVIOUS';
    if globalinfo.rotsign(index)<0
    l{7}=['ROTATED (0,rot' num ',0) RELATIVE PREVIOUS'];
    else
    l{7}=['ROTATED (0,-rot' num ',0) RELATIVE PREVIOUS'];
    end
else
    l{1}='// Vertical kink module of the length in the following arm';
    l{2}=['COMPONENT CenterArm_' num '= Arm()'];
    if globalinfo.rotsign(index)<0
    l{3}=['AT (0,translation' num ',length' num ' ) RELATIVE PREVIOUS'];
    else
    l{3}=['AT (0,-translation' num ',length' num ' ) RELATIVE PREVIOUS'];    
    end
    l{4}='';
    l{5}=['COMPONENT EndOfelement_' num '= Arm()'];
    l{6}='AT (0,0,0) RELATIVE PREVIOUS';
    if globalinfo.rotsign(index)>0
    l{7}=['ROTATED (rot' num ',0,0) RELATIVE PREVIOUS'];
    else
    l{7}=['ROTATED (-rot' num ',0,0) RELATIVE PREVIOUS'];
    end
end

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];
McStasStr.trace_seg=[McStasStr.trace_seg '\n\n' tracestring];
%McStasStr.trace=[tracestring '\n\n' McStasStr.trace ];
clear l;

%     l{1}='';
%     l{end+1} =['startxpoint[' num '][1][1] = endPoint' num ' - length' num ';'];
%     l{end+1} =['startxpoint[' num '][2][1] = -0.5*startx' num ';'];
%     l{end+1} =['startxpoint[' num '][1][2] = endPoint' num ' - length' num ';'];
%     l{end+1} =['startxpoint[' num '][2][2] = 0.5*startx' num ';'];
% 
%     l{end+1} =['startypoint[' num '][1][1] = endPoint' num ' - length' num ';'];
%     l{end+1} =['startypoint[' num '][2][1] = -0.5*starty' num ';'];
%     l{end+1} =['startypoint[' num '][1][2] = endPoint' num ' - length' num ';'];
%     l{end+1} =['startypoint[' num '][2][2] = 0.5*starty' num ';'];
% 
%     l{end+1}='';
%     
%     if globalinfo.rotlogic(index)>0
%     l{end+1} =['    kinkendpointx[1] = endPoint' num ';'];
%     if globalinfo.rotsign(index)>0
%     l{end+1} =['    kinkendpointx[2] = translation' num ';'];
%     else
%     l{end+1} =['    kinkendpointx[2] = -translation' num ';'];
%     end
%     %l{end+1} =['    kinkendpointy[1] = endPoint' num ';'];
%     %l{end+1} =['    kinkendpointy[2] = 0;'];
%     else
%     l{end+1} =['    kinkendpointy[1] = endPoint' num ';'];
%     if globalinfo.rotsign(index)>0
%     l{end+1} =['    kinkendpointy[2] = translation ' num ';'];
%     else
%     l{end+1} =['    kinkendpointy[2] = -translation' num ';'];
%     end
%     %l{end+1} =['    kinkendpointx[1] = endPoint' num ';'];
%     %l{end+1} =['    kinkendpointx[2] = 0;'];
%     end
%     
%     initstring='';
%     for i=1:length(l)
%         initstring=[initstring l{i} '\n'];
%     end
%     clear l;
%     McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    


% start
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['startx' num];
    McStasStr.declare{ldec+2}=['starty' num];

    % Here a certain solution is chossen. Same var_div_req, but larger end
    % of the former guide element which enlarges the phase-space
    % requirement. But it could be done with larger var_div_req and then
    % not iluminating the entire guide start of the latter element.
    
    if globalinfo.rotlogic(index)>0
    l{1} = ['startx' num ' = endx' num ' + length' num '*( tan((var_divreq_x+rot' num ')*DEG2RAD)+tan((fabs(var_divreq_x-rot' num '))*DEG2RAD) );'];
    l{2} = ['starty' num ' = endy' num ' + 2*length' num '*tan(var_divreq_y*DEG2RAD);'];   
    else
    l{1} = ['startx' num ' = endx' num ' + 2*length' num '*tan(var_divreq_x*DEG2RAD);'];   
    l{2} = ['starty' num ' = endy' num ' + length' num '*( tan((var_divreq_y+rot' num ')*DEG2RAD)+tan((fabs(var_divreq_y-rot' num '))*DEG2RAD) );'];
    end
        
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
    

% Calculate the rot to avoid line of sight (initialize in instr file)

% McStasStr.declare{end+1}=['translation' num];
% if (globalinfo.kinkoverwrite==0)
%     McStasStr.declare{end+1}=['rot' num];
%     if globalinfo.rotlogic(index)>0
%         l{1} = '// This calculation of rot is an approximation and only true if there is only one kink!s';
%         l{end+1} = ['rot' num '= RAD2DEG*atan(0.5*(startx' numM1 ' + 0.06)*(Mod_sample - guide_start - sample_dist)/(endPoint' numM1 ' - length' numM1 ')/( Mod_sample - guide_start - sample_dist - endPoint' numM1 ' + length' numM1 '));'];% TEMP CODE FOR CALCULATING ROT
%         l{end+1} =['translation' num '= 0.5 * length' num '*(tan((var_divreq_x+rot' num ')*DEG2RAD)-tan((fabs(var_divreq_x-rot' num '))*DEG2RAD) );'];
%     else
%         l{1} = '// This calculation of rot is an approximation and only true if there is only one kink!s';
%         l{end+1} = ['rot' num '= RAD2DEG*atan(0.5*(starty' numM1 ' + 0.06)*(Mod_sample - guide_start - sample_dist)/(endPoint' numM1 ' - length' numM1 ')/( Mod_sample - guide_start - sample_dist - endPoint' numM1 ' + length' numM1 '));'];% TEMP CODE FOR CALCULATING ROT
%         l{end+1} =['translation' num '= 0.5 * length' num '*(tan((var_divreq_y+rot' num ')*DEG2RAD)-tan((fabs(var_divreq_y-rot' num '))*DEG2RAD) );'];
%     end
% else
%     linp=length(McStasStr.input);
%     McStasStr.input{linp+1}=['rot' num];
%     McStasStr.inputvalue(linp+1)=rot;
%     if globalinfo.rotlogic(index)>0
%         l{1} = '// Calculation of kink angle overwritten, manual value added as input';    
%         l{end+1} =['translation' num '= 0.5 * length' num '*(tan((var_divreq_x+rot' num ')*DEG2RAD)-tan((fabs(var_divreq_x-rot' num '))*DEG2RAD) );'];
%     else
%         l{1} = '// Calculation of kink angle overwritten, manual value added as input';    
%         l{end+1} =['translation' num '= 0.5 * length' num '*(tan((var_divreq_y+rot' num ')*DEG2RAD)-tan((fabs(var_divreq_y-rot' num '))*DEG2RAD) );'];
%     end
% end

McStasStr.declare{end+1}=['translation' num];
if (globalinfo.kinkoverwrite(index)==0)
    
    current_los_logic = false(1,last);
    current_los_logic(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section})) = 1;
    current_raytracers_logic = current_los_logic.*globalinfo.raytracers;
    current_raytracers_index = find(current_raytracers_logic);
    master_rot_index = current_raytracers_index(1);
    
    if index == master_rot_index
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
        l{end+1}=['rot' num ' = rot' num ' + 0.002;'];
        l{end+1}=[''];
        if globalinfo.rotlogic(index)>0
        l{end+1}=['var_divreq_x = var_divreq_x_protected+rot' num ';'];
        l{end+1}=['var_divreq_y = var_divreq_y_protected;'];
        else
        l{end+1}=['var_divreq_x = var_divreq_x_protected;'];
        l{end+1}=['var_divreq_y = var_divreq_y_protected+rot' num ';'];    
        end
        
        l{end+1}=[''];
    else
    % This los breaker is a slave    
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
        if globalinfo.rotlogic(index)>0
        l{end+1}=['var_divreq_x = var_divreq_x+rot' num ';'];
        else
        l{end+1}=['var_divreq_y = var_divreq_y+rot' num ';'];    
        end
        l{end+1}=[''];
    end
else
    linp=length(McStasStr.input);
    McStasStr.input{linp+1}=['rot' num];
    McStasStr.inputvalue(linp+1)=rot;
    l{1}='';
end
    

if globalinfo.rotlogic(index)>0
        %l{end+1} = '// Calculation of kink angle overwritten, manual value added as input';    
        l{end+1} =['translation' num '= 0.5 * length' num '*(tan((var_divreq_x+rot' num ')*DEG2RAD)-tan((fabs(var_divreq_x-rot' num '))*DEG2RAD) );'];
else
        %l{end+1} = '// Calculation of kink angle overwritten, manual value added as input';    
        l{end+1} =['translation' num '= 0.5 * length' num '*(tan((var_divreq_y+rot' num ')*DEG2RAD)-tan((fabs(var_divreq_y-rot' num '))*DEG2RAD) );'];
end




initstring='';
for i=1:length(l)
    initstring=[initstring l{i} '\n'];
end
clear l;
McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    
    
    
%     
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
%    
%     McStasStr.declare{end+1}=['endx' num];
%     McStasStr.declare{end+1}=['endy' num];
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
% 
% 
% % length
%     % fraction of remaning length given as optimized input
%     linp=length(McStasStr.input);
%     McStasStr.input{linp+1}=['lengthfrac' num];    
%     McStasStr.optimize(linp+1)=1;
%     McStasStr.optimvals.min(linp+1)=0.01;
%     McStasStr.optimvals.max(linp+1)=0.30; % hard to determine
%     McStasStr.optimvals.guess(linp+1)=1./(1+last-index);
%     
%     ldec=length(McStasStr.declare);
%     McStasStr.declare{ldec+1}=['length' num];
%     
%     %strings=['length' num ' = (Mod_sample - sample_dist - guide_start'];
%     %for i=1:index-1
%     %    strings = [strings ' - length' num2str(i)];
%     %end
%     %strings = [ strings '- u)*lengthfrac' num ';'];
%     
%     strings=['length' num ' = (endPoint' num ' - guide_start '];
%     strings = [ strings ')*lengthfrac' num ';'];
%     
%     McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
%     
% 
% % EndPoint
%     ldec=length(McStasStr.declare);
%     McStasStr.declare{ldec+1}=['endPoint' num];
% 
%     % Udregning af hvor endpoint num skal være, baseret på længderne af
%     % alle andre elementer.
%     if (index==1)
%         strings=['endPoint' num ' = Mod_sample - sample_dist;'];    
%     else 
%         strings=['endPoint' num ' = endPoint' numM1 ' - length' numM1 ' - u;'];
%     end    
%     McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
    
% end
McStasStr=guide_writer_end(McStasStr,index,last);

% length 
McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,[0.01 0.25 0.02]);

end % end of overall logic, runs if no error.




end

